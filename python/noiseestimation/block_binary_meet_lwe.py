"""Combinatorial attack costs on anchored block-binary secrets (ordinary Python).

TUTTI (ePrint 2026/2848) switches to a block-binary secret before evaluating
the decryption equation: N coefficients split into h blocks of B = N/h, one
coefficient equal to 1 per block, and s_0 = 1 public. A subring key
s(X) = s'(X^R) is an instance of degree N/R (their Proposition 1), so every
function here takes that reduced degree.

Implemented, following their Section 3.2 and Appendix D:

* Odlyzko's meet-in-the-middle over whole blocks (Algorithm 9);
* May's Meet-LWE with Rep-1 representations for a one-signed secret and a
  bit-granular filter, with the error guessed from a bounded alphabet
  (Algorithms 10 and 11, depths 2 to 4);
* the depth-2 guessing-free tree of Hhan-Kim-Lee-Son (ePrint 2022/1473) with
  the cost model of Bi et al. (ePrint 2026/1060, Section 3.3), using the box
  and bucket grids of the authors' public estimator.

Costs are log2 list elements, without polynomial factors, as in May's
tables; `*_ops` adds log2(log2 L) for hashing an element of the largest list.
The RLWE sample is read as uniform LWE equations: the ring is not used and
hybrids of these lists with lattice reduction are not estimated. A cost here
is a model-dependent screen, not a security proof.
"""

from dataclasses import dataclass
import math

NEG_INF = float("-inf")
_LN2 = math.log(2.0)


def lg(x):
    return math.log2(x) if x > 0 else NEG_INF


def log2_factorial(n):
    return 0.0 if n <= 1 else math.lgamma(n + 1) / _LN2


def log2_multinomial(n, *parts):
    """log2 of n! / (p_1! ... p_k! (n - sum p)!), or -inf when it is zero."""
    used = sum(parts)
    if any(p < 0 for p in parts) or used > n:
        return NEG_INF
    return log2_factorial(n) - sum(map(log2_factorial, parts)) - log2_factorial(n - used)


@dataclass(frozen=True)
class AnchoredBlockBinary:
    """The instance an attacker faces at degree `degree` = N/R with h blocks."""
    degree: int
    blocks: int

    def __post_init__(self):
        if self.degree < 2 or self.degree & (self.degree - 1):
            raise ValueError("degree must be a power of two")
        if self.blocks < 2 or self.blocks & (self.blocks - 1) or self.blocks > self.degree:
            raise ValueError("blocks must be a power of two between 2 and the degree")

    @property
    def block_length(self):
        return self.degree // self.blocks

    @property
    def dimension(self):
        return self.degree - self.block_length

    @property
    def weight(self):
        return self.blocks - 1

    @property
    def log2_search_space(self):
        return self.weight * math.log2(self.block_length)


def odlyzko_bits(instance):
    """Larger list of the split over whole blocks: ceil(h'/2) * log2 B."""
    if instance.block_length == 1:
        return 0.0
    return math.ceil(instance.weight / 2) * math.log2(instance.block_length)


def rep1_counts(dimension, weight, epsilons):
    """Algorithm 10: (log2 S_j, log2 R_j) per level, or None when rejected."""
    if not epsilons or any(e < 0 for e in epsilons):
        raise ValueError("at least one nonnegative level parameter is required")
    first = epsilons[0]
    if 2 * first > dimension - weight:
        return None
    ones = math.ceil(weight / 2)
    plus, minus = ones + first, first
    counts = [(log2_multinomial(dimension, plus, minus),
               log2_multinomial(weight, ones) + log2_multinomial(dimension - weight, first, first))]
    for eps in epsilons[1:]:
        zeros = dimension - plus - minus
        if 2 * eps > zeros:
            return None
        representations = (log2_multinomial(plus, math.ceil(plus / 2)) +
                           log2_multinomial(minus, math.ceil(minus / 2)) +
                           log2_multinomial(zeros, eps, eps))
        plus, minus = math.ceil(plus / 2) + eps, math.ceil(minus / 2) + eps
        counts.append((log2_multinomial(dimension, plus, minus), representations))
    return counts


def guess_alphabet_size(sigma):
    """Symmetric alphabet covering about +-2.5 standard deviations."""
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be positive and finite")
    return 2 * math.ceil(2.5 * sigma) + 1


def tree_bits(dimension, log_q, counts, alphabet):
    """Algorithm 11: (T_lists, T_guess, guessed coordinates), or None."""
    representations = [r for _, r in counts]
    if any(r <= 0 for r in representations) or any(
            a < b for a, b in zip(representations, representations[1:])):
        return None
    lists = [s - r for s, r in counts]
    if any(size < 0 for size in lists):
        return None
    lists.append(counts[-1][0] / 2)
    representations.append(0.0)
    cost = lists[-1]
    for j in range(len(counts) - 1, -1, -1):
        cost = max(cost, lists[j],
                   2 * lists[j + 1] - (representations[j] - representations[j + 1]))
    coordinates = math.ceil(representations[0] / log_q)
    cost = max(cost, 2 * lists[0] - (dimension - coordinates))
    return cost, coordinates / 2 * math.log2(alphabet), coordinates


def meet_lwe_guessing(instance, log_q, sigma=3.2, max_depth=4):
    """Cheapest Rep-1 tree with error guessing over May's depths 2..max_depth."""
    dimension, weight = instance.dimension, instance.weight
    alphabet = guess_alphabet_size(sigma)
    limits = (40, 24, 12)
    best = None

    def visit(prefix):
        nonlocal best
        counts = rep1_counts(dimension, weight, prefix)
        if counts is not None:
            cost = tree_bits(dimension, log_q, counts, alphabet)
            if cost is not None and (best is None or sum(cost[:2]) < best["bits"]):
                best = {"bits": cost[0] + cost[1], "list_bits": cost[0],
                        "guess_bits": cost[1], "guessed_coordinates": cost[2],
                        "depth": len(prefix) + 1, "epsilons": list(prefix)}
        if len(prefix) + 1 < max_depth:
            bound = min(prefix[-1], limits[len(prefix)])
            for eps in range(bound + 1):
                visit(prefix + (eps,))

    for first in range(limits[0] + 1):
        visit((first,))
    return best


_BOX_GRID = sorted(set(
    [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 256,
     320, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6144, 8192, 12288, 16384, 32768, 65536]
    + [2**k for k in range(17, 200, 2)]))
_BUCKET_GRID = [x / 8 for x in range(4, 65)] + [16, 32, 64, 128]


def _gaussian(sigma):
    span = max(60, math.ceil(12 * sigma))
    weights = {e: math.exp(-e * e / (2 * sigma * sigma)) for e in range(-span, span + 1)}
    total = sum(weights.values())
    return {e: w / total for e, w in weights.items()}


def _guessing_free_top(dimension, log_q, log_s, log_r, law, repetitions):
    leaf = log_s / 2
    boxes = [l for l in _BOX_GRID if 2 * l + 1 < 2 ** (log_q - 1)]

    def survive(l):
        return sum(p * max(0.0, 1 - abs(e) / (2 * l + 1)) for e, p in law.items())

    def same_bucket(l, width):
        return 1 - l * (l + 1) / ((2 * l + 1) * width) if width >= l else width / (2 * l + 1)

    consumed = {l: (log_q - lg(2 * l + 1) - lg(survive(l)), lg(2 * l + 1)) for l in boxes}
    hashes = {l: [None if l * c >= 2**log_q else
                  (-lg(same_bucket(l, l * c)), lg(l * c) - log_q) for c in _BUCKET_GRID]
              for l in boxes}
    best = None
    for l in boxes:
        bits = consumed[l][0]
        if bits <= 0:
            continue
        full = int(math.floor(log_r / bits))
        rest = log_r - full * bits
        for extra in [None] + [x for x in boxes if 0 < consumed[x][0] <= rest + 1e-9]:
            coordinates = full + (extra is not None)
            if coordinates == 0:
                continue
            top = log_s - full * (log_q - consumed[l][1])
            if extra is not None:
                top -= log_q - consumed[extra][1]
            merge = max(top, 2 * top - (dimension - coordinates))
            for index, factor in enumerate(_BUCKET_GRID):
                a = hashes[l][index] if full else (0.0, 0.0)
                b = hashes[extra][index] if extra is not None else (0.0, 0.0)
                if a is None or b is None:
                    continue
                repeats = full * a[0] + b[0] + repetitions
                search = max(leaf, 2 * leaf + full * a[1] + b[1]) + repeats
                total = max(search, merge, top, leaf)
                if best is None or total < best["bits"]:
                    best = {"bits": total, "top_list_bits": top, "leaf_list_bits": leaf,
                            "repetition_bits": repeats, "box_coordinates": coordinates,
                            "box_half_width": l, "extra_box_half_width": extra,
                            "bucket_factor": factor}
    return best


def meet_lwe_guessing_free(instance, log_q, sigma=3.2, ten_repetitions=False, max_epsilon=30):
    """Depth-2 guessing-free tree; an upper bound where deeper trees may win."""
    law = _gaussian(sigma)
    repetitions = math.log2(10) if ten_repetitions else 0.0
    best = None
    for eps in range(max_epsilon):
        counts = rep1_counts(instance.dimension, instance.weight, (eps,))
        if counts is None:
            continue
        top = _guessing_free_top(instance.dimension, log_q, counts[0][0], counts[0][1],
                                 law, repetitions)
        if top is not None and (best is None or top["bits"] < best["bits"]):
            best = dict(top, epsilon=eps)
    if best is None:
        return None
    best["ops_bits"] = best["bits"] + lg(max(best["top_list_bits"], best["leaf_list_bits"]))
    return best


def combinatorial_report(degree, blocks, log_q, sigma=3.2):
    """All combinatorial screens for one (N/R, h, log q0) row."""
    instance = AnchoredBlockBinary(degree, blocks)
    odlyzko = odlyzko_bits(instance)
    guessing = meet_lwe_guessing(instance, log_q, sigma)
    free = meet_lwe_guessing_free(instance, log_q, sigma)
    conservative = meet_lwe_guessing_free(instance, log_q, sigma, ten_repetitions=True)
    search = instance.log2_search_space
    costs = {"odlyzko_ops": odlyzko + lg(odlyzko), "guessing_free_ops": free["ops_bits"],
             "meet_lwe_guessing_lists": guessing["bits"]}
    return {
        "scope": __doc__.splitlines()[0],
        "degree": degree, "blocks": blocks, "block_length": instance.block_length,
        "anchored_dimension": instance.dimension, "anchored_weight": instance.weight,
        "log_q": log_q, "sigma": sigma, "guess_alphabet": guess_alphabet_size(sigma),
        "log2_search_space": search,
        "odlyzko_list_bits": odlyzko, "odlyzko_ops_bits": costs["odlyzko_ops"],
        "meet_lwe_guessing": guessing,
        "guessing_free": free,
        "guessing_free_ten_repetitions": conservative,
        "closed_exponents_not_attacks": {"S^0.25": 0.25 * search, "S^0.28": 0.28 * search,
                                         "S^0.4": 0.4 * search},
        "cheapest_attack": min(costs, key=costs.get),
        "cheapest_bits": min(costs.values()),
        "passes_128": min(costs.values()) >= 128,
        "excluded": ["ring structure", "lattice hybrids of the representation lists",
                     "van Oorschot-Wiener variant", "polynomial factors", "memory cost"],
    }
