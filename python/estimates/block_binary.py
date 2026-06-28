#!/bin/python3
"""Block-binary LWE security estimates.

This script keeps the vendored lattice-estimator unchanged.  It wraps the
public estimator APIs and restricts hybrid guessing to whole block-binary
blocks, as used in Lee et al., "Faster TFHE Bootstrapping with Block Binary
Keys" (ePrint 2023/958).
"""

import argparse
import importlib
import math
import os
import sys
from functools import partial


sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

estimator = importlib.import_module(".estimator", "lattice-estimator")
lwe_primal = importlib.import_module(".estimator.lwe_primal", "lattice-estimator")
lwe_dual = importlib.import_module(".estimator.lwe_dual", "lattice-estimator")
lwe_guess = importlib.import_module(".estimator.lwe_guess", "lattice-estimator")
conf = importlib.import_module(".estimator.conf", "lattice-estimator")
util = importlib.import_module(".estimator.util", "lattice-estimator")
cost_module = importlib.import_module(".estimator.cost", "lattice-estimator")
errors = importlib.import_module(".estimator.errors", "lattice-estimator")

try:
    from sage.all import RR, ceil, exp, log, oo, pi, sqrt
except ImportError as exc:  # pragma: no cover - this script is for sage-python.
    raise SystemExit("Run this script with sage -python.") from exc


TFHE636_ALPHA = 0.000_092_511_997_467_675_6
PAPER_PARAMS = (
    # ell, n, N, paper Esti, paper Modified Est, paper MitM
    (2, 630, 1024, 128.8, 128.8, 139.793),
    (3, 687, 1024, 128.3, 126.7, 128.240),
    (4, 788, 1024, 128.6, 127.4, 128.078),
    (5, 885, 1024, 127.8, 126.2, 128.111),
    (6, 978, 1024, 127.2, 125.4, 128.128),
)


class BlockBinary:
    """Distribution B_{ell,k}, represented locally for estimator wrappers.

    The estimator only needs aggregate moments and support sizes here.  We mark
    the distribution as non-sparse to avoid lattice-estimator's SparseTernary-
    specific split logic; this wrapper handles block-aligned splits itself.
    """

    def __init__(self, ell, n=None):
        self.ell = int(ell)
        self.n = None if n is None else int(n)
        if self.ell < 1:
            raise ValueError("ell must be positive")
        if self.n is not None and self.n % self.ell:
            raise ValueError(f"n={self.n} is not divisible by ell={self.ell}")

        p_one = RR(1) / RR(self.ell + 1)
        self.mean = p_one
        self.stddev = sqrt(p_one * (1 - p_one))
        self.bounds = (0, 1)
        self._density = p_one
        self.is_Gaussian_like = False

    def resize(self, new_n):
        return BlockBinary(self.ell, new_n)

    def __len__(self):
        if self.n is None:
            raise ValueError("Distribution has no length.")
        return self.n

    @property
    def hamming_weight(self):
        return int(round(float(len(self) * self._density)))

    @property
    def is_bounded(self):
        return True

    @property
    def is_sparse(self):
        return False

    def support_size(self, fraction=1.0):
        if self.n is None:
            raise ValueError("Distribution has no length.")
        blocks = ceil(RR(self.n) / RR(self.ell))
        return ceil(RR(fraction) * RR(self.ell + 1) ** blocks)

    def split_probability(self, new_n):
        if new_n % self.ell:
            return RR(0)
        return RR(1)

    def __lt__(self, other):
        try:
            return self.stddev < other.stddev
        except AttributeError:
            return self.stddev < other

    def __le__(self, other):
        try:
            return self.stddev <= other.stddev
        except AttributeError:
            return self.stddev <= other

    def __hash__(self):
        return hash(("BlockBinary", self.ell, self.n))

    def __repr__(self):
        if self.n is None:
            return f"BlockBinary(ell={self.ell})"
        return f"BlockBinary(ell={self.ell}, n={self.n})"

    __str__ = __repr__


def log2_rr(x):
    if x == oo:
        return math.inf
    return float(log(x, 2))


def cost_bits(cost):
    return log2_rr(cost["rop"])


def is_finite_cost(cost):
    return cost is not None and cost.get("rop", oo) < oo


def mitm_bits(ell, n):
    blocks = RR(n) / RR(ell)
    return float(RR("0.28") * blocks * log(ell + 1, 2))


def paper_profile():
    alpha = RR(2) ** -15
    beta = RR(2) ** -25
    return {
        "name": "paper",
        "q": 2 ** 30,
        "err_stddev": alpha * 2 ** 30,
        "baseline_secret": estimator.ND.Binary,
        "rlwe_q": 2 ** 50,
        "rlwe_err_stddev": beta * 2 ** 50,
        "N": 1024,
        "note": "TLWE alpha=2^-15 scaled to q=2^30; TRLWE beta=2^-25 scaled to q=2^50",
    }


def tfhe636_profile():
    q = 2 ** 16
    return {
        "name": "tfhe636",
        "q": q,
        "err_stddev": RR(str(TFHE636_ALPHA)) * q,
        "baseline_secret": estimator.ND.Binary,
        "rlwe_q": 2 ** 50,
        "rlwe_err_stddev": RR(2) ** 25,
        "N": 1024,
        "note": "TLWE q=2^16, error stddev=0.0000925119974676756*q, baseline Xs=Binary",
    }


def params_for(n, q, err_stddev, xs, tag):
    return estimator.lwe_parameters.LWEParameters(
        n=int(n),
        q=int(q),
        Xs=xs,
        Xe=estimator.ND.DiscreteGaussian(stddev=err_stddev),
        tag=tag,
    )


def binary_params(n, q, err_stddev, tag, xs=estimator.ND.Binary):
    # This is the paper's "Esti" baseline: run the unmodified estimator as if
    # the secret were ordinary dense binary, without exploiting block structure.
    return params_for(n, q, err_stddev, xs, tag)


def block_binary_params(ell, n, q, err_stddev, tag):
    return params_for(n, q, err_stddev, BlockBinary(ell), tag)


def min_cost(costs):
    finite = [c for c in costs if is_finite_cost(c)]
    if not finite:
        return cost_module.Cost(rop=oo, tag="none")
    return min(finite, key=lambda c: c["rop"])


def with_tag(cost, tag):
    if cost is not None:
        cost["tag"] = tag
    return cost


def estimate_original(params, red_cost_model, quiet=True):
    deny = ("arora-gb",)
    return estimator.LWE.estimate(
        params,
        red_cost_model=red_cost_model,
        deny_list=deny,
        quiet=quiet,
        catch_exceptions=True,
    )


def estimate_original_bits(params, red_cost_model, quiet=True):
    costs = estimate_original(params, red_cost_model, quiet=quiet)
    return min((cost_bits(c) for c in costs.values() if is_finite_cost(c)), default=math.inf)


def optimize_block_primal_hybrid(params, ell, red_cost_model, red_shape_model, mitm, babai):
    """Primal hybrid with block-aligned guesses and explicit search space."""

    params = params.normalize()
    m = params.m + params.n if params.Xs <= params.Xe else params.m
    baseline = estimator.LWE.primal_usvp(
        params,
        red_cost_model=red_cost_model,
        red_shape_model=red_shape_model,
        optimize_d=False,
        log_level=0,
    )
    max_beta = conf.max_beta if baseline["rop"] == oo else min(conf.max_beta, int(baseline["beta"]))
    precision = 4
    best = cost_module.Cost(rop=oo)

    for guessed_blocks in range(0, params.n // ell + 1):
        zeta = guessed_blocks * ell
        search_space = RR(ell + 1) ** guessed_blocks
        hit_probability = RR(1)
        f = partial(
            lwe_primal.PrimalHybrid.cost,
            params=params,
            zeta=zeta,
            babai=babai,
            mitm=mitm,
            m=m,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            search_space=search_space,
            hit_probability=hit_probability,
            log_level=0,
        )
        with util.local_minimum(40, max_beta + precision, precision=precision, log_level=0) as it:
            for beta in it:
                it.update(f(beta=beta))
            for beta in it.neighborhood:
                it.update(f(beta=beta))
            candidate = it.y

        if is_finite_cost(candidate) and candidate["rop"] < best["rop"]:
            best = candidate
            best["block_guesses"] = guessed_blocks
        # Once the whole secret is enumerated, larger zeta is impossible.
        if zeta == params.n:
            break

    tag = "bdd_mitm_block_hybrid" if mitm else "bdd_block_hybrid"
    return with_tag(best, tag)


def block_exhaustive_solver(params, success_probability=0.99):
    params = params.normalize()
    probability = sqrt(success_probability)
    try:
        size = params.Xs.support_size(probability)
    except NotImplementedError:
        return cost_module.Cost(rop=oo, mem=oo, m=1)

    sigma = params.Xe.stddev / params.q
    m_required = RR(8 * exp(4 * pi * pi * sigma * sigma)) * (
        log(size) - log(log(1 / probability))
    )
    m_required = RR(m_required)
    if params.m < m_required:
        raise errors.InsufficientSamplesError(
            f"Exhaustive search: Need {m_required} samples but only {params.m} available."
        )
    rop = 2 * size * m_required
    return cost_module.Cost(rop=rop, mem=rop / 2, m=m_required)


block_exhaustive_solver.__name__ = "block_exhaustive_solver"


def optimize_block_dual_hybrid(params, ell, red_cost_model):
    """Dual hybrid where the solved part is an integer number of blocks."""

    params = params.normalize()
    max_blocks = params.n // ell

    def estimate_for_blocks(solved_blocks):
        zeta = solved_blocks * ell
        try:
            return lwe_dual.DH.optimize_blocksize(
                solver=block_exhaustive_solver,
                params=params,
                zeta=zeta,
                h1=0,
                red_cost_model=red_cost_model,
                log_level=20,
                opt_step=8,
            )
        except Exception:
            return cost_module.Cost(rop=oo)

    with util.local_minimum(
        0,
        max_blocks + 1,
        precision=4,
        suppress_bounds_warning=True,
        log_level=20,
    ) as it:
        for solved_blocks in it:
            it.update(estimate_for_blocks(solved_blocks))
        for solved_blocks in it.neighborhood:
            it.update(estimate_for_blocks(solved_blocks))
        best = it.y

    if is_finite_cost(best):
        best["solved_blocks"] = it.x

    return with_tag(best, "dual_block_hybrid")


def optimize_modified_dual(params, ell, red_cost_model, block_n=None, tag="modified_dual"):
    """Modified dual attack from the block-binary key paper.

    Guess t whole blocks, paying (ell+1)^t, and run the usual dual attack on
    the remaining dimension with expected Hamming weight (n - t*ell)/(ell+1).
    This is the model stated in the public FHE.org 2023 presentation for the
    paper's modified estimator.
    """

    params = params.normalize()
    block_n = params.n if block_n is None else int(block_n)
    if block_n % ell:
        raise ValueError(f"block_n={block_n} is not divisible by ell={ell}")
    if block_n > params.n:
        raise ValueError(f"block_n={block_n} exceeds dimension {params.n}")
    binary_tail = params.n - block_n
    max_blocks = block_n // ell

    def estimate_for_blocks(guessed_blocks):
        remaining_n = params.n - guessed_blocks * ell
        if remaining_n <= 0:
            return cost_module.Cost(rop=oo)

        remaining_block_n = block_n - guessed_blocks * ell
        hw = int(round(remaining_block_n / (ell + 1) + binary_tail / 2))
        reduced = params.updated(
            n=remaining_n,
            Xs=estimator.ND.SparseBinary(hw),
            tag=f"{params.tag}-{tag}-t{guessed_blocks}",
        )
        try:
            cost = estimator.LWE.dual(reduced, red_cost_model=red_cost_model)
        except Exception:
            return cost_module.Cost(rop=oo)

        search_space = RR(ell + 1) ** guessed_blocks
        cost = cost.repeat(search_space)
        cost["guessed_blocks"] = guessed_blocks
        cost["zeta"] = guessed_blocks * ell
        cost["|S|"] = search_space
        return with_tag(cost, tag)

    with util.local_minimum(
        0,
        max_blocks,
        precision=4,
        suppress_bounds_warning=True,
        log_level=20,
    ) as it:
        for guessed_blocks in it:
            it.update(estimate_for_blocks(guessed_blocks))
        for guessed_blocks in it.neighborhood:
            it.update(estimate_for_blocks(guessed_blocks))
        best = it.y

    return with_tag(best, tag)


def estimate_modified(
    params,
    ell,
    red_cost_model,
    red_shape_model,
    block_n=None,
    include_primal_block=False,
    include_experimental_hybrid=False,
    lwe_only=False,
):
    costs = []
    for name, attack in (
        ("usvp", estimator.LWE.primal_usvp),
        ("bdd", estimator.LWE.primal_bdd),
        ("dual", estimator.LWE.dual),
    ):
        try:
            cost = attack(
                params,
                red_cost_model=red_cost_model,
                red_shape_model=red_shape_model,
            ) if name in ("usvp", "bdd") else attack(params, red_cost_model=red_cost_model)
            costs.append(with_tag(cost, name))
        except Exception:
            pass

    costs.append(
        optimize_modified_dual(
            params,
            ell,
            red_cost_model=red_cost_model,
            block_n=block_n,
            tag="modified_dual",
        )
    )

    if include_primal_block:
        costs.append(
            optimize_block_primal_hybrid(
                params,
                ell,
                red_cost_model=red_cost_model,
                red_shape_model=red_shape_model,
                mitm=False,
                babai=False,
            )
        )
        costs.append(
            optimize_block_primal_hybrid(
                params,
                ell,
                red_cost_model=red_cost_model,
                red_shape_model=red_shape_model,
                mitm=True,
                babai=True,
            )
        )
    if include_experimental_hybrid:
        costs.append(optimize_block_dual_hybrid(params, ell, red_cost_model=red_cost_model))
    return costs, min_cost(costs)


def estimate_row(
    ell,
    n,
    profile,
    red_cost_model,
    red_shape_model,
    quiet=True,
    include_primal_block=False,
    include_experimental_hybrid=False,
    lwe_only=False,
):
    tlwe_baseline_params = binary_params(
        n,
        profile["q"],
        profile["err_stddev"],
        f"binary-baseline-ell{ell}-n{n}",
        profile.get("baseline_secret", estimator.ND.Binary),
    )
    tlwe_block_params = block_binary_params(
        ell, n, profile["q"], profile["err_stddev"], f"block-binary-ell{ell}-n{n}"
    )
    N = int(profile.get("N", 1024))
    tlwe_original = estimate_original_bits(tlwe_baseline_params, red_cost_model, quiet=quiet)
    tlwe_modified_costs, tlwe_modified_best = estimate_modified(
        tlwe_block_params,
        ell,
        red_cost_model=red_cost_model,
        red_shape_model=red_shape_model,
        block_n=n,
        include_primal_block=include_primal_block,
        include_experimental_hybrid=include_experimental_hybrid,
    )

    if lwe_only:
        rlwe_original = math.inf
        rlwe_modified_best = cost_module.Cost(rop=oo, tag="skipped")
        rlwe_modified_costs = []
    else:
        rlwe_baseline_params = binary_params(
            N,
            profile.get("rlwe_q", profile["q"]),
            profile.get("rlwe_err_stddev", profile["err_stddev"]),
            f"binary-trlwe-baseline-ell{ell}-n{n}-N{N}",
            profile.get("baseline_secret", estimator.ND.Binary),
        )
        rlwe_mixed_params = binary_params(
            N,
            profile.get("rlwe_q", profile["q"]),
            profile.get("rlwe_err_stddev", profile["err_stddev"]),
            f"mixed-trlwe-ell{ell}-n{n}-N{N}",
            profile.get("baseline_secret", estimator.ND.Binary),
        )
        rlwe_original = estimate_original_bits(rlwe_baseline_params, red_cost_model, quiet=quiet)
        rlwe_modified_costs, rlwe_modified_best = estimate_modified(
            rlwe_mixed_params,
            ell,
            red_cost_model=red_cost_model,
            red_shape_model=red_shape_model,
            block_n=n,
            include_primal_block=False,
            include_experimental_hybrid=False,
        )

    original = min(tlwe_original, rlwe_original)

    if cost_bits(tlwe_modified_best) <= cost_bits(rlwe_modified_best):
        modified_best = tlwe_modified_best
        modified_source = "TLWE"
    else:
        modified_best = rlwe_modified_best
        modified_source = "TRLWE"

    mitm = mitm_bits(ell, n)
    modified_bits = cost_bits(modified_best)
    total = min(original, modified_bits, mitm)
    if total == original:
        limiting_tag = "Esti"
    elif total == modified_bits:
        limiting_tag = f"{modified_source}:{modified_best.get('tag', 'none')}"
    else:
        limiting_tag = "MitM"
    return {
        "ell": ell,
        "n": n,
        "original_bits": original,
        "tlwe_original_bits": tlwe_original,
        "rlwe_original_bits": rlwe_original,
        "modified_bits": modified_bits,
        "tlwe_modified_bits": cost_bits(tlwe_modified_best),
        "rlwe_modified_bits": cost_bits(rlwe_modified_best),
        "modified_tag": f"{modified_source}:{modified_best.get('tag', 'none')}",
        "limiting_tag": limiting_tag,
        "mitm_bits": mitm,
        "total_bits": total,
        "modified_costs": tlwe_modified_costs + rlwe_modified_costs,
    }


def print_row_table(rows, include_paper=False):
    if include_paper:
        print(
            "ell     n  Esti(calc)  Esti(paper)  Mod(calc)  Mod(paper)  "
            "MitM(calc)  MitM(paper)  limiting"
        )
        for row in rows:
            print(
                f"{row['ell']:>3} {row['n']:>5} "
                f"{row['original_bits']:>11.1f} {row['paper_original']:>12.1f} "
                f"{row['modified_bits']:>10.1f} {row['paper_modified']:>11.1f} "
                f"{row['mitm_bits']:>11.3f} {row['paper_mitm']:>12.3f} "
                f"{row.get('limiting_tag', row['modified_tag'])}"
            )
    else:
        print("ell     n  Esti(calc)  Mod(calc)  MitM(calc)  total  limiting")
        for row in rows:
            print(
                f"{row['ell']:>3} {row['n']:>5} "
                f"{row['original_bits']:>11.1f} "
                f"{row['modified_bits']:>10.1f} "
                f"{row['mitm_bits']:>11.3f} "
                f"{row['total_bits']:>6.1f} "
                f"{row.get('limiting_tag', row['modified_tag'])}"
            )


def run_paper_table(args):
    if args.quiet:
        estimator.Logging.set_level(estimator.Logging.QUIET)
    profile = paper_profile()
    print(f"profile: {profile['name']} ({profile['note']})")
    rows = []
    for ell, n, _N, paper_original, paper_modified, paper_mitm in PAPER_PARAMS:
        row = estimate_row(
            ell,
            n,
            profile,
            red_cost_model=estimator.RC.BDGL16,
            red_shape_model=conf.red_shape_model,
            quiet=args.quiet,
            include_primal_block=args.include_primal_block,
            include_experimental_hybrid=args.include_experimental_hybrid,
            lwe_only=False,
        )
        row.update(
            {
                "paper_original": paper_original,
                "paper_modified": paper_modified,
                "paper_mitm": paper_mitm,
            }
        )
        rows.append(row)
    print_row_table(rows, include_paper=True)


def run_single(args):
    if args.quiet:
        estimator.Logging.set_level(estimator.Logging.QUIET)
    profile = tfhe636_profile() if args.profile == "tfhe636" else paper_profile()
    row = estimate_row(
        args.ell,
        args.n,
        profile,
        red_cost_model=estimator.RC.BDGL16,
        red_shape_model=conf.red_shape_model,
        quiet=args.quiet,
        include_primal_block=args.include_primal_block,
        include_experimental_hybrid=args.include_experimental_hybrid,
        lwe_only=args.lwe_only,
    )
    print(f"profile: {profile['name']} ({profile['note']})")
    print_row_table([row])
    if args.verbose_costs:
        print()
        for cost in row["modified_costs"]:
            if is_finite_cost(cost):
                print(cost)
                print()


def run_search(args):
    estimator.Logging.set_level(estimator.Logging.QUIET)
    profile = tfhe636_profile() if args.profile == "tfhe636" else paper_profile()
    print(f"profile: {profile['name']} ({profile['note']})")
    rows = []
    for ell in args.ells:
        min_blocks_mitm = int(math.ceil(args.target / (0.28 * math.log2(ell + 1))))
        n = max(args.n_min, ell * min_blocks_mitm)
        if n % ell:
            n += ell - (n % ell)

        found = None
        while n <= args.n_max:
            row = estimate_row(
                ell,
                n,
                profile,
                red_cost_model=estimator.RC.BDGL16,
                red_shape_model=conf.red_shape_model,
                quiet=True,
                include_primal_block=args.include_primal_block,
                include_experimental_hybrid=args.include_experimental_hybrid,
                lwe_only=args.lwe_only,
            )
            if row["total_bits"] >= args.target:
                found = row
                break
            n += ell
        if found is None:
            found = {
                "ell": ell,
                "n": None,
                "original_bits": math.inf,
                "modified_bits": math.inf,
                "modified_tag": "not found",
                "mitm_bits": math.inf,
                "total_bits": math.inf,
            }
        rows.append(found)
    print_row_table(rows)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    paper = sub.add_parser("paper-table", help="reproduce the paper Table 1 parameter rows")
    paper.add_argument("--quiet", action="store_true", help="suppress lattice-estimator attack logs")
    paper.add_argument("--include-primal-block", action="store_true", help="also run the slow primal block-hybrid scan")
    paper.add_argument("--include-experimental-hybrid", action="store_true", help="also run an aggressive block-dual-hybrid model not used for paper reproduction")
    paper.set_defaults(func=run_paper_table)

    single = sub.add_parser("single", help="estimate one block-binary LWE parameter set")
    single.add_argument("--profile", choices=("paper", "tfhe636"), default="tfhe636")
    single.add_argument("--ell", type=int, required=True)
    single.add_argument("--n", type=int, required=True)
    single.add_argument("--quiet", action="store_true")
    single.add_argument("--verbose-costs", action="store_true")
    single.add_argument("--include-primal-block", action="store_true", help="also run the slow primal block-hybrid scan")
    single.add_argument("--include-experimental-hybrid", action="store_true", help="also run an aggressive block-dual-hybrid model not used for paper reproduction")
    single.add_argument("--lwe-only", action="store_true", help="skip the TRLWE-as-LWE estimate")
    single.set_defaults(func=run_single)

    search = sub.add_parser("search", help="find the smallest n for each ell at a target bit security")
    search.add_argument("--profile", choices=("paper", "tfhe636"), default="tfhe636")
    search.add_argument("--target", type=float, default=128.0)
    search.add_argument("--ells", type=int, nargs="+", default=[2, 3, 4, 5, 6])
    search.add_argument("--n-min", type=int, default=1)
    search.add_argument("--n-max", type=int, default=1400)
    search.add_argument("--include-primal-block", action="store_true", help="also run the slow primal block-hybrid scan")
    search.add_argument("--include-experimental-hybrid", action="store_true", help="also run an aggressive block-dual-hybrid model not used for paper reproduction")
    search.add_argument("--lwe-only", action="store_true", help="skip the TRLWE-as-LWE estimate")
    search.set_defaults(func=run_search)

    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
