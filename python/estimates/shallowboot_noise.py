#!/usr/bin/env python3
"""Average-case noise estimator/search for ePrint 2026/1730 Algorithm 3.

It follows PBC aggregation, every balanced QH-RLWE product layer, and the
complete modulus chain. For monomial messages and independent coefficient
errors, Vproduct = Vleft + Vright + N*Vleft*Vright.
"""
from __future__ import annotations
import argparse
import json
import math
from dataclasses import asdict, dataclass

@dataclass
class Layer:
    layer: int
    outputs: int
    modulus_bits: int
    log2_variance_before_switch: float
    switched_to_bits: int | None
    log2_variance_after_switch: float
    wrap_headroom_bits: float

def multiply_variance(left: float, right: float, degree: int) -> float:
    return left + right + degree * left * right

def balanced_layer(values: list[float], degree: int) -> list[float]:
    output = [multiply_variance(values[i], values[i + 1], degree)
              for i in range(0, len(values) - 1, 2)]
    if len(values) & 1:
        output.append(values[-1])
    return output

def bucket_entry_counts(n: int, copies: int, buckets: int) -> list[int]:
    quotient, remainder = divmod(n * copies, buckets)
    return [quotient + (i < remainder) + 1 for i in range(buckets)]

def tail_log2(margin: float, variance: float) -> float:
    probability = math.erfc(margin / math.sqrt(2 * variance))
    return math.log2(probability) if probability else -math.inf

def estimate(args: argparse.Namespace, middle_bits: int | None = None) -> dict:
    buckets = args.lwe_h + args.pbc_slack
    counts = bucket_entry_counts(args.lwe_n, args.pbc_copies, buckets)
    # RLWE_lsb phase is m+t*e, with t the plaintext modulus.
    fresh_variance = (args.plaintext_modulus * args.rlwe_sigma) ** 2
    values = [count * fresh_variance for count in counts]
    chain = list(args.modulus_chain)
    if middle_bits is not None:
        chain[-2] = middle_bits
    switches = dict(zip(args.switch_after, chain[1:]))
    current_bits = chain[0]
    layers: list[Layer] = []
    depth = 0
    while len(values) > 1:
        depth += 1
        input_bits = current_bits
        values = balanced_layer(values, args.ring_n)
        before = max(values)
        headroom = current_bits - 1 - math.log2(before) / 2
        target = switches.get(depth)
        if target is not None:
            ratio = 2.0 ** (target - current_bits)
            # Phase rounding is r_b-r_a*s for a coefficient-small secret.
            rounding = (1 + args.secret_second_moment * args.ring_n) / 12
            values = [value * ratio * ratio + rounding for value in values]
            current_bits = target
        layers.append(Layer(depth, len(values), input_bits, math.log2(before),
                            target, math.log2(max(values)), headroom))

    final_ratio = 2.0 ** (args.lwe_qbits - current_bits)
    selector_variance = values[0] * final_ratio * final_ratio
    # Multiplication by t^{-1} converts LSB to MSB and divides the error by t.
    # A dense binary sign LUT then sums approximately density*N independently
    # shifted error coefficients into each accumulator coefficient.
    final_variance = (selector_variance / args.plaintext_modulus**2 *
                      args.lut_density * args.ring_n)
    margin = 2 ** (args.lwe_qbits - 3)
    coefficient_failure = tail_log2(margin, final_variance)
    union_failure = min(0.0, coefficient_failure + math.log2(args.ring_n))
    return {
        "lwe_h": args.lwe_h, "pbc_buckets": buckets,
        "entries_per_bucket_min": min(counts),
        "entries_per_bucket_max": max(counts), "ring_n": args.ring_n,
        "rlwe_sigma": args.rlwe_sigma,
        "plaintext_modulus": args.plaintext_modulus,
        "modulus_chain": chain, "switch_after": args.switch_after,
        "selector_log2_variance": math.log2(selector_variance),
        "selector_log2_stddev": math.log2(selector_variance) / 2,
        "final_log2_variance": math.log2(final_variance),
        "final_log2_stddev": math.log2(final_variance) / 2,
        "coefficient_log2_failure": coefficient_failure,
        "ring_union_log2_failure": union_failure,
        "minimum_wrap_headroom_bits": min(x.wrap_headroom_bits for x in layers),
        "passes_128bit_noise_screen": union_failure <= -128,
        "layers": [asdict(layer) for layer in layers],
    }

def parse_range(text: str) -> range:
    fields = [int(value) for value in text.split(":")]
    if len(fields) != 2 or fields[0] >= fields[1]:
        raise argparse.ArgumentTypeError("range must be START:STOP")
    return range(fields[0], fields[1] + 1)

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lwe-n", type=int, default=1450)
    parser.add_argument("--lwe-h", type=int, default=37)
    parser.add_argument("--lwe-qbits", type=int, default=9)
    parser.add_argument("--pbc-copies", type=int, default=3)
    parser.add_argument("--pbc-slack", type=int, default=3)
    parser.add_argument("--ring-n", type=int, default=8192)
    parser.add_argument("--rlwe-sigma", type=float, default=0.75)
    parser.add_argument("--plaintext-modulus", type=int, default=4)
    parser.add_argument("--lut-density", type=float, default=0.5)
    parser.add_argument("--secret-second-moment", type=float, default=2 / 3)
    parser.add_argument("--modulus-chain", type=int, nargs="+",
                        default=[151, 69, 52, 36])
    parser.add_argument("--switch-after", type=int, nargs="+",
                        default=[3, 4, 5])
    parser.add_argument("--search-middle", type=parse_range)
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()

def main() -> None:
    args = arguments()
    if len(args.modulus_chain) != len(args.switch_after) + 1:
        raise SystemExit("modulus-chain needs one more entry than switch-after")
    if any(a <= b for a, b in zip(args.modulus_chain, args.modulus_chain[1:])):
        raise SystemExit("modulus-chain must be strictly decreasing")
    if args.search_middle:
        candidates = [estimate(args, bits) for bits in args.search_middle
                      if args.modulus_chain[-1] < bits < args.modulus_chain[-3]]
        candidates.sort(key=lambda item: item["ring_union_log2_failure"])
        if args.json:
            print(json.dumps(candidates, indent=2)); return
        print("middle  log2(final sigma)  log2(ring failure)  wrap headroom  pass")
        for item in candidates:
            print(f"{item['modulus_chain'][-2]:6d}  "
                  f"{item['final_log2_stddev']:17.2f}  "
                  f"{item['ring_union_log2_failure']:18.1f}  "
                  f"{item['minimum_wrap_headroom_bits']:13.1f}  "
                  f"{'yes' if item['passes_128bit_noise_screen'] else 'no'}")
        return
    result = estimate(args)
    if args.json:
        print(json.dumps(result, indent=2)); return
    print("Algorithm 3 PBC/QH end-to-end accumulator noise estimate")
    print(f"  h={args.lwe_h}, buckets={result['pbc_buckets']}, entries/factor="
          f"{result['entries_per_bucket_min']}..{result['entries_per_bucket_max']}")
    print(f"  modulus chain={result['modulus_chain']}, switches after layers={args.switch_after}")
    for layer in result["layers"]:
        switch = (f" -> 2^{layer['switched_to_bits']}"
                  if layer["switched_to_bits"] else "")
        print(f"  layer {layer['layer']}: outputs={layer['outputs']}, "
              f"log2(V)={layer['log2_variance_before_switch']:.2f}{switch}, "
              f"post-switch={layer['log2_variance_after_switch']:.2f}, "
              f"wrap headroom={layer['wrap_headroom_bits']:.1f} bits")
    print(f"  selector log2(sigma) at q=2^{args.lwe_qbits}: "
          f"{result['selector_log2_stddev']:.2f}")
    print(f"  post-conversion/LUT log2(sigma): {result['final_log2_stddev']:.2f}")
    print(f"  ring-union failure estimate: 2^{result['ring_union_log2_failure']:.1f}")
    print("  end-to-end accumulator 128-bit noise screen: "
          f"{'PASS' if result['passes_128bit_noise_screen'] else 'FAIL'}")

if __name__ == "__main__":
    main()
