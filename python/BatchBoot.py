#!/usr/bin/env python3
"""Conservative BatchBoot noise and LWE-security input generator.

This script follows the operation count in TFHEpp/include/tfhe/batchboot.hpp,
not merely the asymptotic count in the BatchBoot paper.  It is deliberately a
screening bound: it does not replace a distributional proof or measurements.

Run without Sage for a functional-BatchBoot noise screen:

  python3 BatchBoot.py --profile tfhepp-lvl1 --slots 1024 --hamming-weight 34

Run the matching LWE estimates with Sage/lattice-estimator:

  sage -python BatchBoot.py --profile tfhepp-lvl1 --security

The source and target instances are separate.  In particular, the BSK/EMP
selectors and automorphism keys are RLWE samples under the *target* secret;
the source sparse secret only protects the input ciphertext.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Instance:
    name: str
    n: int
    q_bits: int
    levels: int
    base_bits: int
    alpha_bits: int
    secret: str
    hamming_weight: int | None = None

    @property
    def alpha(self) -> float:
        return 2.0 ** -self.alpha_bits


@dataclass(frozen=True)
class Profile:
    name: str
    source: Instance
    target: Instance
    output_bits: int


# These are the native default structs selected by TFHEpp/include/params.hpp.
# `source` is made sparse binary because BatchBootKeyGen requires that shape;
# callers should override --hamming-weight with their public fixed weight.
PROFILES = {
    "tfhepp-lvl1": Profile(
        "tfhepp-lvl1",
        Instance("source", 1024, 32, 3, 6, 25, "sparse_binary", 34),
        Instance("target", 1024, 32, 3, 6, 25, "ternary"),
        3,
    ),
    "tfhepp-lvl2": Profile(
        "tfhepp-lvl2",
        Instance("source", 2048, 64, 4, 9, 51, "sparse_binary", 34),
        Instance("target", 2048, 64, 4, 9, 51, "ternary"),
        3,
    ),
    # This is the target used by TestD1EMP.  It is an EMP-only test profile;
    # it is not the default level-2 circuit-Bootstrap parameter set.
    "tfhepp-d1-emp": Profile(
        "tfhepp-d1-emp",
        Instance("source", 2048, 64, 1, 23, 53, "sparse_binary", 34),
        Instance("target", 2048, 64, 1, 23, 53, "ternary"),
        3,
    ),
}


def log2_or_neg_inf(value: float) -> float:
    return math.log2(value) if value else float("-inf")


def secret_second_moment(instance: Instance) -> float:
    if instance.secret == "ternary":
        # Uniform {-1,0,1}, as used by TFHEpp's default lvl1/lvl2 params.
        return 2.0 / 3.0
    if instance.secret == "sparse_binary":
        assert instance.hamming_weight is not None
        return instance.hamming_weight / instance.n
    raise ValueError(f"unknown secret distribution: {instance.secret}")


def decomposition_round_variance(instance: Instance) -> float:
    """Normalized variance from the unrepresented low gadget bits.

    This is the usual uniform-rounding screen.  It is zero when the primary
    decomposition covers the torus word; it intentionally does not model
    floating-point FFT error.
    """
    missing = instance.q_bits - instance.levels * instance.base_bits
    if missing <= 0:
        return 0.0
    return 2.0 ** (-2 * missing) / 12.0


def external_product_increment(instance: Instance, message_bound_sq: float) -> float:
    """A normalized, conservative external-product variance increment.

    The fresh-key term has one nonce and one body decomposition for k=1.  The
    rounding term bounds both decompositions and weights nonce rounding by the
    target secret's second moment.  This matches the standard-decomposition
    path used by BatchEMP/BatchExternalProductAddFD.
    """
    b_sq_plus_two = 2.0 ** (2 * instance.base_bits) + 2.0
    fresh = (
        2.0
        * instance.levels
        * instance.n
        * b_sq_plus_two
        / 12.0
        * instance.alpha**2
    )
    rounding = decomposition_round_variance(instance) * (
        instance.n * secret_second_moment(instance) + 1.0
    )
    return fresh + message_bound_sq * rounding


def automorphism_increment(instance: Instance, message_bound_sq: float) -> float:
    # EvalAuto uses a HalfTRGSW key (one body decomposition).  Counting it as
    # a full EP is deliberately conservative and avoids claiming a tighter
    # bound for the current coefficient-domain implementation.
    return external_product_increment(instance, message_bound_sq)


def normal_tail_log2(z: float) -> float:
    """log2(Pr[|X| >= z sigma]) without underflow for large z."""
    if z <= 0:
        return 0.0
    p = math.erfc(z / math.sqrt(2.0))
    if p:
        return math.log2(p)
    # Mills-ratio upper approximation, sufficient for a reportable screen.
    return (-z * z / 2.0 - math.log(z) - 0.5 * math.log(2.0 * math.pi)) / math.log(2.0)


def noise_screen(profile: Profile, slots: int, components: int) -> dict[str, float | int]:
    target = profile.target
    source = profile.source
    if slots < 2 or slots & (slots - 1) or slots > target.n:
        raise ValueError("slots must be a power of two in [2, target.n]")
    if components < 1:
        raise ValueError("components must be positive")

    # BatchEMPKeyGen stores ceil(log_4(slots)) digits and BatchEMP evaluates
    # four selector external products per digit.  This is exact for the C++
    # code, including odd log2(slots), rather than the paper's 2*ell shorthand.
    emp_digits = (slots.bit_length() - 1 + 1) // 2
    eps_per_emp = 4 * emp_digits
    h = source.hamming_weight or 0
    emp_calls = h + components  # one final_positive key per non-empty component
    emp_increment = external_product_increment(target, 1.0)
    auto_increment = automorphism_increment(target, 1.0)

    # BatchRepack uses log2(slots) automorphisms and BatchHomomorphicTrace
    # supplies the remaining log2(N/slots), for exactly log2(N) in total.
    automorphisms = target.n.bit_length() - 1
    output_variance = emp_calls * eps_per_emp * emp_increment + automorphisms * auto_increment
    output_std = math.sqrt(output_variance)
    decoding_margin = 2.0 ** (-(profile.output_bits + 1))
    failure_log2 = normal_tail_log2(decoding_margin / output_std)

    # Normalized source error after ModSwitch(q_source -> 2N_target), measured
    # in exponent bins.  Values much smaller than 1/2 leave a rounding margin;
    # this is a diagnostic, not a complete input-correctness proof.
    source_modswitch_sigma_bins = 2.0 * target.n * source.alpha

    return {
        "emp_digits": emp_digits,
        "external_products_per_emp": eps_per_emp,
        "emp_calls": emp_calls,
        "automorphisms": automorphisms,
        "ep_increment": emp_increment,
        "auto_increment": auto_increment,
        "output_variance": output_variance,
        "output_std": output_std,
        "failure_log2": failure_log2,
        "source_modswitch_sigma_bins": source_modswitch_sigma_bins,
    }


def print_screen(profile: Profile, slots: int, components: int) -> None:
    r = noise_screen(profile, slots, components)
    source, target = profile.source, profile.target
    print(f"profile: {profile.name}")
    print(
        f"source: n={source.n}, q=2^{source.q_bits}, {source.secret}, "
        f"h={source.hamming_weight}, alpha=2^-{source.alpha_bits}"
    )
    print(
        f"target: n={target.n}, q=2^{target.q_bits}, {target.secret}, "
        f"l={target.levels}, Bg=2^{target.base_bits}, alpha=2^-{target.alpha_bits}"
    )
    print(
        "BatchEMP: "
        f"{r['emp_digits']} radix-4 digits x 4 EP = {r['external_products_per_emp']} EP/call; "
        f"{r['emp_calls']} calls; repack/trace={r['automorphisms']} automorphisms"
    )
    print(f"EP variance increment (normalized): 2^{log2_or_neg_inf(r['ep_increment']):.2f}")
    print(f"output variance bound (normalized): 2^{log2_or_neg_inf(r['output_variance']):.2f}")
    print(f"output standard deviation bound: 2^{log2_or_neg_inf(r['output_std']):.2f}")
    print(f"per-coefficient Gaussian failure screen: <= 2^{r['failure_log2']:.1f}")
    print(
        "source ModSwitch sigma: "
        f"{r['source_modswitch_sigma_bins']:.3g} exponent bins (rounding threshold is 0.5)"
    )


def run_security(profile: Profile, role: str) -> None:
    """Run Sage's lattice-estimator for both security-relevant secrets."""
    try:
        import importlib

        estimator = importlib.import_module(".estimator", "lattice-estimator")
    except Exception as error:
        raise SystemExit(
            "Sage/lattice-estimator is required for --security. Run with "
            "`sage -python BatchBoot.py ... --security`.\n"
            f"Import error: {error}"
        )

    nd = estimator.nd
    LWEParameters = estimator.lwe_parameters.LWEParameters
    instances = (profile.source, profile.target)
    if role != "both":
        instances = (profile.source,) if role == "source" else (profile.target,)
    for instance in instances:
        if instance.secret == "ternary":
            xs = nd.UniformMod(3)
        else:
            h = instance.hamming_weight
            assert h is not None
            # SparseBinary has exactly h ones and n-h zeroes.
            xs = nd.SparseTernary(h, 0, instance.n)
        param = LWEParameters(
            n=instance.n,
            q=2**instance.q_bits,
            Xs=xs,
            Xe=nd.DiscreteGaussian(stddev=2 ** (instance.q_bits - instance.alpha_bits)),
            m=instance.n,
            tag=f"BatchBoot-{profile.name}-{instance.name}",
        )
        print("=" * 72)
        print(param)
        results = estimator.LWE.estimate(param, red_cost_model=estimator.RC.BDGL16)
        for attack, result in results.items():
            rop = result.get("rop") if hasattr(result, "get") else None
            if rop is not None:
                # `rop` is an absolute operation count, often much larger
                # than a Python float.  Take its logarithm before coercion.
                print(f"{attack}: log2(rop)={float(math.log(rop, 2)):.2f}")
            else:
                print(f"{attack}: {result}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=PROFILES, default="tfhepp-lvl1")
    parser.add_argument("--slots", type=int, help="Active BatchBoot slots (default: target N)")
    parser.add_argument("--hamming-weight", type=int, help="Public source sparse-binary weight")
    parser.add_argument("--source-alpha-bits", type=int,
                        help="Use source fresh noise alpha=2^-bits")
    parser.add_argument("--target-alpha-bits", type=int,
                        help="Use target/evaluation-key fresh noise alpha=2^-bits")
    parser.add_argument("--components", type=int, default=1, help="Non-empty source MLWE components")
    parser.add_argument("--output-bits", type=int, help="Decoded output width")
    parser.add_argument("--security", action="store_true", help="Run lattice-estimator (requires Sage)")
    parser.add_argument("--security-role", choices=("both", "source", "target"),
                        default="both", help="LWE instance(s) to estimate")
    args = parser.parse_args()

    profile = PROFILES[args.profile]
    source = profile.source
    if args.hamming_weight is not None:
        if not 1 <= args.hamming_weight <= source.n:
            parser.error("--hamming-weight must be in [1, source.n]")
        source = replace(source, hamming_weight=args.hamming_weight)
    if args.source_alpha_bits is not None:
        if args.source_alpha_bits < 1:
            parser.error("--source-alpha-bits must be positive")
        source = replace(source, alpha_bits=args.source_alpha_bits)
    target = profile.target
    if args.target_alpha_bits is not None:
        if args.target_alpha_bits < 1:
            parser.error("--target-alpha-bits must be positive")
        target = replace(target, alpha_bits=args.target_alpha_bits)
    if args.output_bits is not None:
        if args.output_bits < 1:
            parser.error("--output-bits must be positive")
        profile = replace(profile, output_bits=args.output_bits)
    profile = replace(profile, source=source, target=target)
    slots = args.slots or profile.target.n

    print_screen(profile, slots, args.components)
    if args.security:
        run_security(profile, args.security_role)


if __name__ == "__main__":
    main()
