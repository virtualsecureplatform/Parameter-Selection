from __future__ import annotations

from pathlib import Path
import re
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compact_cover_bgv_certificate import (  # noqa: E402
    PLAINTEXT_SQUARE,
    RNS_PRIMES_AND_ROOTS,
    build_certificate,
)
from noiseestimation.bfv_polynomial import (  # noqa: E402
    lowest_digit_removal_polynomial_over_range,
)


class CompactCoverBGVCertificateTest(unittest.TestCase):
    def test_certificate(self) -> None:
        certificate = build_certificate()
        self.assertTrue(certificate.correctness_certified)
        self.assertTrue(certificate.estimated_security_meets_target)
        self.assertEqual(certificate.modulus_bits, 1402)
        self.assertEqual(certificate.digit_error_bound, 23)
        self.assertEqual(certificate.digit_polynomial_degree, 93)
        self.assertEqual(
            certificate.digit_polynomial_coefficients,
            tuple(lowest_digit_removal_polynomial_over_range(65_537, 23)),
        )
        self.assertEqual(
            certificate.trace_exponents,
            (5, 25, 625, 128_481, 28_609, 61_313, 7_937, 81_409,
             31_745, 63_489, 126_977, 122_881, 114_689, 98_305,
             65_537, 131_071),
        )
        self.assertEqual(certificate.output_limbs, 13)
        self.assertEqual(certificate.accepted_input_error_bound, 3_215_720_448)
        self.assertEqual(certificate.projected_error_bound, 70_304_015_244_276)
        self.assertEqual(certificate.addition_output_error_bound, 36)
        self.assertEqual(certificate.multiplication_output_error_bound, 18)
        self.assertLessEqual(certificate.addition_output_error_bound,
                             certificate.accepted_input_error_bound)
        self.assertLessEqual(certificate.multiplication_output_error_bound,
                             certificate.accepted_input_error_bound)
        self.assertLess(certificate.output_error_log2_bound,
                        certificate.output_capacity_log2)
        self.assertLess(certificate.multiplication_output_error_log2_bound,
                        certificate.multiplication_capacity_log2)
        self.assertGreater(certificate.cycle_contraction_bits, 0)
        self.assertGreaterEqual(certificate.retained_security_proxy_bits, 128)

    def test_rns_congruences(self) -> None:
        modulus = 2 * 65_536 * PLAINTEXT_SQUARE
        for prime, _ in RNS_PRIMES_AND_ROOTS:
            self.assertEqual((prime - 1) % modulus, 0)

    def test_certificate_hash_is_stable(self) -> None:
        self.assertEqual(
            build_certificate().sha256(),
            "f3b1e9f169d152bdab7d17305d54e881d20fd022aed273cbc0640afe946a4e73",
        )

    def test_cross_repository_parameter_alignment(self) -> None:
        workspace = Path(__file__).resolve().parents[3]
        tfhepp = workspace / "TFHEpp"
        formal = workspace / "FormalProof4FHE"
        if not tfhepp.is_dir() or not formal.is_dir():
            self.skipTest("sibling TFHEpp/FormalProof4FHE repositories unavailable")

        ntt_source = (tfhepp / "include/modular_ntt.hpp").read_text()
        ntt_block = re.search(
            r"degree65536_primes\{\{(.*?)\}\};", ntt_source, re.DOTALL
        )
        self.assertIsNotNone(ntt_block)
        cpp_pairs = tuple(
            (int(value), int(root))
            for value, root in re.findall(
                r"\{(\d+)ULL,\s*(\d+)\}", ntt_block.group(1)
            )
        )
        self.assertEqual(cpp_pairs, RNS_PRIMES_AND_ROOTS)
        cpp_primes = tuple(value for value, _ in cpp_pairs)

        lean_source = (formal / "FormalProof4FHE/RLWE/CompactCoverBGVExactNoise.lean").read_text()
        lean_block = re.search(
            r"def rnsPrimes : List ℕ :=\s*\[(.*?)\]",
            lean_source,
            re.DOTALL,
        )
        self.assertIsNotNone(lean_block)
        lean_primes = tuple(int(value) for value in re.findall(r"\d+", lean_block.group(1)))
        self.assertEqual(lean_primes, cpp_primes)
        lean_prime_data = tuple(
            (int(value), int(root))
            for value, root in re.findall(
                r"\{ value := (\d+), generator := (\d+),", lean_source
            )
        )
        self.assertEqual(lean_prime_data, cpp_pairs)

        instantiation_source = (
            formal / "FormalProof4FHE/RLWE/CompactCoverBGV65536Instantiation.lean"
        ).read_text()
        trace_block = re.search(
            r"def concreteTraceExponents : List ℕ :=\s*\[(.*?)\]",
            instantiation_source,
            re.DOTALL,
        )
        self.assertIsNotNone(trace_block)
        lean_trace = tuple(
            int(value) for value in re.findall(r"\d+", trace_block.group(1))
        )
        self.assertEqual(lean_trace, build_certificate().trace_exponents)

        coefficient_block = re.search(
            r"def digitRemovalCoefficients : List ℕ :=\s*\[(.*?)\]\s*\n\n/-!",
            lean_source,
            re.DOTALL,
        )
        self.assertIsNotNone(coefficient_block)
        lean_coefficients = tuple(
            int(value) for value in re.findall(r"\d+", coefficient_block.group(1))
        )
        self.assertEqual(
            lean_coefficients,
            tuple(lowest_digit_removal_polynomial_over_range(65_537, 23)),
        )

        bootstrap_source = (
            tfhepp / "include/bfv/compact-cover-bgv-bootstrap.hpp"
        ).read_text()
        self.assertIn('certificate_version = 5', bootstrap_source)
        self.assertIn(build_certificate().sha256(), bootstrap_source)
        polynomial_hash = 1_469_598_103_934_665_603
        for coefficient in lean_coefficients:
            for byte in coefficient.to_bytes(8, "little"):
                polynomial_hash ^= byte
                polynomial_hash = (
                    polynomial_hash * 1_099_511_628_211
                ) & ((1 << 64) - 1)
        self.assertIn(f"0x{polynomial_hash:016x}", bootstrap_source)


if __name__ == "__main__":
    unittest.main()
