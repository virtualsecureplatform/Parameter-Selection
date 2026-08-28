from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compact_cover_bgv_certificate import (  # noqa: E402
    ACCEPTED_INPUT_LOG2_VARIANCE,
    FAILURE_TARGET_LOG2,
    PLAINTEXT_SQUARE,
    RNS_PRIMES_AND_ROOTS,
    build_certificate,
)


class CompactCoverBGVCertificateTest(unittest.TestCase):
    def test_certificate(self) -> None:
        certificate = build_certificate()
        self.assertTrue(certificate.certified)
        self.assertEqual(certificate.modulus_bits, 915)
        self.assertEqual(certificate.digit_error_bound, 43)
        self.assertEqual(certificate.digit_polynomial_degree, 173)
        self.assertLessEqual(
            certificate.rounding_failure_log2_bound, FAILURE_TARGET_LOG2
        )
        self.assertLess(
            certificate.fresh_output_log2_variance,
            ACCEPTED_INPUT_LOG2_VARIANCE,
        )
        self.assertGreaterEqual(certificate.certified_security_bits, 128)
        self.assertLess(
            certificate.deterministic_output_error_log2_bound,
            certificate.modulus_log2 - 1,
        )
        self.assertEqual(certificate.bootstrap_key_bytes, 75 * (1 << 20) + 48)

    def test_rns_congruences(self) -> None:
        modulus = 2 * 65_536 * PLAINTEXT_SQUARE
        for prime, _ in RNS_PRIMES_AND_ROOTS:
            self.assertEqual((prime - 1) % modulus, 0)

    def test_certificate_hash_is_stable(self) -> None:
        self.assertEqual(
            build_certificate().sha256(),
            "143be2a21f31cd5bb2bb3d6359b81f1103a05329639fd42152c3503669dd39b1",
        )


if __name__ == "__main__":
    unittest.main()
