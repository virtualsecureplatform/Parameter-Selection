from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compact_cover_bgv_certificate import (  # noqa: E402
    PLAINTEXT_SQUARE,
    RNS_PRIMES_AND_ROOTS,
    build_certificate,
)


class CompactCoverBGVCertificateTest(unittest.TestCase):
    def test_certificate(self) -> None:
        certificate = build_certificate()
        self.assertTrue(certificate.certified)
        self.assertEqual(certificate.modulus_bits, 1402)
        self.assertEqual(certificate.digit_error_bound, 23)
        self.assertEqual(certificate.digit_polynomial_degree, 93)
        self.assertEqual(certificate.output_limbs, 13)
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
            "09c36cf48f9ac9b51cbd0a454b109e563663d3dc48a3f141fac7637cff13f16c",
        )


if __name__ == "__main__":
    unittest.main()
