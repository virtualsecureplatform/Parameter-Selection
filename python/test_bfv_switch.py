import math
import random
import unittest
from dataclasses import replace

from noiseestimation.bfv_switch import (
    Parameters, analyse, envelope, product_envelope, reverse_reference, search,
    parameters_for_width, tune_16bit, key_bytes,
)


class BFVSwitchTests(unittest.TestCase):
    def test_paired_forward(self):
        p = parameters_for_width(16)
        report = analyse(p)
        old = analyse(replace(p, forward_digit_bits=1))
        self.assertEqual(report["events"]["pbs3"], 16 * p.products)
        self.assertEqual(report["events"]["pbs1"], 80 * p.products)
        self.assertEqual(report["events"]["ks1h"], 81 * p.products)
        self.assertLess(report["envelope"]["forward_error"], old["envelope"]["forward_error"])
        # Integer negacyclic LUT reference, including both ends of every
        # strict address interval, for all eight weighted radix-4 digits.
        n, q = p.ring_n, 1 << p.qbits
        for pair in range(8):
            half_weight = (q >> 33) << (2 * pair)
            for digit in range(4):
                center = (2 * digit - 3) * n // 8
                for displacement in (-n // 8 + 1, 0, n // 8 - 1):
                    index = (center + displacement) % (2 * n)
                    coefficient = index % n
                    value = half_weight * (1 if coefficient < n // 4 or coefficient >= 3 * n // 4 else 3)
                    if index >= n:
                        value = -value
                    self.assertEqual(value + 3 * half_weight, digit * 2 * half_weight)
        with self.assertRaises(ValueError):
            analyse(replace(p, forward_digit_bits=3))

    def test_all_uint16_reverse_values(self):
        for message in range(65536):
            self.assertEqual(reverse_reference(message), message)

    def test_signed_product_error_and_carry_boundaries(self):
        params = Parameters()
        report = analyse(params)
        radius = int(report["envelope"]["product_error"] * 2**128)
        values = {0, 1, 255, 65025, 65535}
        for bit in range(1, 16):
            values.update({2**bit - 1, 2**bit, 2**bit + 1})
        for value in sorted(values):
            for error in (-radius, 0, radius):
                self.assertEqual(reverse_reference(value, error), value)

    def test_primitive_error_injections(self):
        p = Parameters()
        errors = analyse(p)["envelope"]["errors"]
        randomizer = random.Random(20260919)
        mapping = {"ks31": "ks31", "ks1h": "ks1h", "extract_ks": "ks1h",
                   "br1": "br1", "low_br": "br1", "high_br": "br1", "pbs1": "pbs1"}
        for _ in range(2048):
            value = randomizer.randrange(65536)
            injections = {(stage, digit): randomizer.choice((-1, 1)) * int(errors[name] * 2**128)
                          for stage, name in mapping.items() for digit in range(9)}
            self.assertEqual(reverse_reference(value, injections=injections), value)

    def test_target_and_union_bound(self):
        report = analyse()
        self.assertTrue(report["meets_target_conditionally"])
        expected = math.log2(2 * report["total_events"]) - report["cutoff_with_slack"]**2 / (2 * math.log(2))
        self.assertAlmostEqual(report["conditional_log2_failure"], expected)
        doubled = analyse(Parameters(products=512))
        self.assertAlmostEqual(doubled["conditional_log2_failure"], expected + 1)
        self.assertFalse(envelope(Parameters(), report["maximum_uniform_cutoff"] * 1.01)["passes"])

    def test_fail_closed(self):
        self.assertFalse(analyse(Parameters(boolean_fft_error=0.25))["meets_target_conditionally"])
        self.assertFalse(analyse(Parameters(to_io_levels=1, to_io_basebit=1))["meets_target_conditionally"])
        with self.assertRaises(ValueError):
            analyse(Parameters(products=0))
        with self.assertRaises(ValueError):
            analyse(Parameters(to_io_levels=10, to_io_basebit=8))
        self.assertTrue(math.isinf(product_envelope(1, 0, Parameters())))

    def test_search_and_multiplication_lift_terms(self):
        self.assertEqual(search()["parameters"], analyse()["parameters"])
        p = Parameters()
        e = 2.0**-70
        bound = product_envelope(e, 0, p)
        self.assertGreater(bound, 2 * 255 * e + p.ring_n * 65536 * e * e)

    def test_wider_profiles_and_reference(self):
        rng = random.Random(20260920)
        for bits in (9, 10, 12, 16, 24, 32):
            p = Parameters(plaintext_bits=2 * bits)
            report = analyse(p)
            self.assertEqual(report["events"]["input"], 2 * bits * p.products)
            self.assertEqual(report["events"]["ks31"], (bits + 1) * p.products)
            self.assertEqual(report["evaluation_key_bytes"], analyse()["evaluation_key_bytes"])
            values = {0, 1, (1 << (2 * bits)) - 1}
            for bit in range(1, 2 * bits):
                values.update({(1 << bit) - 1, 1 << bit, (1 << bit) + 1})
            values.update(rng.randrange(1 << (2 * bits)) for _ in range(256))
            # A small signed error tests semantic guard/carry handling even
            # when the cryptographic envelope cannot certify this width.
            radius = 1 << (128 - 2 * bits - 9)
            for value in values:
                for error in (-radius, 0, radius):
                    self.assertEqual(reverse_reference(value, error, plaintext_bits=2 * bits), value)
        self.assertTrue(analyse(Parameters(plaintext_bits=18))["meets_target_conditionally"])
        self.assertFalse(analyse(Parameters(plaintext_bits=20))["meets_target_conditionally"])
        for invalid in (0, 3, 66):
            with self.assertRaises(ValueError):
                analyse(Parameters(plaintext_bits=invalid))

    def test_tuned_16bit_profile(self):
        self.assertEqual(parameters_for_width(8), Parameters())
        self.assertFalse(analyse(parameters_for_width(16, legacy=True))["meets_target_conditionally"])
        p = parameters_for_width(16)
        report = analyse(p)
        self.assertTrue(report["meets_target_conditionally"])
        self.assertLess(report["conditional_log2_failure"], -48)
        self.assertLess(key_bytes(p) + 4 * 1024**3, 24 * 1024**3)
        self.assertEqual(tune_16bit()["parameters"], analyse(replace(p, forward_digit_bits=1))["parameters"])
        radius = int(report["envelope"]["product_error"] * 2**128)
        rng = random.Random(20260921)
        values = {0, 1, 65535**2, 2**32 - 1}
        for bit in range(1, 32):
            values.update({2**bit - 1, 2**bit, 2**bit + 1})
        values.update(rng.randrange(2**32) for _ in range(2048))
        for value in values:
            for error in (-radius, 0, radius):
                self.assertEqual(reverse_reference(value, error, plaintext_bits=32), value)


if __name__ == "__main__":
    unittest.main()
