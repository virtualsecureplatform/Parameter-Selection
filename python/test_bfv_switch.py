import math
import random
import unittest
from dataclasses import replace

from noiseestimation.bfv_switch import (
    Parameters, analyse, envelope, product_envelope, reverse_reference, search,
)


class BFVSwitchTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
