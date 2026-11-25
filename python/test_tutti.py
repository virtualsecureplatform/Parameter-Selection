import math
import unittest
from dataclasses import replace

from noiseestimation.block_binary_meet_lwe import (
    AnchoredBlockBinary, meet_lwe_guessing, meet_lwe_guessing_free, odlyzko_bits)
from noiseestimation.tutti import (
    PROFILES, TuttiParams, budget, estimate, minimax_sine_scale, published_caps, search,
    truncated_normal_rms)


class BlockBinaryCombinatorialTests(unittest.TestCase):
    def test_reproduces_published_costs(self):
        # ePrint 2026/2848, Tables 3 and 10: (N/R, h, log q0) ->
        # Odlyzko lists, guessing-free lists, guessing-free operations, Meet-LWE.
        rows = {(2048, 64, 27): (160, 157.96, 165.26, 167.12),
                (8192, 32, 27): (128, 125.12, 132.09, 132.83)}
        for (degree, blocks, log_q), expected in rows.items():
            instance = AnchoredBlockBinary(degree, blocks)
            free = meet_lwe_guessing_free(instance, log_q)
            measured = (odlyzko_bits(instance), free["bits"], free["ops_bits"],
                        meet_lwe_guessing(instance, log_q)["bits"])
            for got, want in zip(measured, expected):
                self.assertAlmostEqual(got, want, delta=0.01)

    def test_anchor_removes_one_block(self):
        instance = AnchoredBlockBinary(2048, 64)
        self.assertEqual((instance.dimension, instance.weight), (2016, 63))


class TuttiModelTests(unittest.TestCase):
    def test_truncated_sampler_rms(self):
        self.assertAlmostEqual(truncated_normal_rms(3.2), 2.8283017, places=6)
        self.assertGreaterEqual(truncated_normal_rms(3.6), 3.2)

    def test_minimax_sine_equioscillates(self):
        gap, bound = 16.0, 1.0
        scale, error = minimax_sine_scale(gap, bound)
        g = gap / (2 * math.pi)
        worst = max(abs(scale * g * math.sin(v / g) - v)
                    for v in (bound * i / 4096 for i in range(4097)))
        self.assertAlmostEqual(worst, error, delta=1e-6)
        plain = g * (bound / g - math.sin(bound / g))
        self.assertLess(error, plain / 3)

    def test_ll13_budget_matches_tfhepp_schedule(self):
        p = PROFILES["ll13"]
        exposure = budget(p)
        self.assertEqual((p.candidates, p.indicators, p.depth), (32, 2016, 6))
        self.assertEqual(exposure["input_bits"], 38)
        self.assertEqual(exposure["factor_bits"], 192)
        self.assertEqual(exposure["selector_bits"], 212)
        self.assertEqual(exposure["maximum_dense_key_bits"], 212)
        self.assertEqual(exposure["dense_key_exposures"]["output_conjugation"], 78)
        self.assertEqual(exposure["sparse_key_exposure_bits"], 24)
        self.assertTrue(exposure["one_multiplication_returns_to_input"])
        caps = published_caps(p)
        self.assertTrue(caps["dense_within_cap"] and caps["sparse_within_cap"])
        self.assertTrue(caps["integer_noise_rms_matches_model"])

    def test_boolean_refresh_budget(self):
        p = PROFILES["ll13_boolean"]
        exposure = budget(p)
        self.assertEqual(exposure["selector_bits"], 200)
        self.assertEqual(exposure["maximum_dense_key_bits"], 200)
        self.assertEqual(exposure["dense_key_exposures"]["s2c_rotations"], 66)
        self.assertEqual(exposure["dense_key_exposures"]["tlwe_packing"], 48)
        self.assertEqual(exposure["sparse_key_exposure_bits"], 24)
        self.assertEqual(exposure["multiplications_returning_to_input"], 0)
        self.assertTrue(published_caps(p)["dense_within_cap"])
        with self.assertRaises(ValueError):
            budget(replace(p, output_bits=58))

    def test_ll13_precision_screen(self):
        result = estimate(PROFILES["ll13"])
        # TFHEpp ckks_tutti_ll13 measured 5.27 bits (largest slot error).
        self.assertAlmostEqual(result["expected_worst_slot_precision_bits"], 5.42, delta=0.05)
        self.assertAlmostEqual(result["rms_precision_bits"], 7.97, delta=0.05)
        self.assertEqual(result["accmul"]["selector_digits"], 4)
        self.assertFalse(result["security_validated"])

    def test_indicator_scale_is_cheaper_than_tree_scale(self):
        p = PROFILES["ll13"]
        flat = replace(p, indicator_scale_bits=p.scale_bits)
        self.assertEqual(p.selector_bits - flat.selector_bits, 2)
        self.assertGreater(estimate(p)["rms_precision_bits"],
                           estimate(flat)["rms_precision_bits"] + 0.5)

    def test_dense_weight_and_auxiliary_lift_matter(self):
        p = PROFILES["ll13"]
        self.assertLess(estimate(replace(p, dense_weight=0))["rms_precision_bits"],
                        estimate(p)["rms_precision_bits"])
        unlifted = estimate(replace(p, auxiliary_bits=0))
        self.assertLess(unlifted["rms_precision_bits"], 0)

    def test_validation(self):
        p = PROFILES["ll13"]
        for bad in (replace(p, stride=256), replace(p, output_bits=38),
                    replace(p, indicator_scale_bits=21), replace(p, message_bound=4.0),
                    replace(p, output_scale_bits=23)):
            with self.assertRaises(ValueError):
                bad.validate()
        TuttiParams().validate()

    def test_search_respects_budget_and_recovers_profile(self):
        p = PROFILES["ll13"]
        found = search(p, 212, gap_bits=(4,), s2c_bits=(14, 15), plain_bits=(17, 18, 19),
                       indicator_gain_bits=range(4))
        self.assertTrue(found)
        for candidate in found:
            self.assertLessEqual(candidate["selector_bits"], 212)
        best = found[0]["parameters"]
        self.assertEqual({k: best[k] for k in ("scale_bits", "indicator_scale_bits",
                                               "plain_bits", "s2c_bits", "output_bits")},
                         {"scale_bits": 22, "indicator_scale_bits": 24, "plain_bits": 18,
                          "s2c_bits": 14, "output_bits": 58})


if __name__ == "__main__":
    unittest.main()
