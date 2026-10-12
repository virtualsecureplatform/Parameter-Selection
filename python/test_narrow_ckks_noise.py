import math
import unittest
from dataclasses import replace
from noiseestimation.narrow_ckks import NarrowParams, estimate, switch_slot_variance, retained_square_range, minimum_square_output_bits, retained_unit_product_chain


class NarrowNoiseTests(unittest.TestCase):
    def test_ll13_fbs_requires_expanded_setup_and_product_normalization(self):
        p = NarrowParams(n=8192, selector_bits=217, scale_bits=27, plain_bits=24,
                         source_bits=25, input_scale_bits=21, s2c_bits=20, message_bound=4)
        p.validate()
        self.assertEqual(p.output_bits, 58)
        self.assertEqual(p.selector_bits - p.plain_bits + 24, 217)
        self.assertGreater(p.packed_bits, 218)
        self.assertFalse(retained_square_range(p.output_bits, p.scale_bits, 4.01)['guarded_range_passed'])
        for bound, drop, delta in [(1.01, 3, 24), (4.01, 6, 21)]:
            square = retained_square_range(p.output_bits - drop, p.scale_bits - drop, bound)
            self.assertTrue(square['guarded_range_passed'])
            self.assertEqual(square['square_log_q'], 31)
            self.assertEqual(square['square_log_delta'], delta)

    def test_ll14_fbs_plain37_budget_preserves_retained_square(self):
        p = NarrowParams(n=16384, selector_bits=385, scale_bits=48, plain_bits=37,
                         source_bits=48, input_scale_bits=44, s2c_bits=37, message_bound=4)
        p.validate()
        self.assertEqual(p.packed_bits, 436)
        self.assertEqual(p.output_bits, 108)
        model = estimate(p, streamed=True, relinearize_before_rescale=True,
                         auxiliary_switch_bits=28, selector_chunk_size=40, batch_columns=24)
        self.assertLessEqual(model["maximum_evaluation_key_bits"], 437)
        self.assertEqual(model["prepared_fft_memory"]["factor_automorphism_bytes"], 2316304384)
        self.assertEqual(model["prepared_fft_memory"]["dd_chain_bytes"], 5579735040)
        self.assertEqual(model["selector_memory"]["one_group_torus_bytes"], 970719232)
        square = retained_square_range(p.output_bits, p.scale_bits, 4.01)
        self.assertEqual(square["square_log_q"], 60)
        self.assertEqual(square["minimum_output_bits"], 106)
        self.assertTrue(square["guarded_range_passed"])

    def test_retained_nine_level_unit_chain(self):
        p = NarrowParams(n=16384, selector_bits=383, scale_bits=25, plain_bits=23,
                         source_bits=28, input_scale_bits=26, s2c_bits=23)
        chain = retained_unit_product_chain(p, 9, 23, 0.19, 58)
        self.assertEqual(chain["final_output_bits"], 26)
        self.assertEqual(chain["scale_reduction_bits"], 2)
        self.assertEqual(chain["minimum_refresh_output_bits"], 235)
        self.assertTrue(chain["guarded_range_passed"])
        self.assertLess(chain["eight_rms_plus_sine_screen"], 2**-4)
        self.assertFalse(retained_unit_product_chain(replace(p, selector_bits=382),
                         9, 23, 0.19, 58)["guarded_range_passed"])
        with self.assertRaises(ValueError):
            retained_unit_product_chain(p, 10, 23, 0.19, 58)
        with self.assertRaises(ValueError):
            retained_unit_product_chain(p, 9, 26, 0.19, 58)

    def test_file_backed_selector_memory_excludes_expanded_storage(self):
        p = NarrowParams(n=32768, selector_bits=669, scale_bits=46, plain_bits=39,
                         source_bits=50, input_scale_bits=43, s2c_bits=39)
        result = estimate(p, streamed=True, relinearize_before_rescale=True, auxiliary_switch_bits=212)
        memory = result["selector_memory"]
        self.assertEqual(memory["inferred_torus_storage_bits"], 896)
        self.assertGreater(memory["expanded_torus_bytes"], 250 * 2**30)
        self.assertLess(memory["file_backed_torus_and_ntt_peak_bytes"], 64 * 2**30)
        self.assertEqual(memory["file_backed_torus_and_ntt_peak_bytes"],
                         memory["one_group_torus_bytes"] + result["ntt_cache_peak_bytes"])
        self.assertGreater(result["raw_dd_chain_torus_bytes"], 600 * 2**30)
        self.assertLess(result["compact_dd_chain_digit_bytes"], 6 * 2**30)
        self.assertEqual(result["raw_dd_chain_torus_bytes"],
                         112 * result["compact_dd_chain_digit_bytes"])
        self.assertEqual(result["dense_s2c_plaintext_polynomial_bytes"], 112 * 2**30)
        self.assertEqual(result["compact_s2c_plaintext_coefficient_bytes"], 8 * 2**30)
        fft = result["prepared_fft_memory"]
        self.assertEqual(fft["factor_automorphism_bytes"], 23341826048)
        self.assertEqual(fft["output_automorphism_bytes"], 12273057792)
        self.assertEqual(fft["dd_chain_bytes"], 48327819264)
        self.assertEqual(fft["tail_rotation_binary_keys"], 9)
        self.assertEqual(fft["tail_rotation_bytes"], 210076434432)
        self.assertEqual(fft["bounded_tail_rotation_peak_bytes"], fft["factor_automorphism_bytes"])
        self.assertEqual(fft["staged_factor_phase_bytes"], 46683652096)
        self.assertEqual(fft["staged_product_tree_phase_bytes"], 48327819264)
        self.assertEqual(fft["staged_output_phase_bytes"], 12273057792)
        self.assertEqual(fft["staged_peak_bytes"], 48327819264)

    def test_ll14_fft_storage_matches_cpp_layout(self):
        p = NarrowParams(n=16384, selector_bits=387, scale_bits=49, plain_bits=34,
                         source_bits=48, input_scale_bits=44, s2c_bits=37)
        fft = estimate(p, streamed=True, relinearize_before_rescale=True,
                       auxiliary_switch_bits=28)["prepared_fft_memory"]
        self.assertEqual(fft["tail_rotation_binary_keys"], 6)
        self.assertEqual(fft["bounded_tail_rotation_peak_bytes"], 2415919104)
        self.assertEqual(fft["tail_rotation_bytes"], 14495514624)
        self.assertEqual(fft["dd_chain_bytes"], 5792595968)

    def test_chunked_ntt_storage_preserves_noise_screen(self):
        p = NarrowParams(n=32768, selector_bits=669, scale_bits=46, plain_bits=39,
                         source_bits=50, input_scale_bits=43, s2c_bits=39)
        full = estimate(p, streamed=True, relinearize_before_rescale=True, auxiliary_switch_bits=212)
        chunked = estimate(p, streamed=True, relinearize_before_rescale=True,
                           auxiliary_switch_bits=212, selector_chunk_size=40, batch_columns=2)
        memory = chunked["chunked_selector_memory"]
        self.assertEqual(memory["ntt_chunk_bytes"] * 1058, full["ntt_cache_peak_bytes"] * 40)
        self.assertEqual(memory["modular_sum_bytes"], 176160768)
        self.assertLess(memory["file_backed_group_chunk_and_sums_bytes"], 10 * 2**30)
        self.assertEqual(chunked["eight_rms_plus_sine_screen"], full["eight_rms_plus_sine_screen"])
        for size, columns in [(-1, 1), (1.5, 1), (1, 0)]:
            with self.assertRaises(ValueError):
                estimate(p, selector_chunk_size=size, batch_columns=columns)

    def test_retained_square_range_includes_constant_complex_slots(self):
        old = retained_square_range(100, 49, 4)
        self.assertFalse(old["amplitude_fits"])
        self.assertEqual(old["centered_coefficient_capacity_log2"], 1)
        self.assertEqual(old["square_coefficient_bound_log2"], 5)
        self.assertEqual(old["minimum_output_bits"], 108)
        self.assertTrue(retained_square_range(108, 49, 4)["guarded_range_passed"])
        self.assertFalse(retained_square_range(104, 49, 4, 0)["amplitude_fits"])
        self.assertTrue(retained_square_range(105, 49, 4, 0)["amplitude_fits"])
        self.assertTrue(retained_square_range(100, 49, 0.125)["guarded_range_passed"])
        for bound in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                minimum_square_output_bits(49, bound)
        with self.assertRaises(ValueError):
            minimum_square_output_bits(49, 4, -1)

    def test_auxiliary_scale_preserves_fresh_noise_and_counts_key_moduli(self):
        p = NarrowParams(n=8192, selector_bits=187, scale_bits=24, plain_bits=17,
                         source_bits=25, input_scale_bits=24, s2c_bits=18)
        ordinary = estimate(p, relinearize_before_rescale=True)
        raised = estimate(p, relinearize_before_rescale=True, auxiliary_switch_bits=28)
        self.assertEqual(raised["s2c_log2_variances"]["fresh"],
                         ordinary["s2c_log2_variances"]["fresh"])
        self.assertEqual(raised["structured_switch_key_bits"], 53)
        self.assertEqual(raised["maximum_evaluation_key_bits"], 217)
        self.assertGreater(estimate(p, auxiliary_switch_bits=128)["maximum_evaluation_key_bits"], 217)
        self.assertLess(raised["eight_rms_plus_sine_screen"], 2**-4)
        # Even noiseless raised keys incur rounding when divided back.
        noiseless = replace(p, sigma=0)
        self.assertEqual(switch_slot_variance(noiseless, 170, 4), float("-inf"))
        self.assertTrue(math.isfinite(switch_slot_variance(noiseless, 170, 4, 28)))
        with self.assertRaises(ValueError):
            estimate(p, auxiliary_switch_bits=-1)

    def test_raw_switch_has_two_ring_gains(self):
        p = NarrowParams(n=512)
        # 41 rows, balanced base-16 second moment 21.5, sigma=3.2.
        expected = 41 * 512**2 * 21.5 * 3.2**2 / 2**40
        self.assertAlmostEqual(2**switch_slot_variance(p, 164, 4), expected)
        self.assertAlmostEqual(switch_slot_variance(replace(p, n=1024), 164, 4) -
                               switch_slot_variance(p, 164, 4), 2)

    def test_rotation_chain_divides_once(self):
        p = NarrowParams(n=8192, scale_bits=24)
        single = 2**switch_slot_variance(p, 170, 4, 28)
        chain = 2**switch_slot_variance(p, 170, 4, 28, 6)
        rounding = 2**switch_slot_variance(replace(p, sigma=0), 170, 4, 28)
        self.assertAlmostEqual(chain - rounding, 6 * (single - rounding), delta=1e-20)
        self.assertEqual(switch_slot_variance(p, 170, 4, 28, 0), float("-inf"))
        self.assertEqual(estimate(NarrowParams(n=8192))["maximum_tail_rotation_switches"], 6)
        with self.assertRaises(ValueError):
            switch_slot_variance(p, 170, 4, 28, -1)

    def test_zero_noise_retains_tensor_rounding(self):
        p = NarrowParams(n=512, sigma=0)
        result = estimate(p)
        self.assertGreater(result["complex_factor_tree_rms"], 0)
        self.assertTrue(math.isfinite(result["levels"][0]["tensor_rounding_log2_variance"]))
        self.assertIsNone(result["levels"][0]["relin_log2_variance"])
        self.assertGreater(estimate(replace(p, sigma=3.2))["complex_factor_tree_rms"],
                           result["complex_factor_tree_rms"])

    def test_schedule_and_cache_match_small_executable(self):
        result = estimate(NarrowParams(n=512))
        self.assertEqual(result["output_bits"], 64)
        self.assertEqual(result["packed_bits"], 213)
        self.assertEqual(result["ntt_cache_bytes_upper_bound"], 58490880)
        with self.assertRaises(ValueError):
            estimate(NarrowParams(selector_bits=100))

    def test_complex_sine_bound_and_streamed_memory(self):
        p = NarrowParams(n=16384, source_bits=48, input_scale_bits=41, s2c_bits=37,
                         selector_bits=379, scale_bits=49, plain_bits=34)
        result = estimate(p, streamed=True)
        bound = result["real_sine_error_bound"]
        gamma = 128
        for i in range(1001):
            x = p.message_bound * i / 1000
            error = abs(x - gamma / (2 * math.pi) * math.sin(2 * math.pi * x / gamma))
            self.assertLessEqual(error, bound)
        self.assertAlmostEqual(result["complex_sine_error_bound"], math.sqrt(2) * bound)
        self.assertEqual(result["ntt_cache_peak_bytes"], 6656360448)
        # Counting the previously omitted tail rotation chains changes the
        # old schedule's screen; preserve the 16-bit gate rather than relax it.
        self.assertGreater(result["eight_rms_plus_sine_screen"], 2**-16)
        raised = estimate(p, streamed=True, relinearize_before_rescale=True,
                          auxiliary_switch_bits=28)
        self.assertLess(raised["eight_rms_plus_sine_screen"], 2**-16)


if __name__ == "__main__":
    unittest.main()
