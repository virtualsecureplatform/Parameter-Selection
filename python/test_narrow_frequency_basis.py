import cmath
import itertools
import unittest
from dataclasses import replace

from noiseestimation.narrow_frequency_basis import minimum_positive_frequency_basis, frequency_basis_schedule_budget
from noiseestimation.narrow_ckks import NarrowParams, estimate


class FrequencyBasisTests(unittest.TestCase):
    def test_lean_basis_layout_and_normalization(self):
        p = NarrowParams(n=32768, selector_bits=518, scale_bits=47, plain_bits=39,
                         source_bits=48, input_scale_bits=44, s2c_bits=37,
                         primary_bits=18, bbar_bits=8, message_bound=4)
        budget = frequency_basis_schedule_budget(p, 40, coefficient_bits=31, auxiliary_bits=44)
        self.assertEqual((budget["extra_product_input_bits"], budget["output_bits"],
                          budget["output_scale_bits"]), (244, 197, 78))
        self.assertEqual(budget["retained_square"]["square_log_q"], 119)
        self.assertEqual(budget["normalization_scale_drop"], 34)
        self.assertTrue(budget["retained_square_can_normalize_for_refresh"])
        self.assertFalse(budget["retained_square_fits_refresh_phase_domain"])
        self.assertEqual(budget["maximum_evaluation_key_bits"], 523)
        self.assertEqual(budget["maximum_key_bits"], 572)
        self.assertEqual(budget["minimum_native_storage_bits_including_optional_packed_key"], 576)
        self.assertLess(budget["coefficient_quantization_uniform_bound_at_exact_roots"], 2**-23)
        model = estimate(p, streamed=True, relinearize_before_rescale=True,
                         auxiliary_switch_bits=44, selector_chunk_size=40, batch_columns=9)
        self.assertEqual(model["prepared_fft_memory"]["factor_automorphism_bytes"], 2283798528)
        self.assertEqual(model["prepared_fft_memory"]["dd_chain_bytes"], 2852126720)
        self.assertEqual(model["selector_memory"]["expanded_torus_bytes"], 179721732096)
        self.assertEqual(model["chunked_selector_memory"]["ntt_chunk_bytes"], 1384120320)

    def test_full_basis_memory_matches_compiled_cpp_payloads(self):
        p = NarrowParams(n=32768, selector_bits=536, scale_bits=49, plain_bits=39,
                         source_bits=48, input_scale_bits=44, s2c_bits=37)
        model = estimate(p, streamed=True, relinearize_before_rescale=True,
                         auxiliary_switch_bits=44, selector_chunk_size=40, batch_columns=9)
        self.assertEqual(model["prepared_fft_memory"]["factor_automorphism_bytes"], 9697230848)
        self.assertEqual(model["prepared_fft_memory"]["dd_chain_bytes"], 27099922432)
        self.assertEqual(model["chunked_selector_memory"]["ntt_chunk_bytes"], 1426063360)
        self.assertEqual(model["selector_memory"]["one_group_torus_bytes"], 5546967040)

    def test_basis_budget_matches_compiled_gaussian_fbs_profile(self):
        p = NarrowParams(n=512, selector_bits=440, scale_bits=40, plain_bits=34,
                         source_bits=48, input_scale_bits=44, s2c_bits=37)
        budget = frequency_basis_schedule_budget(p, 40, auxiliary_bits=16)
        self.assertEqual(budget["extra_product_input_bits"], 206)
        self.assertEqual((budget["output_bits"], budget["output_scale_bits"]), (166, 72))
        self.assertEqual(budget["retained_square"]["square_log_q"], 94)
        self.assertTrue(budget["retained_square"]["guarded_range_passed"])
        self.assertFalse(budget["retained_square_can_normalize_for_refresh"])
        self.assertLess(budget["coefficient_quantization_uniform_bound_at_exact_roots"], 2**-23)

    def test_full_basis_guard_and_repeatability_are_separate(self):
        p = NarrowParams(n=32768, selector_bits=504, scale_bits=49, plain_bits=39,
                         source_bits=48, input_scale_bits=44, s2c_bits=37)
        budget = frequency_basis_schedule_budget(p, 40, auxiliary_bits=44)
        self.assertFalse(budget["retained_square"]["guarded_range_passed"])
        self.assertEqual(budget["minimum_selector_bits_for_guarded_square"], 505)
        guarded = frequency_basis_schedule_budget(replace(p, selector_bits=505), 40, auxiliary_bits=44)
        self.assertTrue(guarded["retained_square"]["guarded_range_passed"])
        self.assertTrue(guarded["output_can_normalize_for_refresh"])
        self.assertFalse(guarded["retained_square_can_normalize_for_refresh"])
        repeatable = frequency_basis_schedule_budget(replace(p, selector_bits=536), 40, auxiliary_bits=44)
        self.assertTrue(repeatable["retained_square_can_normalize_for_refresh"])
        self.assertFalse(repeatable["retained_square_fits_refresh_phase_domain"])
        self.assertFalse(repeatable["retained_square_fits_declared_message_domain"])
        self.assertEqual(repeatable["maximum_key_bits"], 590)
        self.assertEqual(repeatable["minimum_native_storage_bits_including_optional_packed_key"], 640)
        self.assertFalse(repeatable["noise_validated"])
        self.assertFalse(repeatable["security_validated"])

    def test_invalid_budget_and_insufficient_basis_level(self):
        p = NarrowParams(n=512, selector_bits=440, scale_bits=40, plain_bits=34,
                         source_bits=48, input_scale_bits=44, s2c_bits=37)
        for bits in [0, 41, True, 32.0]:
            with self.assertRaises(ValueError):
                frequency_basis_schedule_budget(p, 12, coefficient_bits=bits)
        for auxiliary in [-1, True, 16.0]:
            with self.assertRaises(ValueError):
                frequency_basis_schedule_budget(p, 12, auxiliary_bits=auxiliary)
        direct = replace(p, selector_bits=387, scale_bits=49)
        with self.assertRaises(ValueError):
            frequency_basis_schedule_budget(direct, 12)

    def test_small_minima_against_unpruned_enumeration(self):
        for degree in range(1, 13):
            size = len(minimum_positive_frequency_basis(degree).basis)
            required = set(range(1, degree + 1))
            for candidate in itertools.combinations(range(1, degree + 1), size - 1):
                covered = set(candidate) | {a + b for a in candidate for b in candidate}
                self.assertFalse(required <= covered)

    def test_all_supported_degrees_cover_every_frequency(self):
        previous_size = 0
        for degree in range(1, 41):
            plan = minimum_positive_frequency_basis(degree)
            self.assertEqual(len(plan.reconstruction), degree)
            self.assertEqual(tuple(sorted(set(plan.basis))), plan.basis)
            self.assertGreaterEqual(len(plan.basis), previous_size)
            previous_size = len(plan.basis)
            for k, expression in enumerate(plan.reconstruction, 1):
                self.assertIn(len(expression), (1, 2))
                self.assertEqual(sum(expression), k)
                self.assertTrue(all(a in plan.basis for a in expression))

    def test_profile_sizes_and_post_reconstruction_complex_weights(self):
        for degree, expected_size in [(12, 4), (40, 9)]:
            plan = minimum_positive_frequency_basis(degree)
            self.assertEqual(len(plan.basis), expected_size)
            for phase in [-3.0, -0.2, 0.0, 0.7, 3.1]:
                z = cmath.exp(1j * phase)
                roots = {a: z**a for a in plan.basis}
                direct = reconstructed = 0.125
                for k, expression in enumerate(plan.reconstruction, 1):
                    coefficient = complex(0.3 / k, -0.2 / (k + 1))
                    root = roots[expression[0]]
                    if len(expression) == 2:
                        root *= roots[expression[1]]
                    weighted = coefficient * root
                    reconstructed += weighted + weighted.conjugate()
                    weighted_direct = coefficient * z**k
                    direct += weighted_direct + weighted_direct.conjugate()
                self.assertAlmostEqual(abs(direct - reconstructed), 0, places=13)

    def test_invalid_degrees(self):
        for degree in [0, 41, -1, 12.0, True, None]:
            with self.assertRaises(ValueError):
                minimum_positive_frequency_basis(degree)


if __name__ == "__main__":
    unittest.main()
