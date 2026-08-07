#!/usr/bin/env python3

import math
import unittest

from python.noiseestimation.clpx import (
    _failure_log2,
    estimate_clpx_to_tlwes,
    estimate_clpx_multiplication_depth,
    estimate_tlwes_to_clpx,
    paper_clpx_fixed_noise_threshold,
    pbs_input_margin_log2,
)
from python.noiseestimation.params import clpx as params


class CLPXNoiseTest(unittest.TestCase):
    def switch_estimate(self, w: int):
        return estimate_tlwes_to_clpx(
            validbit=128,
            num_multi=4,
            shift=w // 4 - 1,
            w=w,
            iksP=params.lvl1hparam,
            bkP=params.SS2CLPXlvlh2param,
            sskP=params.SS2CLPXlvl22param,
        )

    def multiplication_estimate(self, w: int):
        switch = self.switch_estimate(w)
        return estimate_clpx_multiplication_depth(
            initial_coefficient_variance=switch.packed_variance,
            validbit=128,
            num_multi=4,
            shift=w // 4 - 1,
            w=w,
            max_multiplications=1,
            P=params.SS2CLPXlvl2param,
            iksP=params.lvl1hparam,
            bkP=params.SS2CLPXlvlh2param,
            sskP=params.SS2CLPXlvl22param,
            model="paper",
        )

    def test_delta_b_approximation_variance(self):
        estimate = self.switch_estimate(20)
        self.assertEqual(math.log2(estimate.approximation_variance), 88.0)

    def test_paper_w_transition(self):
        self.assertEqual(self.multiplication_estimate(16).steps[1].status, "FAIL")
        self.assertEqual(self.multiplication_estimate(20).steps[1].status, "OK")
        self.assertLess(
            self.multiplication_estimate(20).steps[1].log2_failure, -64.0
        )

    def test_relinearization_width_is_not_delta_b_width(self):
        expected = paper_clpx_fixed_noise_threshold(
            P=params.SS2CLPXlvl2param,
            relinearization_base=params.SS2CLPXlvl2param.B,
        )
        self.assertEqual(
            paper_clpx_fixed_noise_threshold(P=params.SS2CLPXlvl2param),
            expected,
        )

    def test_optimized_scheme_switch_decomposition(self):
        self.assertEqual(params.SS2CLPXlvl2param.l, 3)
        self.assertEqual(params.SS2CLPXlvl2param.Bbit, 13)

        class OriginalSS2CLPXParam(params.SS2CLPXlvl2param):
            l = 4
            lₐ = l
            ℬbit = 9
            ℬₐbit = ℬbit
            ℬ = 2**ℬbit
            ℬₐ = ℬ

        class OriginalBkParam:
            domainP = params.lvlhalfparam
            targetP = OriginalSS2CLPXParam

        class OriginalSSKParam:
            t = params.lvl22param.t
            basebit = params.lvl22param.basebit
            domainP = OriginalSS2CLPXParam
            targetP = OriginalSS2CLPXParam

        original_switch = estimate_tlwes_to_clpx(
            validbit=128,
            num_multi=4,
            shift=4,
            w=20,
            iksP=params.lvl1hparam,
            bkP=OriginalBkParam,
            sskP=OriginalSSKParam,
        )
        original_mult = estimate_clpx_multiplication_depth(
            initial_coefficient_variance=original_switch.packed_variance,
            validbit=128,
            max_multiplications=1,
            P=OriginalSS2CLPXParam,
            model="paper",
        )
        optimized_log_failure = (
            self.multiplication_estimate(20).steps[1].log2_failure
        )
        original_log_failure = original_mult.steps[1].log2_failure
        self.assertLess(optimized_log_failure, -64.0)
        self.assertLess(abs(optimized_log_failure - original_log_failure), 0.2)

    def test_many_lut_uses_semantic_input_interval(self):
        estimate = self.switch_estimate(20)
        margin = pbs_input_margin_log2(
            params.SS2CLPXlvlh2param,
            estimate.iks_variance,
            num_out=4,
            input_plain_modulus=params.lvl1param.plain_modulus,
        )
        self.assertGreater(margin, 12.0)

    def test_failure_probability_stays_in_log_domain(self):
        log_probability = _failure_log2(1.0, 1e-6)
        self.assertTrue(math.isfinite(log_probability))
        self.assertLess(log_probability, -1000.0)

    def test_clpx_to_tfhe_block_matches_data_digits(self):
        estimate = estimate_clpx_to_tlwes(validbit=13, numdigit=9, basebit=2)
        self.assertEqual(estimate.produced_tlwes, 13)
        with self.assertRaises(ValueError):
            estimate_clpx_to_tlwes(
                validbit=13, batch_size=8, numdigit=9, basebit=2
            )


if __name__ == "__main__":
    unittest.main()
