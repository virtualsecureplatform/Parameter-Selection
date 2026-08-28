#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compact_cover_schedule import (
    general_baby_giant,
    index_to_sequence,
    scalar_direct_gate_manifest_json,
    scalar_direct_gate_manifest_sha256,
    thin_power_of_two_gate_manifest,
    thin_power_of_two_schedule,
)


class CompactCoverScheduleTest(unittest.TestCase):
    def test_mixed_radix(self) -> None:
        self.assertEqual(index_to_sequence(1, (2, 3)), (0, 0))
        self.assertEqual(index_to_sequence(6, (2, 3)), (1, 2))

    def test_general_baby_giant_matches_magma_heuristic(self) -> None:
        self.assertEqual(general_baby_giant((2, 16_384)), ((2, 91), (1, 181)))

    def test_target_schedule(self) -> None:
        schedule = thin_power_of_two_schedule()
        self.assertEqual(schedule.frobenius_order, 2)
        self.assertEqual(schedule.slot_count, 32_768)
        self.assertEqual(schedule.dimension_sizes, (16_384, 2))
        self.assertEqual(schedule.baby_product, 182)
        self.assertEqual(schedule.giant_product, 181)
        self.assertEqual(len(schedule.switch_key_exponents), 362)
        self.assertEqual(schedule.generated_group_size, 65_536)
        self.assertEqual(schedule.peak_live_ciphertexts, 368)

    def test_canonical_gate_manifest(self) -> None:
        manifest = thin_power_of_two_gate_manifest()
        self.assertEqual(len(manifest.baby_switches), 181)
        self.assertEqual(len(manifest.giant_switches), 180)
        self.assertEqual(manifest.backward_exponent, 65_537)
        self.assertEqual(manifest.trace_exponent, 65_537)
        self.assertEqual(len(manifest.stages), 6)
        used = {
            use.exponent
            for use in manifest.baby_switches + manifest.giant_switches
        }
        used.update((manifest.backward_exponent, manifest.trace_exponent))
        self.assertEqual(used, manifest.schedule.switch_key_exponents)
        self.assertEqual(len(manifest.canonical_json()), 29_433)
        self.assertEqual(
            manifest.sha256(),
            "6d376c090ad655badcb23b00da1e16c4429c83d3f20c3195df93eef5bf94049e",
        )

    def test_scalar_direct_manifest(self) -> None:
        self.assertEqual(len(scalar_direct_gate_manifest_json()), 640)
        self.assertEqual(
            scalar_direct_gate_manifest_sha256(),
            "209b8826908383c92bce2ea41f27eda9febbf250e69e2c5541f5c18b76e454f0",
        )


if __name__ == "__main__":
    unittest.main()
