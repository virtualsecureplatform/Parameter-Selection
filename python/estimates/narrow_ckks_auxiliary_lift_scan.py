#!/usr/bin/env python3
"""Screen LL14 DD auxiliary lifts using the factor/tree RMS model.

This excludes nonlinear coefficient weights, maximum-slot tails and FFT error.
It neither validates encrypted accuracy nor qualifies structured-key security.
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from noiseestimation import narrow_ckks


def scan():
    params = narrow_ckks.NarrowParams(
        n=16384, selector_bits=385, scale_bits=48, plain_bits=37,
        source_bits=48, input_scale_bits=44, s2c_bits=37, message_bound=4)
    rows = []
    input_bits = params.source_bits + params.s2c_bits
    factor_bits = params.selector_bits - params.plain_bits
    for auxiliary in range(0, 33, 4):
        model = narrow_ckks.estimate(
            params, streamed=True, relinearize_before_rescale=True,
            auxiliary_switch_bits=auxiliary, selector_chunk_size=40,
            batch_columns=12)
        rms = model['complex_factor_tree_rms']
        if not math.isfinite(rms) or rms <= 0:
            raise ValueError('invalid model RMS')
        rows.append(dict(
            auxiliary_bits=auxiliary,
            structured_key_bits=params.source_bits + auxiliary,
            factor_key_bits=factor_bits + auxiliary,
            s2c_input_key_bits=input_bits + auxiliary,
            output_conjugation_key_bits=params.output_bits + auxiliary,
            root_model_rms=rms,
            root_model_below_23bit_tolerance=rms < 2**-23,
            model_maximum_key_bits=model['maximum_evaluation_key_bits'],
            packed_budget_used_by_full_uncompressed_run=False,
            actual_uncompressed_max_dense_modulus_bits=max(
                params.selector_bits, factor_bits + auxiliary,
                input_bits + auxiliary, params.output_bits + auxiliary)))
    model_source = Path(narrow_ckks.__file__)
    return dict(
        utc=datetime.now(timezone.utc).isoformat(), scope=__doc__,
        parameters=asdict(params),
        model_source_sha256=hashlib.sha256(model_source.read_bytes()).hexdigest(),
        generator_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        rows=rows,
        auxiliary8_to16_root_rms_ratio=rows[2]['root_model_rms']/rows[4]['root_model_rms'],
        encrypted_accuracy_established_by_model=False, security_validated=False,
        limitations=[
            'Root RMS excludes nonlinear coefficient weights, max-slot tails and FFT error',
            'Unused packed modulus436 remains a conservative dense proxy, not an exposed key in the uncompressed run',
            'C++ implementation and encrypted validation are assessed separately from this model'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name(
        'narrow_ckks_ll14_auxiliary_lift_scan.json'))
    args = parser.parse_args()
    report = scan()
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    for row in report['rows']:
        print('auxiliary', row['auxiliary_bits'], 'root RMS', row['root_model_rms'],
              'dense modulus', row['actual_uncompressed_max_dense_modulus_bits'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
