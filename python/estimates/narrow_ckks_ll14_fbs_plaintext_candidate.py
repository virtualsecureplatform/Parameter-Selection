#!/usr/bin/env python3
"""LL14 FBS candidate improving table encoding precision within log(PQ)=437.

The factor/tree model excludes the trigonometric coefficients and FFT error;
its RMS is not a complete nonlinear-function accuracy certificate.
"""
import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from noiseestimation.narrow_ckks import NarrowParams, estimate, retained_square_range


def validation_report(path, auxiliary_bits=28):
    """Independently check a full-ring report; this does not certify security."""
    result = json.loads(path.read_text())
    name = result.get('function')
    if name not in {'sigmoid', 'gelu', 'tanh'}:
        raise ValueError('unsupported validation function: ' + str(name))
    if auxiliary_bits not in (16, 28):
        raise ValueError('unsupported LL14 plain37 auxiliary scale')
    profile = 'll14-fbs-plain37-aux16-dd' if auxiliary_bits == 16 else 'll14-fbs-plain37-dd'
    expected = dict(profile=profile, n=16384, h=31,
                    source_log_q=48, selector_log_q=385, packed_log_q=436,
                    input_log_q=85, input_log_delta=44, output_log_q=108,
                    output_log_delta=48, root_plain_log_delta=37,
                    s2c_plain_log_delta=37, storage_bits=448, gaussian_sigma=3.2,
                    gadget_bits=4, secondary_digit_bits=4, auxiliary_switch_bits=auxiliary_bits,
                    factor_key_log_q=348+auxiliary_bits, structured_switch_key_log_q=48+auxiliary_bits,
                    s2c_input_key_log_q=85+auxiliary_bits, message_bound=4,
                    input_pattern='random-with-endpoints', series_tile_size=12,
                    degree=40 if name == 'tanh' else 12,
                    approximation_method='windowed-remez' if name == 'tanh' else 'direct-remez',
                    packed_selectors=False, file_backed_selectors=True,
                    compact_dd=True, compact_s2c=True, lean_dd=False,
                    cache_mode='streamed', bounded_rotation_fft=True,
                    staged_fft_caches=True, selector_ntt_chunk_size=40,
                    tree_depth=5, retained_unit_product_levels=0,
                    relinearize_before_rescale=True,
                    post_refresh_product_relinearize_before_rescale=True)
    for key, value in expected.items():
        actual = result.get(key)
        if type(actual) is not type(value) and not (type(value) in (int, float) and
                                                  type(actual) in (int, float)):
            raise ValueError('validation profile type mismatch: ' + key)
        if actual != value:
            raise ValueError('validation profile mismatch: ' + key)
    if type(result.get('accuracy_passed')) is not bool:
        raise ValueError('accuracy_passed must be boolean')
    tolerance = 2**-23
    for key in ('evaluation_tolerance', 'target_error_tolerance'):
        value = result.get(key)
        if type(value) not in (int, float) or not math.isclose(value, tolerance,
                                                            rel_tol=1e-14, abs_tol=0):
            raise ValueError('validation uses a different precision gate: ' + key)
    fields = ('max_error', 'evaluation_max_error', 'fresh_input_max_error',
              'approximation_uniform_bound', 'post_refresh_product_max_error',
              'post_refresh_product_evaluation_max_error', 'rms_error')
    for key in fields:
        value = result.get(key)
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError('invalid validation error measurement: ' + key)
    endpoint = (1/(1+math.exp(-4)) if name == 'sigmoid' else
                2*(1+math.erf(4/math.sqrt(2))) if name == 'gelu' else math.tanh(4))
    amplitude = math.sqrt(2)*endpoint
    uniform = math.sqrt(2)*result['approximation_uniform_bound']
    product_tolerance = 2*amplitude*tolerance + tolerance*tolerance + tolerance
    evaluation_product_tolerance = (2*max(0, amplitude-uniform)*tolerance +
                                    tolerance*tolerance + tolerance)
    checked = (all(result[key] < tolerance for key in fields[:3]) and
               uniform < tolerance and result['rms_error'] < tolerance and
               result['post_refresh_product_max_error'] < product_tolerance and
               result['post_refresh_product_evaluation_max_error'] < evaluation_product_tolerance)
    if result['accuracy_passed'] and not checked:
        raise ValueError('reported success fails independent accuracy gates')
    return name, dict(report=str(path), target_precision_bits=23,
                      accuracy_passed=result['accuracy_passed'] and checked,
                      **{key: result[key] for key in fields},
                      target_tolerance=tolerance, product_tolerance=product_tolerance,
                      product_evaluation_tolerance_lower_bound=evaluation_product_tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation-report', type=Path, action='append', default=[],
                        help='check a completed full-ring report; repeat for each function')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--auxiliary-switch-bits', type=int, choices=(16, 28), default=28)
    args = parser.parse_args()
    if args.output is None:
        source = Path(__file__)
        args.output = (source.with_name(source.stem + '_aux16.json')
                       if args.auxiliary_switch_bits == 16 else source.with_suffix('.json'))
    params = NarrowParams(n=16384, selector_bits=385, scale_bits=48,
                         plain_bits=37, source_bits=48, input_scale_bits=44,
                         s2c_bits=37, message_bound=4)
    model = estimate(params, streamed=True, relinearize_before_rescale=True,
                     auxiliary_switch_bits=args.auxiliary_switch_bits, selector_chunk_size=40,
                     batch_columns=12)
    square = retained_square_range(params.output_bits, params.scale_bits, 4.01)
    budget_passed = model["maximum_evaluation_key_bits"] <= 437 and square["guarded_range_passed"]
    report = {"scope": __doc__, "profile": 'll14-fbs-plain37-aux16-dd' if args.auxiliary_switch_bits == 16 else 'll14-fbs-plain37-dd',
              "parameters": asdict(params), "auxiliary_switch_bits": args.auxiliary_switch_bits,
              "maximum_key_modulus_budget": 437, "target_precision_bits": 23,
              "output_bits": params.output_bits, "retained_square": square,
              "budget_passed": budget_passed, "factor_tree_model": model,
              "encrypted_validation_passed": False, "security_validated": False}
    validated = {}
    try:
        for path in args.validation_report:
            name, checked = validation_report(path, args.auxiliary_switch_bits)
            if name in validated:
                raise ValueError('duplicate validation function: ' + name)
            validated[name] = checked
    except (ValueError, KeyError, OSError) as error:
        parser.error(str(error))
    report['full_encrypted_validation'] = validated
    report['missing_full_validation_functions'] = sorted({'sigmoid', 'gelu', 'tanh'} - set(validated))
    report['encrypted_validation_passed'] = (not report['missing_full_validation_functions'] and
                                            all(v['accuracy_passed'] for v in validated.values()))
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print("budget", budget_passed, "packed bits", params.packed_bits,
          "output bits", params.output_bits, "root RMS", model["complex_factor_tree_rms"])
    if validated:
        print('checked accuracy', {name: v['accuracy_passed'] for name, v in validated.items()},
              'missing', report['missing_full_validation_functions'],
              'all functions passed', report['encrypted_validation_passed'])
    return 0 if budget_passed and all(v['accuracy_passed'] for v in validated.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
