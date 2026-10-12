#!/usr/bin/env python3
"""LL13 expanded-selector FBS candidate; no packed setup within key budget.

The unweighted factor/tree RMS excludes function coefficients and FFT error.
Published precision targets are experimental gates, not predicted accuracy.
"""
import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from noiseestimation.narrow_ckks import NarrowParams, estimate, retained_square_range


def memory_budget(model, compiled):
    """Combine compiled payloads with disjoint workspaces; not measured RSS."""
    if (compiled['n'], compiled['selector_bits'], compiled['root_scale_bits'],
        compiled['root_plain_bits'], compiled['primary_bits'], compiled['secondary_bits']) != (8192, 217, 27, 24, 4, 4):
        raise ValueError('compiled storage belongs to a different LL13 schedule')
    fft = model['prepared_fft_memory']
    for key, expected in [('prepared_factor_auto_bytes', fft['factor_automorphism_bytes']),
                          ('prepared_dd_chain_bytes', fft['dd_chain_bytes']),
                          ('one_selector_group_bytes', model['selector_memory']['one_group_torus_bytes'])]:
        if compiled[key] != expected:
            raise ValueError('compiled storage and independent model disagree: ' + key)
    columns, n, torus_bytes = 12, 8192, 32
    ciphertext_bytes = 2 * n * torus_bytes
    persistent = compiled['compact_dd_setup_bytes'] + compiled['compact_s2c_term_bytes']
    digit_workspace = ((217 + 15) // 16) * columns * 2 * ciphertext_bytes
    factor_capacity = columns * 32 * 2 * ciphertext_bytes
    chunk = model['chunked_selector_memory']
    factor = (persistent + 2 * compiled['prepared_factor_auto_bytes'] +
              compiled['one_selector_group_bytes'] + chunk['ntt_chunk_bytes'] +
              chunk['modular_sum_bytes'] + digit_workspace + factor_capacity)
    product = persistent + compiled['prepared_dd_chain_bytes'] + factor_capacity
    extra, permutation, reserve = 2 * columns * ciphertext_bytes, 64 * 2**20, 10 * 2**30
    envelope = max(factor, product) + extra + permutation + reserve
    limit = 16 * 2**30
    if envelope >= limit:
        raise ValueError('compiled memory envelope exceeds proposed limit')
    return dict(scope='Conservative compiled payloads plus disjoint workspaces and runtime reserve; not measured RSS',
                profile=compiled['profile'], n=n, batch_columns=columns, selector_chunk_size=40,
                persistent_setup_and_compact_s2c_bytes=persistent,
                conservative_factor_phase_payload_bytes=factor, product_phase_payload_bytes=product,
                digit_accumulation_workspace_bytes=digit_workspace,
                head_tail_factor_capacity_bytes=factor_capacity, additional_factor_accum_bytes=extra,
                permutation_bytes=permutation, temporary_and_runtime_reserve_bytes=reserve,
                estimated_process_envelope_bytes=envelope, proposed_virtual_address_limit_bytes=limit,
                file_selector_payload_bytes=compiled['expanded_selector_bytes'],
                full_size_job_launched=False, full_size_validation_passed=False)


def validation_report(path):
    """Check one full-ring result against this candidate's fixed accuracy gates."""
    result = json.loads(path.read_text())
    name = result.get('function')
    precisions = {'sigmoid': 9.04, 'gelu': 7.17, 'tanh': 7.89}
    if name not in precisions:
        raise ValueError('unsupported validation function: ' + str(name))
    drop = 6 if name == 'gelu' else 3
    expected = dict(profile='ll13-fbs-key217-dd', n=8192, h=31, source_log_q=25,
                    selector_log_q=217, packed_log_q=254, input_log_q=45,
                    input_log_delta=21, output_log_q=58, output_log_delta=27,
                    root_plain_log_delta=24, s2c_plain_log_delta=20, storage_bits=256,
                    gaussian_sigma=3.2, gadget_bits=4, secondary_digit_bits=4,
                    auxiliary_switch_bits=24, factor_key_log_q=217,
                    structured_switch_key_log_q=49, s2c_input_key_log_q=69,
                    message_bound=4, input_pattern='random-with-endpoints',
                    degree=40 if name == 'tanh' else 12,
                    approximation_method='windowed-remez' if name == 'tanh' else 'direct-remez',
                    packed_selectors=False, file_backed_selectors=True,
                    compact_dd=True, compact_s2c=True, lean_dd=False, cache_mode='streamed',
                    bounded_rotation_fft=True, staged_fft_caches=True, selector_ntt_chunk_size=40,
                    maximum_used_evaluation_key_bits=217, optional_packed_setup_within_key_budget=False,
                    retained_unit_product_levels=0, relinearize_before_rescale=True,
                    post_refresh_product_relinearize_before_rescale=True,
                    retained_square_normalization_bits=drop, retained_square_log_q=31,
                    retained_square_log_delta=27-drop)
    for key, value in expected.items():
        if result.get(key) != value or (isinstance(value, bool) and result.get(key) is not value):
            raise ValueError('validation profile mismatch: ' + key)
    tolerance = 2**(-precisions[name])
    for key in ['evaluation_tolerance', 'target_error_tolerance']:
        if not math.isclose(result[key], tolerance, rel_tol=1e-14, abs_tol=0):
            raise ValueError('validation uses a different precision gate: ' + key)
    fields = ['max_error', 'evaluation_max_error', 'fresh_input_max_error',
              'approximation_uniform_bound', 'post_refresh_product_max_error',
              'post_refresh_product_evaluation_max_error']
    if any(type(result[k]) not in (float, int) or not math.isfinite(result[k]) or
           result[k] < 0 for k in fields):
        raise ValueError('invalid validation error measurement')
    endpoint = (1/(1+math.exp(-4)) if name == 'sigmoid' else
                2*(1+math.erf(4/math.sqrt(2))) if name == 'gelu' else math.tanh(4))
    amplitude = math.sqrt(2)*endpoint
    approximation_bound = math.sqrt(2)*result['approximation_uniform_bound']
    square_tolerance = 2*amplitude*tolerance + tolerance*tolerance + tolerance
    # The approximant's measured maximum is at least its corner amplitude
    # minus the uniform bound. Use that lower bound for a stricter check.
    evaluation_square_tolerance = (2*max(0, amplitude-approximation_bound)*tolerance +
                                   tolerance*tolerance + tolerance)
    checked = (all(result[k] < tolerance for k in fields[:3]) and
               approximation_bound < tolerance and
               result['post_refresh_product_max_error'] < square_tolerance and
               result['post_refresh_product_evaluation_max_error'] < evaluation_square_tolerance)
    if result.get('accuracy_passed') is True and not checked:
        raise ValueError('reported success fails independent accuracy gates')
    return name, dict(report=str(path), n=8192, target_precision_bits=precisions[name],
                      accuracy_passed=result.get('accuracy_passed') is True and checked,
                      **{k: result[k] for k in fields},
                      target_tolerance=tolerance, product_tolerance=square_tolerance,
                      product_evaluation_tolerance_lower_bound=evaluation_square_tolerance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--storage', type=Path, help='TFHEpp compiled --storage JSON report')
    parser.add_argument('--memory-report', type=Path, help='write the full process envelope JSON')
    parser.add_argument('--validation-report', type=Path, action='append', default=[],
                        help='check a completed full-ring function report; repeat for all functions')
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    if bool(args.storage) != bool(args.memory_report):
        parser.error('--storage and --memory-report must be supplied together')
    params = NarrowParams(n=8192, selector_bits=217, scale_bits=27,
                         plain_bits=24, source_bits=25, input_scale_bits=21,
                         s2c_bits=20, message_bound=4)
    model = estimate(params, streamed=True, relinearize_before_rescale=True,
                     auxiliary_switch_bits=24, selector_chunk_size=40,
                     batch_columns=12)
    actual_key_bits = max(params.selector_bits,
                          params.selector_bits - params.plain_bits + 24,
                          params.source_bits + params.s2c_bits + 24,
                          params.output_bits + 24)
    functions = {}
    for name, precision, bound, drop in [('sigmoid', 9.04, 1.01, 3),
                                         ('gelu', 7.17, 4.01, 6),
                                         ('tanh', 7.89, 1.01, 3)]:
        square = retained_square_range(params.output_bits - drop,
                                       params.scale_bits - drop, bound)
        functions[name] = dict(target_precision_bits=precision,
                               target_tolerance=2**(-precision),
                               normalization_bits=drop, retained_square=square)
    passed = actual_key_bits <= 218 and all(
        f['retained_square']['guarded_range_passed'] for f in functions.values())
    report = dict(scope=__doc__, profile='ll13-fbs-key217-dd',
                  parameters=asdict(params), auxiliary_switch_bits=24,
                  maximum_used_evaluation_key_bits=actual_key_bits,
                  maximum_key_modulus_budget=218,
                  optional_packed_setup_bits=params.packed_bits,
                  optional_packed_setup_within_key_budget=params.packed_bits <= 218,
                  factor_tree_model=model, functions=functions,
                  budget_passed=passed, encrypted_validation_passed=False,
                  security_validated=False)
    validated = {}
    try:
        for path in args.validation_report:
            name, checked = validation_report(path)
            if name in validated:
                raise ValueError('duplicate validation function: ' + name)
            validated[name] = checked
    except (ValueError, KeyError, OSError) as error:
        parser.error(str(error))
    report['full_encrypted_validation'] = validated
    report['missing_full_validation_functions'] = sorted(set(functions)-set(validated))
    report['encrypted_validation_passed'] = (set(validated) == set(functions) and
                                            all(v['accuracy_passed'] for v in validated.values()))
    if args.storage:
        profiles = json.loads(args.storage.read_text())['profiles']
        compiled = next(p for p in profiles if p['profile'] == report['profile'])
        args.memory_report.write_text(json.dumps(memory_budget(model, compiled), indent=2) + '\n')
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print('expanded budget', passed, 'used key bits', actual_key_bits,
          'optional packed bits', params.packed_bits, 'root RMS', model['complex_factor_tree_rms'])
    if validated:
        print('checked accuracy', {name: v['accuracy_passed'] for name, v in validated.items()},
              'missing', report['missing_full_validation_functions'],
              'all functions passed', report['encrypted_validation_passed'])
    return 0 if passed and all(v['accuracy_passed'] for v in validated.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
