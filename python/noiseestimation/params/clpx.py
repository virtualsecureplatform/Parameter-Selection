#!/usr/bin/env python3
"""TFHEpp parameter presets used by the CLPX scheme-switch estimator."""

from __future__ import annotations


def _secret_stats(key_min: int, key_max: int) -> tuple[float, float]:
    values = list(range(key_min, key_max + 1))
    mean = sum(values) / len(values)
    second = sum(v * v for v in values) / len(values)
    return mean, second - mean * mean


class lvl0param:
    key_value_min = 0
    key_value_max = 1
    n = 630
    k = 1
    q = 2**16
    alpha = 0.000_092_511_997_467_675_6 * q
    α = alpha
    σ = alpha**2
    expectation_key_coefficient, variance_key_coefficient = _secret_stats(
        key_value_min, key_value_max
    )
    plain_modulus = 8


class lvlhalfparam:
    key_value_min = 0
    key_value_max = 1
    n = 760
    k = 1
    q = 2**32
    alpha = (2**-17) * q
    α = alpha
    σ = alpha**2
    expectation_key_coefficient, variance_key_coefficient = _secret_stats(
        key_value_min, key_value_max
    )
    plain_modulus = 32


class lvl1param:
    key_value_min = -1
    key_value_max = 1
    nbit = 10
    n = 2**nbit
    k = 1
    l = 3
    lₐ = l
    ℬbit = 6
    ℬₐbit = ℬbit
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    l̅ = 1
    l̅ₐ = l̅
    B̅gbit = 32
    B̅gₐbit = B̅gbit
    q = 2**32
    alpha = (2**-25) * q
    α = alpha
    σ = alpha**2
    expectation_key_coefficient, variance_key_coefficient = _secret_stats(
        key_value_min, key_value_max
    )
    plain_modulus = 8


class lvl2param:
    key_value_min = -1
    key_value_max = 1
    nbit = 11
    n = 2**nbit
    k = 1
    l = 4
    lₐ = l
    ℬbit = 9
    ℬₐbit = ℬbit
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    l̅ = 1
    l̅ₐ = l̅
    B̅gbit = 64
    B̅gₐbit = B̅gbit
    q = 2**64
    alpha = (2**-51) * q
    α = alpha
    σ = alpha**2
    expectation_key_coefficient, variance_key_coefficient = _secret_stats(
        key_value_min, key_value_max
    )
    plain_modulus = 8


class lvl10param:
    t = 7
    basebit = 2
    domainP = lvl1param
    targetP = lvl0param


class lvl1hparam:
    t = 10
    basebit = 3
    domainP = lvl1param
    targetP = lvlhalfparam


class lvl21param:
    t = 8
    basebit = 3
    domainP = lvl2param
    targetP = lvl1param


class lvl22param:
    t = 38
    basebit = 1
    domainP = lvl2param
    targetP = lvl2param


class lvl2hparam:
    t = 7
    basebit = 2
    domainP = lvl2param
    targetP = lvlhalfparam


class lvl02param:
    domainP = lvl0param
    targetP = lvl2param


class SS2CLPXlvl2param(lvl2param):
    l = 3
    lₐ = l
    ℬbit = 13
    ℬₐbit = ℬbit
    ℬ = 2**ℬbit
    ℬₐ = 2**ℬₐbit
    plain_modulus = 2


class SS2CLPXlvl02param:
    domainP = lvl0param
    targetP = SS2CLPXlvl2param


class SS2CLPXlvlh2param:
    domainP = lvlhalfparam
    targetP = SS2CLPXlvl2param


class SS2CLPXlvl22param:
    t = lvl22param.t
    basebit = lvl22param.basebit
    domainP = SS2CLPXlvl2param
    targetP = SS2CLPXlvl2param


class lvlh1param:
    domainP = lvlhalfparam
    targetP = lvl1param


class lvlh2param:
    domainP = lvlhalfparam
    targetP = lvl2param


class CLPX2TFHElvl2param(lvl2param):
    """Reverse-switch-only lvl2 PBS output with balanced gadget noise.

    The ring, secret, fresh-noise, and row count are unchanged.  Searching the
    four-row gadget decomposition over integral base bits selects Bgbit=10;
    the default Bgbit=9 leaves too little HomDecomp margin for basebit=4.
    """

    l = 4
    lₐ = l
    ℬbit = 10
    ℬₐbit = ℬbit
    ℬ = 2**ℬbit
    ℬₐ = ℬ


class CLPX2TFHElvlh2param:
    domainP = lvlhalfparam
    targetP = CLPX2TFHElvl2param
