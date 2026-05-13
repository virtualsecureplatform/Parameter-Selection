#!/usr/bin/env python3
"""BFV plaintext polynomials used by the TFHEpp bootstrap estimator."""

from __future__ import annotations

import math


def _trim(poly: list[int]) -> None:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()


def _poly_mul_raw(a: list[int], b: list[int], mod: int) -> list[int]:
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj != 0:
                res[i + j] = (res[i + j] + ai * bj) % mod
    _trim(res)
    return res


def _reduce_mod_monic(poly: list[int], monic_modulus: list[int], mod: int) -> None:
    degree = len(monic_modulus) - 1
    if not monic_modulus or monic_modulus[-1] != 1:
        raise ValueError("modulus polynomial must be monic")
    while len(poly) > degree:
        shift = len(poly) - degree - 1
        lead = poly.pop()
        if lead == 0:
            continue
        for i in range(degree):
            poly[shift + i] = (poly[shift + i] - lead * monic_modulus[i]) % mod
    _trim(poly)


def _poly_mul_mod_monic(
    a: list[int], b: list[int], monic_modulus: list[int], mod: int
) -> list[int]:
    res = _poly_mul_raw(a, b, mod)
    _reduce_mod_monic(res, monic_modulus, mod)
    return res


def _poly_pow_mod_monic(
    base: list[int], exp: int, monic_modulus: list[int], mod: int
) -> list[int]:
    res = [1 % mod]
    base = list(base)
    _reduce_mod_monic(base, monic_modulus, mod)
    while exp != 0:
        if exp & 1:
            res = _poly_mul_mod_monic(res, base, monic_modulus, mod)
        exp >>= 1
        if exp != 0:
            base = _poly_mul_mod_monic(base, base, monic_modulus, mod)
    return res


def _centered_range_null_poly(bound: int, mod: int) -> list[int]:
    poly = [1]
    for i in range(-bound, bound + 1):
        poly = _poly_mul_raw(poly, [(-i) % mod, 1], mod)
    return poly


def _try_reduce_mod_p_times_monic(
    poly: list[int], monic_poly: list[int], p: int
) -> bool:
    mod = p * p
    degree = len(monic_poly) - 1
    if not monic_poly or monic_poly[-1] != 1:
        raise ValueError("modulus polynomial must be monic")

    while len(poly) > degree:
        shift = len(poly) - degree - 1
        lead = poly.pop()
        if lead == 0:
            continue
        if lead % p != 0:
            return False
        quotient = lead // p
        for i in range(degree):
            poly[shift + i] = (
                poly[shift + i] - quotient * p * monic_poly[i]
            ) % mod
    _trim(poly)
    return True


def lowest_digit_removal_polynomial_over_range(p: int, bound: int) -> list[int]:
    """
    Return TFHEpp's bounded low-digit removal polynomial over Z/(p^2).

    The polynomial satisfies f(p*m + e) = p*m mod p^2 for all e in
    [-bound, bound].  This mirrors
    `GetLowestDigitRemovalPolynomialOverRange(p, B)` in
    `TFHEpp/include/bfv-digitext.hpp`.
    """
    if p < 2:
        raise ValueError("p must be at least 2")
    if bound < 0:
        raise ValueError("bound must be non-negative")
    if 2 * bound + 1 > p:
        raise ValueError("digit-error range must fit in one p-digit")

    mod = p * p
    base_poly = _centered_range_null_poly(bound, mod)
    square_modulus = _poly_mul_raw(base_poly, base_poly, mod)

    result = [0] * (len(square_modulus) - 1)
    for e in range(-bound, bound + 1):
        if e == 0:
            continue
        y_minus_e = [(-e) % mod, 1]
        indicator = _poly_pow_mod_monic(y_minus_e, p, square_modulus, mod)
        indicator = _poly_pow_mod_monic(indicator, p - 1, square_modulus, mod)

        if len(result) < len(indicator):
            result.extend([0] * (len(indicator) - len(result)))
        result[0] = (result[0] + e) % mod
        for i, coeff in enumerate(indicator):
            result[i] = (result[i] - e * coeff) % mod

    removal = [0] * max(2, len(result))
    removal[1] = 1
    for i, coeff in enumerate(result):
        removal[i] = (removal[i] - coeff) % mod

    reduced = list(removal)
    if _try_reduce_mod_p_times_monic(reduced, base_poly, p):
        removal = reduced
    _trim(removal)
    return removal


def prime_from_prime_square(t: int) -> int:
    p = math.isqrt(t)
    if p * p != t:
        raise ValueError("plaintext modulus is not a prime square; pass --digit-prime")
    return p
