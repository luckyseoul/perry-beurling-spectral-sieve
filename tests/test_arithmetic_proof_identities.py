"""Finite exact algebra checks supporting the continuous proof documents.

These checks corroborate differential-operator signs, principal parts, and
recurrence constants. They do not test an all-T estimate or prove RH.
Fractions avoid numerical cancellation in the alternating Legendre sums.
"""
from fractions import Fraction
from math import comb, factorial

import pytest


def _negative_derivative(series, order):
    """Apply (-d/dx)^order to a finite Laurent polynomial."""
    result = dict(series)
    for _ in range(order):
        result = {
            power - 1: -power * coefficient
            for power, coefficient in result.items()
            if power and coefficient
        }
    return result


def _moment_operator_principal_part(coefficients, pole, residue, regular=None):
    """Expand the B-ABS operator at s=pole, with x=s-pole.

    The finite Taylor expansion of 1/s includes every term that can contribute
    to the principal part after at most k derivatives. Regular additions to F
    check that the result really depends on the pole alone.
    """
    degree = len(coefficients) - 1
    transform = {-1: Fraction(residue), **(regular or {})}
    result = {}
    for j, coefficient in enumerate(coefficients):
        differentiated = _negative_derivative(transform, j)
        quotient = {}
        for power, value in differentiated.items():
            for ell in range(degree + 1):
                exponent = power + ell
                term = value * (-1) ** ell / pole ** (ell + 1)
                quotient[exponent] = quotient.get(exponent, Fraction(0)) + term
        for exponent, value in _negative_derivative(quotient, degree - j).items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient * value
    return {power: value for power, value in result.items() if power < 0 and value}


def _endpoint_principal_part(coefficients, pole, residue):
    """Independent endpoint integration-by-parts expansion for Q(y)=r exp(a y).

    The exponentially growing endpoint polynomial has Laplace poles at a;
    the lower-endpoint polynomial has poles only at zero.
    """
    degree = len(coefficients) - 1
    result = {}
    for ell in range(degree + 1):
        derivative_at_one = sum(
            coefficients[j] * factorial(j) // factorial(j - ell)
            for j in range(ell, degree + 1)
        )
        coefficient = (
            residue * (-1) ** ell * factorial(degree - ell)
            * derivative_at_one / pole ** (ell + 1)
        )
        if coefficient:
            result[-(degree - ell + 1)] = coefficient
    return result


@pytest.mark.parametrize("degree", [0, 1, 2, 4, 7])
@pytest.mark.parametrize("pole", [Fraction(1, 3), Fraction(-2, 5)])
def test_legendre_moment_operator_preserves_pole_and_full_principal_part(degree, pole):
    # Coefficients of L_k(2u-1); orthonormal normalization multiplies every
    # principal coefficient by the same nonzero sqrt(2k+1).
    coefficients = [
        (-1) ** (degree + j) * comb(degree, j) * comb(degree + j, j)
        for j in range(degree + 1)
    ]
    residue = Fraction(-3, 2)
    actual = _moment_operator_principal_part(
        coefficients, pole, residue,
        regular={0: Fraction(11, 3), 1: Fraction(-5, 2), 3: Fraction(7, 4)},
    )
    assert actual == _endpoint_principal_part(coefficients, pole, residue)
    assert sum(coefficients) == 1
    assert min(actual) == -(degree + 1)
    assert actual[-(degree + 1)] == factorial(degree) * residue / pole


@pytest.mark.parametrize("degree", [1, 2, 4])
def test_endpoint_zero_really_removes_highest_pole(degree):
    # u^(k-1)(u-1) has p(1)=0, p'(1)=1. The endpoint nonvanishing
    # hypothesis in B-ABS is necessary for its stated highest pole order.
    coefficients = [0] * (degree + 1)
    coefficients[-2:] = [-1, 1]
    pole, residue = Fraction(2, 7), Fraction(3, 5)
    actual = _moment_operator_principal_part(coefficients, pole, residue)
    assert actual == _endpoint_principal_part(coefficients, pole, residue)
    assert -(degree + 1) not in actual
    assert actual[-degree] == -factorial(degree - 1) * residue / pole**2


@pytest.mark.parametrize("degree", [0, 1, 4, 8])
def test_energy_recurrence_constants_have_correct_boundary(degree):
    c = Fraction(1, degree + 2)
    delta = c**2 * (degree + 1)**2
    alpha = 1 - c**2 / 4
    multiplier = 4 * alpha / c**2
    assert 0 < delta <= 1
    assert 0 < alpha < 1 < multiplier
    assert 1 - alpha == (c / 2)**2
    assert alpha / multiplier == (c / 2)**2
    # If V(alpha T)/V(T) were only 1/(4K), the upper bound in (2)
    # would be c/4+c/2 < c, contradicting the positive-floor premise.
    assert alpha / (4 * multiplier) == (c / 4)**2
    assert c / 4 + c / 2 < c


@pytest.mark.parametrize("power", [0, 1, 3, 9])
@pytest.mark.parametrize("degree", [0, 2, 4])
def test_power_law_energy_satisfies_fixed_ratio_and_dyadic_bounds(power, degree):
    # Q(y)=y^power has V(T)=T^(2power+1)/(2power+1) and the exact
    # degree-zero ratio (2power+1)/(power+1)^2. Higher-degree projections
    # capture at least this much energy, making delta a valid lower floor.
    exponent = 2 * power + 1
    delta = Fraction(exponent, 2 * (power + 1)**2)
    c_squared = delta / (degree + 1)**2
    alpha = 1 - c_squared / 4
    multiplier = 4 * alpha / c_squared
    assert multiplier * alpha**exponent >= 1

    contraction = Fraction(1)
    steps = 0
    while contraction > Fraction(1, 2):
        contraction *= alpha
        steps += 1
    assert contraction <= Fraction(1, 2) < contraction / alpha
    assert multiplier**steps >= 2**exponent  # V(2T)/V(T)


@pytest.mark.parametrize("low,captured,outside", [
    (Fraction(2), Fraction(3), Fraction(5)),
    (Fraction(0), Fraction(3), Fraction(5)),
    (Fraction(2), Fraction(3), Fraction(0)),
    (Fraction(1000), Fraction(1, 100), Fraction(7)),
])
def test_removing_lower_degrees_cannot_increase_energy_ratio(low, captured, outside):
    raw = (low + captured) / (low + captured + outside)
    detrended = captured / (captured + outside)
    assert detrended == (raw - low / (low + captured + outside)) / (
        1 - low / (low + captured + outside)
    )
    assert raw - detrended == low * outside / (
        (low + captured + outside) * (captured + outside)
    )
    assert 0 <= detrended <= raw <= 1
