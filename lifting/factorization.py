"""
Pipeline:

P, delays = normalize_det(P)

qs, a = euclidean(P[0, 0], P[1, 0])
P0 = reconstruct(qs, a)
s, b = recover(P, P0)

steps = pack(qs, a, s, b, P)

# Output: steps, delays
"""

from enum import Enum

from sage.all import det, matrix

from lifting.division import ldiv, ldivmod, ldivmod_candidates
from lifting.utils import degree


class LiftingStep(Enum):
    PREDICT = 1
    UPDATE = 2
    SCALE_EVEN = 3
    SCALE_ODD = 4
    SWAP = 5


def normalize_det(P):
    R = P.base_ring()
    z = R.gen()

    delay = det(P).degree()
    delay_up = delay // 2
    delay_down = delay - delay_up

    return (
        matrix(
            R,
            [
                [z**-delay_up, 0],
                [0, z**-delay_down],
            ],
        )
        * P
    ), (delay_up, delay_down)


def euclidean(a, b):
    """
    Perform euclidean algorithm over polynomials a and b
    using symmetric division strategy.

    Returns list of quotients and greatest common divisor.
    """

    # Stores all q factors
    result = []
    if degree(a) < degree(b):
        a, b = b, a
        result.append(0)  # Swap term.

    # a(i+1) = b(i)
    # b(i+1) = a(i) % b(i)
    # q(i+1) = a(i) // b(i)
    while True:
        q, r = ldivmod(a, b)

        if degree(r) >= degree(b):
            raise ValueError("Remainder degree must be less than divisor degree")

        a = b
        b = r
        result.append(q)

        if degree(r) == -1:
            break

    return result, a


def quotient_path_score(qs):
    maximum_log2 = float("-inf")
    accumulated_log2 = 0.0
    total_taps = 0
    for q in qs:
        coefficients = q.monomial_coefficients()
        total_taps += len(coefficients)
        for coefficient in coefficients.values():
            if coefficient == 0:
                continue
            magnitude = abs(coefficient)
            magnitude_log2 = (
                magnitude.numerator().nbits() - magnitude.denominator().nbits()
            )
            maximum_log2 = max(maximum_log2, magnitude_log2)
            accumulated_log2 += max(0.0, magnitude_log2)
    return maximum_log2, accumulated_log2, total_taps, len(qs)


def euclidean_candidates(a, b, beam_width=16):
    if beam_width < 1:
        raise ValueError("Beam width must be positive")

    initial_qs = []
    if degree(a) < degree(b):
        a, b = b, a
        initial_qs.append(a.parent()(0))

    active = [(a, b, initial_qs)]
    completed = []
    while active:
        expanded = []
        for current_a, current_b, qs in active:
            for q, remainder in ldivmod_candidates(current_a, current_b):
                candidate_qs = qs + [q]
                if degree(remainder) == -1:
                    completed.append((candidate_qs, current_b))
                else:
                    expanded.append((current_b, remainder, candidate_qs))

        completed.sort(key=lambda candidate: quotient_path_score(candidate[0]))
        del completed[beam_width:]
        expanded.sort(key=lambda candidate: quotient_path_score(candidate[2]))
        active = expanded[:beam_width]

    if not completed:
        raise ValueError("Laurent Euclidean search produced no complete path")
    return completed


def factorization_candidates(P, beam_width=16):
    candidates = []
    paths = [euclidean(P[0, 0], P[1, 0])]
    for qs, a in euclidean_candidates(P[0, 0], P[1, 0], beam_width):
        if any(qs == existing_qs for existing_qs, _ in paths):
            continue
        paths.append((qs, a))

    for qs, a in paths:
        P0 = reconstruct(qs, a)
        s, b = recover(P, P0)
        steps = pack(qs, a, s, b, P)
        candidates.append(steps)
    return candidates


def reconstruct(qs, a):
    """
    Reconstruct the matrix P0 from the q factors and a, such that P0 has the form:
    P0 = [q_1 1  * ... * [a  0
           1  0]          0 1/a]
    """
    R = a.parent()

    P0 = matrix(R, 2, 2, [[a, 0], [0, 1 / a]])
    for q in reversed(qs):
        Q = matrix(R, 2, 2, [[q, 1], [1, 0]])
        P0 = Q * P0

    return P0


def recover(P, P0):
    """
    Compute the scaling factor s and b such that P = P0 * M, where M is a following matrix:
    [1 s
     0 b]
    """
    d = det(P0)
    b = 1 / d

    he, ho, ge, go = P[0, 0], P[1, 0], P[0, 1], P[1, 1]
    ge0, go0 = P0[0][1], P0[1][1]
    se = ldiv(ge - b * ge0, he)
    so = ldiv(go - b * go0, ho)
    assert se == so

    return se, b


def pack(qs, a, s, b, P):
    """
    Pack the factorization steps into a list of (q, step_type) tuples, where step_type is one of the LiftingStep enum values.
    """
    steps = []
    for i, q in enumerate(qs):
        step_type = LiftingStep.PREDICT if i % 2 == 0 else LiftingStep.UPDATE
        if q == 0:
            continue

        steps.append((q, step_type))

    if len(qs) % 2 == 1:
        R = P.base_ring()
        steps.append((R(0), LiftingStep.SWAP))

    steps.append((a * a * s / b, LiftingStep.PREDICT))
    steps.append((a, LiftingStep.SCALE_EVEN))
    steps.append((b / a, LiftingStep.SCALE_ODD))

    return steps


def canonicalize_scale_delays(steps, delays):
    """Move monomial scale delays into the leading delay matrix."""
    adjusted = list(steps)
    even_exponent = 0
    odd_exponent = 0
    ring = steps[0][0].parent()
    z = ring.gen()

    for index, (polynomial, step_type) in enumerate(adjusted):
        if step_type not in {LiftingStep.SCALE_EVEN, LiftingStep.SCALE_ODD}:
            continue
        coefficients = polynomial.monomial_coefficients()
        if len(coefficients) != 1:
            raise ValueError("Scale step must contain one monomial")
        exponent, coefficient = next(iter(coefficients.items()))
        adjusted[index] = (ring(coefficient), step_type)
        if step_type == LiftingStep.SCALE_EVEN:
            even_exponent += exponent
        else:
            odd_exponent += exponent

    if even_exponent + odd_exponent != 0:
        raise ValueError("Scale monomial delays do not preserve determinant degree")

    for index in range(len(adjusted) - 1, -1, -1):
        polynomial, step_type = adjusted[index]
        if step_type == LiftingStep.PREDICT:
            adjusted[index] = (
                polynomial * z ** (odd_exponent - even_exponent),
                step_type,
            )
        elif step_type == LiftingStep.UPDATE:
            adjusted[index] = (
                polynomial * z ** (even_exponent - odd_exponent),
                step_type,
            )
        elif step_type == LiftingStep.SWAP:
            even_exponent, odd_exponent = odd_exponent, even_exponent

    return adjusted, (
        delays[0] + even_exponent,
        delays[1] + odd_exponent,
    )


def factorize(P, normalize: bool = False):
    if normalize:
        P, delays = normalize_det(P)
    else:
        delays = None

    qs, a = euclidean(P[0, 0], P[1, 0])
    P0 = reconstruct(qs, a)
    s, b = recover(P, P0)

    steps = pack(qs, a, s, b, P)

    return steps, delays
