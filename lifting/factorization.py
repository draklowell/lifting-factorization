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

from sage.all import *

from lifting.division import ldiv, ldivmod
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
