from sage.all import *

from lifting.factorization import LiftingStep
from lifting.utils import max_degree, min_degree


def split(filter_coeffs: list[float], R):
    F = R.base_ring()
    z = R.gen()

    even = 0
    odd = 0
    for i, coeff in enumerate(filter_coeffs):
        if i % 2 == 0:
            even += F(coeff) * z ** Integer(i // 2)
        else:
            odd += F(coeff) * z ** Integer((i + 1) // 2)

    return even, odd


def build_polyphase_matrix(low_pass, high_pass, R):
    he, ho = split(low_pass, R)
    ge, go = split(high_pass, R)
    return matrix(R, [[he, ge], [ho, go]])


def build_delay_matrix(delay_up, delay_down, R):
    z = R.gen()
    return matrix(
        R,
        [
            [z**delay_up, 0],
            [0, z**delay_down],
        ],
    )


def build_step_matrix(q, step, R):
    if step == LiftingStep.PREDICT:
        return matrix(R, [[1, q], [0, 1]])
    if step == LiftingStep.UPDATE:
        return matrix(R, [[1, 0], [q, 1]])
    if step == LiftingStep.SCALE_EVEN:
        return matrix(R, [[q, 0], [0, 1]])
    if step == LiftingStep.SCALE_ODD:
        return matrix(R, [[1, 0], [0, q]])
    if step == LiftingStep.SWAP:
        return matrix(R, [[0, 1], [1, 0]])
    raise ValueError(f"Unknown lifting step: {step}")


def reconstruct(steps, delay, R):
    P = matrix(R, 2, 2, [[1, 0], [0, 1]])
    P *= build_delay_matrix(*delay, R)

    for q, step in steps:
        P *= build_step_matrix(q, step, R)

    return P


def serialize_polynomial(q):
    shift = min_degree(q)
    coeffs = []
    for degree_idx in range(min_degree(q), max_degree(q) + 1):
        coeff = q.coefficient(degree_idx)
        coeffs.append(
            {
                "numerator": int(coeff.numerator()),
                "denominator": int(coeff.denominator()),
            }
        )

    if coeffs[0]["numerator"] == 0:
        return {"shift": 0, "coefficients": []}

    return {"shift": shift, "coefficients": coeffs}


STEP_NAMES = {
    LiftingStep.PREDICT: "predict",
    LiftingStep.UPDATE: "update",
    LiftingStep.SCALE_EVEN: "scale-even",
    LiftingStep.SCALE_ODD: "scale-odd",
    LiftingStep.SWAP: "swap",
}


def serialize_lifting_scheme(steps, delay, tap_size: int, metadata: dict | None = None):
    result = {
        "tap_size": tap_size,
        "delay": {
            "even": int(delay[0]),
            "odd": int(delay[1]),
        },
        "steps": [],
    }

    if metadata is not None:
        result["meta"] = metadata

    for q, step in steps:
        q_serialized = serialize_polynomial(q)
        if len(q_serialized["coefficients"]) == 0 and step in {
            LiftingStep.PREDICT,
            LiftingStep.UPDATE,
        }:
            continue

        result["steps"].append(
            {
                "type": STEP_NAMES[step],
                **q_serialized,
            }
        )

    return result


def matrix_ss(a, b, F=QQ):
    accum = 0
    for poly_a, poly_b in zip(a.list(), b.list()):
        coeffs_a = poly_a.change_ring(F).dict()
        coeffs_b = poly_b.change_ring(F).dict()
        for m in set(coeffs_a.keys()) | set(coeffs_b.keys()):
            c_a = coeffs_a.get(m, 0)
            c_b = coeffs_b.get(m, 0)
            r = c_a - c_b
            accum += r * r

    return accum
