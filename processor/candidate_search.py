from collections.abc import Callable
from dataclasses import dataclass

from processor.utils import serialize_lifting_scheme


@dataclass(frozen=True)
class ScoredFactorization:
    score: tuple[float, ...]
    steps: list
    delays: tuple[int, int]


def coefficient_value(value) -> float:
    if isinstance(value, dict):
        return value["numerator"] / value["denominator"]
    return float(value)


def scheme_complexity(scheme: dict) -> tuple[float, int, int, int, int]:
    coefficients = []
    for step in scheme["steps"]:
        for value in step["coefficients"]:
            coefficient = coefficient_value(value)
            if coefficient != 0.0:
                coefficients.append(abs(coefficient))
    even_dependency = (0, 0)
    odd_dependency = (0, 0)
    maximum_dependency_span = 1
    maximum_stencil_radius = 0

    for step in scheme["steps"]:
        coefficient_count = len(step["coefficients"])
        if coefficient_count:
            support = (
                step["shift"],
                step["shift"] + coefficient_count - 1,
            )
            maximum_stencil_radius = max(
                maximum_stencil_radius,
                abs(support[0]),
                abs(support[1]),
            )
            if step["type"] == "predict":
                source = (
                    even_dependency[0] + support[0],
                    even_dependency[1] + support[1],
                )
                odd_dependency = (
                    min(odd_dependency[0], source[0]),
                    max(odd_dependency[1], source[1]),
                )
            elif step["type"] == "update":
                source = (
                    odd_dependency[0] + support[0],
                    odd_dependency[1] + support[1],
                )
                even_dependency = (
                    min(even_dependency[0], source[0]),
                    max(even_dependency[1], source[1]),
                )
        elif step["type"] == "swap":
            even_dependency, odd_dependency = odd_dependency, even_dependency

        maximum_dependency_span = max(
            maximum_dependency_span,
            even_dependency[1] - even_dependency[0] + 1,
            odd_dependency[1] - odd_dependency[0] + 1,
        )

    return (
        max(coefficients, default=0.0),
        sum(len(step["coefficients"]) for step in scheme["steps"]),
        len(scheme["steps"]),
        maximum_dependency_span,
        maximum_stencil_radius,
    )


def candidate_score(
    scheme: dict,
    score: Callable[[dict], tuple[float, ...]],
) -> tuple[float, ...]:
    return (*score(scheme), *scheme_complexity(scheme))


def select_factorization(
    candidates,
    tap_size: int,
    score: Callable[[dict], tuple[float, ...]],
) -> ScoredFactorization:
    scored = []
    for steps, delays in candidates:
        try:
            scheme = serialize_lifting_scheme(
                steps,
                delays,
                tap_size,
                coefficients_as_float=True,
            )
            scored.append(
                ScoredFactorization(
                    score=candidate_score(scheme, score),
                    steps=steps,
                    delays=delays,
                )
            )
        except (FloatingPointError, OverflowError, ValueError):
            continue

    if not scored:
        raise ValueError("No finite FP32 factorization candidate found")
    return min(scored, key=lambda candidate: candidate.score)
