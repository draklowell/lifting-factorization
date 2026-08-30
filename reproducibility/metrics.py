import math
from collections.abc import Callable

import numpy as np

Polynomial = dict[int, object]
PolyphaseMatrix = list[list[Polynomial]]


def polynomial_add(
    left: Polynomial,
    right: Polynomial,
    cast: Callable,
) -> Polynomial:
    result = dict(left)
    zero = cast(0)
    for degree, coefficient in right.items():
        value = cast(result.get(degree, zero) + coefficient)
        if value == zero:
            result.pop(degree, None)
        else:
            result[degree] = value
    return result


def polynomial_multiply(
    left: Polynomial,
    right: Polynomial,
    cast: Callable,
) -> Polynomial:
    result: Polynomial = {}
    zero = cast(0)
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            product = cast(left_coefficient * right_coefficient)
            result[degree] = cast(result.get(degree, zero) + product)
    return {degree: value for degree, value in result.items() if value != zero}


def multiply_polyphase(
    left: PolyphaseMatrix,
    right: PolyphaseMatrix,
    cast: Callable,
) -> PolyphaseMatrix:
    result = [[{}, {}], [{}, {}]]
    for row in range(2):
        for column in range(2):
            value: Polynomial = {}
            for inner in range(2):
                value = polynomial_add(
                    value,
                    polynomial_multiply(
                        left[row][inner],
                        right[inner][column],
                        cast,
                    ),
                    cast,
                )
            result[row][column] = value
    return result


def analysis_polyphase(low_pass, high_pass, cast=float) -> PolyphaseMatrix:
    result = [[{}, {}], [{}, {}]]
    for column, coefficients in enumerate((low_pass, high_pass)):
        for index, coefficient in enumerate(coefficients):
            row = index % 2
            degree = index // 2 if row == 0 else (index + 1) // 2
            result[row][column][degree] = cast(coefficient)
    return result


def sage_polyphase_to_numeric(polyphase_matrix, cast=float) -> PolyphaseMatrix:
    return [
        [
            {
                int(degree): cast(coefficient)
                for degree, coefficient in polyphase_matrix[row, column]
                .monomial_coefficients()
                .items()
            }
            for column in range(2)
        ]
        for row in range(2)
    ]


def scheme_polyphase(scheme: dict, cast=float) -> PolyphaseMatrix:
    result = [
        [{int(scheme["delay"]["even"]): cast(1)}, {}],
        [{}, {int(scheme["delay"]["odd"]): cast(1)}],
    ]
    for step in scheme["steps"]:
        polynomial = {
            int(step["shift"]) + index: cast(coefficient_value(coefficient))
            for index, coefficient in enumerate(step["coefficients"])
        }
        step_type = step["type"]
        if step_type == "predict":
            step_matrix = [[{0: cast(1)}, polynomial], [{}, {0: cast(1)}]]
        elif step_type == "update":
            step_matrix = [[{0: cast(1)}, {}], [polynomial, {0: cast(1)}]]
        elif step_type == "scale-even":
            step_matrix = [[polynomial, {}], [{}, {0: cast(1)}]]
        elif step_type == "scale-odd":
            step_matrix = [[{0: cast(1)}, {}], [{}, polynomial]]
        elif step_type == "swap":
            step_matrix = [[{}, {0: cast(1)}], [{0: cast(1)}, {}]]
        else:
            raise ValueError(f"Unsupported lifting step type: {step_type}")
        result = multiply_polyphase(result, step_matrix, cast)
    return result


def coefficient_value(value) -> float:
    if isinstance(value, dict):
        return value["numerator"] / value["denominator"]
    return value


def matrix_l2_error(reference: PolyphaseMatrix, actual: PolyphaseMatrix) -> float:
    squared_error = 0.0
    for row in range(2):
        for column in range(2):
            degrees = reference[row][column].keys() | actual[row][column].keys()
            for degree in degrees:
                difference = float(
                    reference[row][column].get(degree, 0)
                    - actual[row][column].get(degree, 0)
                )
                squared_error += difference * difference
    return math.sqrt(squared_error)


def evaluate_polynomial(polynomial: Polynomial, value: complex) -> complex:
    return sum(
        complex(coefficient) * value**degree
        for degree, coefficient in polynomial.items()
    )


def paraunitary_residual(polyphase: PolyphaseMatrix, samples: int = 512) -> float:
    maximum = 0.0
    identity = np.eye(2, dtype=np.complex128)
    for angle in np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False):
        value = np.exp(1j * angle)
        evaluated = np.array(
            [
                [
                    evaluate_polynomial(polyphase[row][column], value)
                    for column in range(2)
                ]
                for row in range(2)
            ],
            dtype=np.complex128,
        )
        maximum = max(
            maximum,
            float(np.linalg.norm(evaluated.conj().T @ evaluated - identity)),
        )
    return maximum
