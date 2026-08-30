from dataclasses import dataclass
from enum import Enum

from sage.all import matrix, vector

from lifting.matrix_estimation import solve
from lifting.utils import max_degree, min_degree


@dataclass(frozen=True)
class VariableEntry:
    row: int
    column: int
    multiplier_row: int
    multiplier_column: int
    sign: int


HIGH_PASS_ENTRIES = (
    VariableEntry(0, 1, 1, 0, -1),
    VariableEntry(1, 1, 0, 0, 1),
)


class ProjectionStrategy(Enum):
    FIXED_ROW = "fixed-row"
    PRESERVE_LOW_PASS = "preserve-low-pass"


def dominant_determinant_degree(polyphase_matrix) -> int:
    coefficients = polyphase_matrix.det().monomial_coefficients()
    if not coefficients:
        raise ValueError("Polyphase determinant is zero")
    return max(coefficients, key=lambda exponent: abs(coefficients[exponent]))


def project_high_pass_to_unit_determinant(polyphase_matrix):
    supports = [
        range(min_degree(polynomial), max_degree(polynomial) + 1)
        for polynomial in (
            polyphase_matrix[entry.row, entry.column] for entry in HIGH_PASS_ENTRIES
        )
    ]

    target_degree = dominant_determinant_degree(polyphase_matrix)
    output_minimum = target_degree
    output_maximum = target_degree
    for entry, support in zip(HIGH_PASS_ENTRIES, supports):
        multiplier = polyphase_matrix[entry.multiplier_row, entry.multiplier_column]
        output_minimum = min(
            output_minimum,
            support.start + min_degree(multiplier),
        )
        output_maximum = max(
            output_maximum,
            support.stop - 1 + max_degree(multiplier),
        )

    output_size = output_maximum - output_minimum + 1
    variable_size = sum(len(support) for support in supports)
    field = polyphase_matrix.base_ring().base_ring()
    system = matrix(field, output_size, variable_size)
    initial = vector(field, variable_size)

    column = 0
    for entry, support in zip(HIGH_PASS_ENTRIES, supports):
        variable_coefficients = polyphase_matrix[
            entry.row, entry.column
        ].monomial_coefficients()
        multiplier_coefficients = polyphase_matrix[
            entry.multiplier_row, entry.multiplier_column
        ].monomial_coefficients()
        for exponent in support:
            initial[column] = variable_coefficients.get(exponent, 0)
            for multiplier_exponent, coefficient in multiplier_coefficients.items():
                row = exponent + multiplier_exponent - output_minimum
                system[row, column] += entry.sign * coefficient
            column += 1

    target = vector(field, output_size)
    target[target_degree - output_minimum] = 1
    particular = system.solve_right(target)
    basis = system.right_kernel().basis_matrix().transpose()
    projected = solve(initial, particular, basis)

    result = matrix(polyphase_matrix)
    ring = polyphase_matrix.base_ring()
    z = ring.gen()
    offset = 0
    for entry, support in zip(HIGH_PASS_ENTRIES, supports):
        polynomial = ring(0)
        for index, exponent in enumerate(support):
            polynomial += projected[offset + index] * z**exponent
        result[entry.row, entry.column] = polynomial
        offset += len(support)

    assert result.det() == z**target_degree
    return result
