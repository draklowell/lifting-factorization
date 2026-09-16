"""
Pipeline:

spec = build_specification(H)

T, x_initial = build_system(H, spec)

y = select_target_monomial(H, spec, eps=eps)

x_particular, x_basis = build_solution_space(T, y)

x = solve(x_initial, x_particular, x_basis)

H_new = reconstruct(H, x, spec)

# Output: H_new
"""

from dataclasses import dataclass
from typing import Any

from sage.all import det, matrix, vector

from lifting.utils import (
    from_vector,
    laurent_toeplitz,
    max_degree,
    min_degree,
    to_vector,
)


@dataclass
class EstimationSpecification:
    R: Any
    lowest_x: int
    highest_x: int
    lowest_y: int
    highest_y: int

    @property
    def size_x(self):
        return self.highest_x - self.lowest_x + 1

    @property
    def size_y(self):
        return self.highest_y - self.lowest_y + 1

    @property
    def size_2x(self):
        return 2 * self.size_x

    @property
    def F(self):
        return self.R.base_ring()


def build_specification(H):
    # lowest degree in 1st vector representation
    alpha = min(
        min_degree(H[1, 0]),
        min_degree(H[1, 1]),
    )
    # highest degree in 1st vector representation
    beta = max(
        max_degree(H[1, 0]),
        max_degree(H[1, 1]),
    )

    # lower bound extension to alpha for the output vector
    l = min(
        min_degree(H[0, 0]),
        min_degree(H[0, 1]),
        0,
    )
    # upper bound extension to beta for the output vector
    k = max(
        max_degree(H[0, 0]),
        max_degree(H[0, 1]),
    )

    return EstimationSpecification(
        R=H.base_ring(),
        lowest_x=alpha,
        highest_x=beta,
        lowest_y=alpha + l,
        highest_y=beta + k,
    )


def to_vector_x(p, spec):
    return to_vector(p, spec.lowest_x, spec.size_x)


def to_vector_2x(p1, p2, spec):
    return vector(
        list(to_vector(p1, spec.lowest_x, spec.size_x))
        + list(to_vector(p2, spec.lowest_x, spec.size_x))
    )


def to_vector_y(p, spec):
    return to_vector(p, spec.lowest_y, spec.size_y)


def from_vector_x(v, spec):
    return from_vector(v, spec.lowest_x, spec.R)


def from_vector_2x(v, spec):
    return (
        from_vector(v[: spec.size_x], spec.lowest_x, spec.R),
        from_vector(v[spec.size_x :], spec.lowest_x, spec.R),
    )


def from_vector_y(v, spec):
    return from_vector(v, spec.lowest_y, spec.R)


def build_system(H, spec=None):
    """
    Build system of linear equations to determine space of matrices with second
    row fixed, second row being variable (which is returned as vector) and y being
    determinant.
    """
    spec = spec or build_specification(H)

    T1 = matrix(spec.F, spec.size_y, spec.size_2x)
    T2 = matrix(spec.F, spec.size_y, spec.size_2x)

    for i in range(spec.size_x):
        T1[i + spec.lowest_x - spec.lowest_y, i] = 1
        T2[i + spec.lowest_x - spec.lowest_y, spec.size_x + i] = 1

    Ta = laurent_toeplitz(H[0, 0], spec.size_y, spec.size_y)
    Tb = laurent_toeplitz(H[0, 1], spec.size_y, spec.size_y)

    T = Ta * T2 - Tb * T1

    v = to_vector_2x(H[1, 0], H[1, 1], spec)

    assert T * v == to_vector_y(det(H), spec)

    return T, v


def select_target_monomial(
    H,
    spec: EstimationSpecification,
    eps=None,
    return_eps: bool = False,
    normalize_to=None,
):
    """
    Selet monomial in determinant for which the system from `build_system`
    will be solved. Returns vector representation of the monomial.
    """
    det_H = det(H)

    coeff_max = None
    degree_max = None
    for degree_cur, coeff_cur in det_H.monomial_coefficients().items():
        if eps is None:
            if degree_max is None or abs(coeff_cur) > abs(coeff_max):
                degree_max = degree_cur
                coeff_max = coeff_cur
            continue

        if abs(coeff_cur) > eps:
            if degree_max is not None:
                raise ValueError("Multiple monomials with non-zero coefficients found.")

            degree_max = degree_cur
            coeff_max = coeff_cur

    R = det_H.parent()

    if eps is None:
        selected = R(0) if degree_max is None else coeff_max * R.gen() ** degree_max
        residual = det_H - selected
        coeffs = list(map(abs, residual.monomial_coefficients().values()))

        eps = max(coeffs) if coeffs else 0

    if degree_max is None:
        value = R(0)
    else:
        coefficient = coeff_max if normalize_to is None else R.base_ring()(normalize_to)
        value = coefficient * R.gen() ** degree_max

    if return_eps:
        return to_vector_y(value, spec), eps

    return to_vector_y(value, spec)


def build_solution_space(A, y):
    """
    Solve system of linear equations from `build_system` and return
    particular solution and basis of the homogenous linear system,
    which basically describe Affine space of solutions.
    """
    return A.solve_right(y), A.right_kernel().basis_matrix().transpose()


def solve(x_initial, x_particular, x_basis):
    """
    Do an orthogonal projection of our original matrix into a space of matrices with
    the same first row and determinant equal to monomial using method
    of least squares (specifically normal equations).
    """
    # No DoF
    if x_basis.ncols() == 0:
        return x_particular

    # x_i = x_p + X_b * c + eps
    # x_i - x_p = X_b * c + eps
    y = x_initial - x_particular
    A = x_basis
    # Method of least squares
    coordinate = (A.transpose() * A).solve_right(A.transpose() * y)

    return x_particular + x_basis * coordinate


def normalize(x, y, normalize_to):
    """
    Normalize second row (x) so that the coefficient near determinant (which
    is monomial) will equal to `normalize_to`.
    """

    value = list(filter(lambda c: c != 0, y))
    if len(value) != 1:
        raise ValueError("Expected exactly one monomial in determinant")

    F = x.base_ring()
    return x * F(normalize_to) / value[0]


def reconstruct(H, x, spec=None):
    """
    Given solution x, reconstruct matrix, which is orthogonal projection
    of H. See `solve`.
    """
    spec = spec or build_specification(H)

    h1, h2 = from_vector_2x(x, spec)

    return matrix(
        spec.R,
        [
            [H[0, 0], H[0, 1]],
            [h1, h2],
        ],
    )


def estimate_matrix(H, normalize_to=None, eps=None, return_eps=False):
    spec = build_specification(H)

    T, x_initial = build_system(H, spec)

    y, eps = select_target_monomial(
        H,
        spec,
        eps=eps,
        return_eps=True,
        normalize_to=normalize_to,
    )

    x_particular, x_basis = build_solution_space(T, y)

    x = solve(x_initial, x_particular, x_basis)

    H_new = reconstruct(H, x, spec)

    if return_eps:
        return H_new, eps

    return H_new
