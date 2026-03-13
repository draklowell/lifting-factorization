from sage.all import *

from lifting.utils import (
    from_vector,
    laurent_toeplitz,
    max_degree,
    min_degree,
    to_vector,
)


class MatrixEstimator:
    def __init__(self, H, output = None):
        self.H = H
        self.R = H.base_ring()
        self.F = self.R.base_ring()
        self.output = output

        # lowest degree in 1st vector representation
        self.alpha = min(
            min_degree(H[1, 0]),
            min_degree(H[1, 1]),
        )
        # highest degree in 1st vector representation
        self.beta = max(
            max_degree(H[1, 0]),
            max_degree(H[1, 1]),
        )

        self.size = self.beta - self.alpha + 1

        # lower bound extension to alpha for the output vector
        self.l = min(
            min_degree(H[0, 0]),
            min_degree(H[0, 1]),
            0,
        )
        # upper bound extension to beta for the output vector
        self.k = max(
            max_degree(H[0, 0]),
            max_degree(H[0, 1]),
        )

        self.size_out = self.k - self.l + self.beta - self.alpha + 1

    def print(self, *args):
        if self.output is not None:
            self.output.write(" ".join(map(str, args)) + "\n")

    # get vector representation of the polynomial
    def isomorphism(self, p):
        return to_vector(p, self.alpha, self.size)

    def isomorphism_inverse(self, v):
        return from_vector(v, self.alpha, self.R)

    def isomorphism2(self, p1, p2):
        return vector(list(self.isomorphism(p1)) + list(self.isomorphism(p2)))

    def isomorphism2_inverse(self, v):
        return (
            self.isomorphism_inverse(v[: self.size]),
            self.isomorphism_inverse(v[self.size :]),
        )

    def isomorphism_out(self, p):
        return to_vector(p, self.alpha + self.l, self.size_out)

    def isomorphism_out_inverse(self, v):
        return from_vector(v, self.alpha + self.l, self.R)

    def get_system(self):
        T1 = matrix(self.F, self.size_out, 2 * self.size)
        T2 = matrix(self.F, self.size_out, 2 * self.size)

        for i in range(self.size):
            T1[i - self.l, i] = 1
            T2[i - self.l, self.size + i] = 1

        Ta = laurent_toeplitz(self.H[0, 0], self.size_out, self.size_out)
        Tb = laurent_toeplitz(self.H[0, 1], self.size_out, self.size_out)

        T = Ta * T2 - Tb * T1

        v = self.isomorphism2(self.H[1, 0], self.H[1, 1])

        return T, v

    def l2_norm2_poly(self, p):
        coeffs = p.monomial_coefficients()
        norm_squared = sum(abs(c) ** 2 for c in coeffs.values())
        return norm_squared

    def l2_norm2_vector(self, v):
        norm_squared = sum(abs(c) ** 2 for c in v)
        return norm_squared

    def get_target(self, eps=None):
        det_H = det(self.H)

        coeff_max = None
        degree_max = None
        for m, c in det_H.monomial_coefficients().items():
            if eps is None:
                if degree_max is None:
                    degree_max = m
                    coeff_max = c
                elif abs(c) > abs(coeff_max):
                    degree_max = m
                    coeff_max = c

                continue

            if abs(c) > eps:
                if degree_max is not None:
                    raise ValueError(
                        "Multiple monomials with non-zero coefficients found."
                    )

                degree_max = m
                coeff_max = c

        z = det_H.parent().gen()

        if degree_max is None:
            value = self.R(0)
        else:
            value = coeff_max * z**degree_max

        removed = det_H - value
        if eps is None:
            eps = None
            for m, c in removed.monomial_coefficients().items():
                if eps is None:
                    eps = abs(c)
                else:
                    eps = max(eps, abs(c))

            if eps is None:
                eps = 0

        self.print(f"Selected monomial: z^{degree_max}, Coefficient: {float(coeff_max)}")
        self.print(f"Estimated epsilon: {float(eps)}")
        self.print(
            f"Determinant approximation L2 distance: {sqrt(float(self.l2_norm2_poly(removed)))}"
        )

        return self.isomorphism_out(value)

    def regression(self, X, y):
        # Solve the least squares problem X * c = y
        return (X.transpose() * X).solve_right(X.transpose() * y)

    def solve_system(self, T, v, d, normalize_to=None):
        v_p = T.solve_right(d)

        # Columns are basis vectors of the kernel
        V_g = T.right_kernel().basis_matrix().transpose()

        # v = v_p + V_g * c + eps
        # v - v_p = V_g * c + eps
        c = self.regression(V_g, v - v_p)

        v_new = v_p + V_g * c
        self.print(f"Kernel DoF: {V_g.ncols()}")
        v_diff = v_new - v
        self.print(
            f"Row approximation L2 distance (non-normalized): {sqrt(float(self.l2_norm2_vector(v_diff)))}"
        )

        if normalize_to is not None:
            norm_coeff = None
            for c in d:
                if c != 0:
                    if norm_coeff is not None:
                        raise ValueError("Multiple non-zero coefficients found in d.")
                    norm_coeff = c

            if norm_coeff is None:
                raise ValueError(
                    "No non-zero coefficients found in d for normalization."
                )

            factor = self.F(normalize_to) / norm_coeff
            v_new *= factor

            self.print(f"Normalization factor: {float(factor)}")
            self.print(f"Row approximation L2 distance (normalized): {sqrt(float(self.l2_norm2_vector(v_new - v)))}")

        return v_new

    @staticmethod
    def solve(H, normalize_to, eps=None, output=None):
        estimator = MatrixEstimator(H, output)
        T, v_old = estimator.get_system()

        # Perform sanity check
        assert T * v_old == estimator.isomorphism_out(
            det(H)
        ), "T * V does not equal isomorphism_out(det(H))"

        d = estimator.get_target(eps=eps)

        v_new = estimator.solve_system(
            T, v_old, d, normalize_to=normalize_to
        )

        ho_new, go_new = estimator.isomorphism2_inverse(v_new)
        return matrix(
            H.base_ring(),
            [
                [H[0, 0], H[0, 1]],
                [ho_new, go_new],
            ],
        )
