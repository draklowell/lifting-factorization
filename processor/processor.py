from io import TextIOBase

from sage.all import *

from lifting.matrix_estimation import from_vector_y
from lifting.utils import min_degree, polynomial_ss, vector_ss
from processor.factorizer import Factorizer
from processor.matrix_estimator import MatrixEstimator
from processor.utils import (
    build_polyphase_matrix,
    matrix_ss,
    reconstruct,
    serialize_lifting_scheme,
)


class Processor:
    def __init__(self, output: TextIOBase):
        self.output = output

    def print(self, *args):
        self.output.write(" ".join(map(str, args)) + "\n")

    def estimate(self, H):
        estimator = MatrixEstimator(H)
        _, x_initial = estimator.build_system()

        y, eps = estimator.select_target_monomial(return_eps=True)

        det_H0 = from_vector_y(y, estimator.spec)
        residual = det_H0 - det(H)

        self.print(f"Selected determinant: {det_H0.change_ring(RR)}")
        self.print(f"Estimated epsilon: {float(eps)}")
        self.print(
            f"Determinant approximation L2 distance: {sqrt(float(polynomial_ss(residual)))}"
        )

        _, x_basis = estimator.build_solution_space()
        self.print(f"Kernel DoF: {x_basis.ncols()}")

        x = estimator.solve()

        x_residual = x - x_initial
        self.print(
            f"Row approximation L2 distance (non-normalized): {sqrt(float(vector_ss(x_residual)))}"
        )

        x = estimator.normalize(1)

        x_residual = x - x_initial
        self.print(
            f"Row approximation L2 distance (normalized): {sqrt(float(vector_ss(x_residual)))}"
        )

        return estimator.reconstruct()

    def factorize(self, P):
        factorizer = Factorizer(P)

        factorizer.normalize()

        _, gcd = factorizer.euclidean()

        self.print("Euclidean pre-factorization computed")
        self.print(f"GCD degree: {min_degree(gcd)}")

        P0 = factorizer.reconstruct()
        self.print(f"P0 matrix reconstructed with det P0 = {det(P0)}")

        factorizer.recover()
        self.print("Recovery step completed")

        return factorizer.pack()

    def process(self, low_pass, high_pass, F=QQ):
        assert len(low_pass) == len(high_pass)

        R = LaurentPolynomialRing(F, names=("z",))

        P = build_polyphase_matrix(low_pass, high_pass, R)

        P_est = self.estimate(P)

        steps, delays = self.factorize(P_est)

        P_rec = reconstruct(steps, delays, R)

        assert (
            P_rec == P_est
        ), "Reconstructed matrix does not match the estimated matrix"

        self.print(f"Sanity check: OK")

        l2 = sqrt(float(matrix_ss(P, P_rec)))
        self.print(f"Reconstruction L2 distance in QQ mode: {l2}")

        R_fp64 = LaurentPolynomialRing(RR, names=("z",))
        steps_fp64 = [(q.change_ring(RR), step) for q, step in steps]
        P_rec_fp64 = reconstruct(steps_fp64, delays, R_fp64)

        l2_fp64 = sqrt(float(matrix_ss(P, P_rec_fp64)))

        self.print(f"Reconstruction L2 distance in FP64 mode: {l2_fp64}")

        return serialize_lifting_scheme(
            steps,
            delays,
            len(low_pass),
            metadata={
                "l2": {
                    "qq": l2,
                    "fp64": l2_fp64,
                }
            },
        )
