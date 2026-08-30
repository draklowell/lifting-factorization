from collections.abc import Callable
from io import TextIOBase

from sage.all import QQ, RDF, LaurentPolynomialRing, det, matrix, sqrt

from lifting.matrix_estimation import from_vector_y
from lifting.projection import ProjectionStrategy, project_high_pass_to_unit_determinant
from lifting.utils import min_degree, polynomial_ss, vector_ss
from processor.factorizer import Factorizer
from processor.matrix_estimator import MatrixEstimator
from processor.utils import (
    build_polyphase_matrix,
    matrix_ss,
    reconstruct,
    serialize_lifting_scheme,
)

Factorization = tuple[list, tuple[int, int]]
CandidateSelector = Callable[[list[Factorization]], Factorization]


class Processor:
    def __init__(self, output: TextIOBase):
        self.output = output

    def print(self, *args):
        self.output.write(" ".join(map(str, args)) + "\n")

    def estimate(
        self,
        H,
        strategy: ProjectionStrategy = ProjectionStrategy.FIXED_ROW,
    ):
        if strategy is ProjectionStrategy.PRESERVE_LOW_PASS:
            estimate = project_high_pass_to_unit_determinant(H)
            self.print("Selected low-pass-preserving projection")
            return estimate
        if strategy is not ProjectionStrategy.FIXED_ROW:
            raise ValueError(f"Unsupported projection strategy: {strategy}")

        estimator = MatrixEstimator(H)
        _, x_initial = estimator.build_system()

        y, eps = estimator.select_target_monomial(
            return_eps=True,
            normalize_to=1,
        )

        det_H0 = from_vector_y(y, estimator.spec)
        residual = det_H0 - det(H)

        self.print(f"Selected determinant: {det_H0.change_ring(RDF)}")
        self.print(f"Estimated epsilon: {float(eps)}")
        self.print(
            f"Determinant approximation L2 distance: {sqrt(float(polynomial_ss(residual)))}"
        )

        _, x_basis = estimator.build_solution_space()
        self.print(f"Kernel DoF: {x_basis.ncols()}")

        estimator.solve()

        x_residual = estimator.x - x_initial
        self.print(
            f"Row approximation L2 distance: {sqrt(float(vector_ss(x_residual)))}"
        )

        self.print("Selected row projection")
        return estimator.reconstruct()

    def factorize(
        self,
        P,
        candidate_selector: CandidateSelector | None = None,
        beam_width: int = 16,
    ):
        factorizer = Factorizer(P)
        factorizer.normalize()

        if candidate_selector is None:
            _, gcd = factorizer.euclidean()
            self.print("Euclidean pre-factorization computed")
            self.print(f"GCD degree: {min_degree(gcd)}")

            P0 = factorizer.reconstruct()
            self.print(f"P0 matrix reconstructed with det P0 = {det(P0)}")

            factorizer.recover()
            self.print("Recovery step completed")
            return factorizer.pack()

        candidates = factorizer.build_candidates(beam_width=beam_width)
        self.print(f"Evaluated {len(candidates)} factorization candidates")
        return candidate_selector(candidates)

    def process(
        self,
        low_pass,
        high_pass,
        F=QQ,
        projection_strategy: ProjectionStrategy = ProjectionStrategy.FIXED_ROW,
        candidate_selector: CandidateSelector | None = None,
        beam_width: int = 16,
    ):
        if len(low_pass) != len(high_pass):
            raise ValueError("Analysis filters must have the same length")

        R = LaurentPolynomialRing(F, names=("z",))

        P = build_polyphase_matrix(low_pass, high_pass, R)

        P_est = self.estimate(P, strategy=projection_strategy)

        steps, delays = self.factorize(
            P_est,
            candidate_selector=candidate_selector,
            beam_width=beam_width,
        )

        P_rec = reconstruct(steps, delays, R)

        assert (
            P_rec == P_est
        ), "Reconstructed matrix does not match the estimated matrix"

        self.print("Sanity check: OK")

        l2 = sqrt(float(matrix_ss(P, P_rec)))
        self.print(f"Reconstruction L2 distance in QQ mode: {l2}")

        R_fp64 = LaurentPolynomialRing(RDF, names=("z",))
        steps_fp64 = [(q.change_ring(RDF), step) for q, step in steps]
        P_rec_fp64 = reconstruct(steps_fp64, delays, R_fp64)
        reference_fp64 = matrix(
            R_fp64,
            2,
            2,
            [polynomial.change_ring(RDF) for polynomial in P.list()],
        )

        l2_fp64 = sqrt(float(matrix_ss(reference_fp64, P_rec_fp64, F=RDF)))

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
