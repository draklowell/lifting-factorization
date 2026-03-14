from sage.all import *

from lifting.factorizer import Factorizer, LiftingStep
from lifting.matrix_estimator import MatrixEstimator
from lifting.utils import max_degree, min_degree
from pywt import Wavelet


class ProxyOutput:
    def __init__(self, prefix: str, output=None):
        self.prefix = prefix
        self.output = output

    def write(self, data: str):
        for line in data.splitlines():
            self.output.write(self.prefix + line + "\n")

def split(filter_coeffs: list[float], R):
    F = R.base_ring()
    (z,) = R.gens()

    even = 0
    odd = 0
    for i, coeff in enumerate(filter_coeffs):
        if i % 2 == 0:
            even += F(coeff) * z ** Integer(i // 2)
        else:
            odd += F(coeff) * z ** Integer(i // 2)

    return even, odd

def normalize(P, R, output=None):
    (z,) = R.gens()
    d = det(P)
    delay = d.degree()
    delay_up = delay // 2
    delay_down = delay - delay_up

    P_norm = matrix(
        R,
        [
            [z ** -delay_up, 0],
            [0, z ** -delay_down],
        ],
    ) * P

    if output is not None:
        output.write(f"Determinant after meta-normalization: {det(P_norm)}\n")

    return P_norm, (delay_up, delay_down)

def get_l2(a, b):
    """Compute coefficient-wise L2 distance between two 2x2 polynomial matrices."""
    accum = 0
    for poly_a, poly_b in zip(a.list(), b.list()):
        coeffs_a = poly_a.dict()
        coeffs_b = poly_b.dict()
        for m in set(coeffs_a.keys()) | set(coeffs_b.keys()):
            c_a = coeffs_a.get(m, 0)
            c_b = coeffs_b.get(m, 0)
            accum += (c_a - c_b) ** 2

    return sqrt(float(accum))

class Processor:
    STEP_NAMES = {
        LiftingStep.PREDICT: "predict",
        LiftingStep.UPDATE: "update",
        LiftingStep.SCALE_EVEN: "scale-even",
        LiftingStep.SCALE_ODD: "scale-odd",
        LiftingStep.SWAP: "swap",
    }

    def __init__(self, output=None):
        self.output = output

    @staticmethod
    def _build_error_matrix(wavelet, ring):
        he, ho = split(wavelet.dec_lo, ring)
        ge, go = split(wavelet.dec_hi, ring)
        return matrix(ring, [[he, ge], [ho, go]])

    @staticmethod
    def _delay_matrix(ring, delay_up, delay_down):
        (z,) = ring.gens()
        return matrix(
            ring,
            [
                [z ** delay_up, 0],
                [0, z ** delay_down],
            ],
        )

    @staticmethod
    def _step_matrix(ring, q, step):
        if step == LiftingStep.PREDICT:
            return matrix(ring, [[1, q], [0, 1]])
        if step == LiftingStep.UPDATE:
            return matrix(ring, [[1, 0], [q, 1]])
        if step == LiftingStep.SCALE_EVEN:
            return matrix(ring, [[q, 0], [0, 1]])
        if step == LiftingStep.SCALE_ODD:
            return matrix(ring, [[1, 0], [0, q]])
        if step == LiftingStep.SWAP:
            return matrix(ring, [[0, 1], [1, 0]])
        raise ValueError(f"Unknown lifting step: {step}")

    @staticmethod
    def _serialize_coeffs(q):
        coeffs = []
        for degree_idx in range(min_degree(q), max_degree(q) + 1):
            coeff = q.coefficient(degree_idx)
            coeffs.append(
                {
                    "numerator": int(coeff.numerator()),
                    "denominator": int(coeff.denominator()),
                }
            )

        if coeffs == [0]:
            return []

        return coeffs

    def _reconstruct(self, steps, ring, delay_up, delay_down):
        rec = matrix(ring, 2, 2, [[1, 0], [0, 1]])
        rec *= self._delay_matrix(ring, delay_up, delay_down)
        for q, step in steps:
            rec *= self._step_matrix(ring, q, step)
        return rec

    def _to_result(self, wavelet, steps, delay_up, delay_down, l2_fp64, l2):
        result = {
            "tap_size": len(wavelet.dec_hi),
            "l2": {
                "fp64": l2_fp64,
                "qq": l2,
            },
            "delay": {
                "even": int(delay_up),
                "odd": int(delay_down),
            },
            "steps": [],
        }

        for q, step in steps:
            q_vector = self._serialize_coeffs(q)
            if len(q_vector) == 0 and step in {LiftingStep.PREDICT, LiftingStep.UPDATE}:
                continue

            result["steps"].append(
                {
                    "type": self.STEP_NAMES[step],
                    "shift": min_degree(q),
                    "coefficients": q_vector,
                }
            )

        return result

    def process(self, name: str, F=QQ):
        wavelet = Wavelet(name)

        R = LaurentPolynomialRing(F, names=("z",))

        output = ProxyOutput(f"[{name}] ", self.output)

        P_err = self._build_error_matrix(wavelet, R)

        P_est = MatrixEstimator.solve(P_err, normalize_to=1, output=output)

        P, (delay_up, delay_down) = normalize(P_est, R, output)

        steps = Factorizer.factorize(P, output=output)

        # Sanity check
        P_rec = self._reconstruct(steps, R, delay_up, delay_down)

        R_fp64 = LaurentPolynomialRing(RR, names=("z",))
        steps_fp64 = [(q.change_ring(RR), step) for q, step in steps]
        P_rec_fp64 = self._reconstruct(steps_fp64, R_fp64, delay_up, delay_down)

        assert P_rec == P_est, "Reconstructed matrix does not match the estimated matrix"
        output.write(f"Sanity check: OK\n")

        l2_fp64 = get_l2(P_err, P_rec_fp64)
        l2 = get_l2(P_err, P_rec)

        output.write(f"Reconstruction L2 distance in FP64 mode: {l2_fp64}\n")
        return self._to_result(wavelet, steps, delay_up, delay_down, l2_fp64, l2)
