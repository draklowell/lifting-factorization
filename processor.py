from sage.all import *
from lifting.matrix_estimator import MatrixEstimator
from lifting.factorizer import Factorizer, LiftingStep
from lifting.utils import max_degree, min_degree
from pywt import Wavelet


class ProxyOutput:
    def __init__(self, prefix: str, output = None):
        self.prefix = prefix
        self.output = output

    def write(self, data: str):
        for line in data.splitlines():
            self.output.write(self.prefix + line + "\n")

def split(filter: list[float], R):
    F = R.base_ring()
    (z,) = R.gens()

    even = 0
    odd = 0
    for i, coeff in enumerate(filter):
        if i % 2 == 0:
            even += F(coeff) * z ** Integer(i // 2)
        else:
            odd += F(coeff) * z ** Integer(i // 2)

    return even, odd

def normalize(P, R, output = None):
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
    def __init__(self, output = None):
        self.output = output

    def process(self, name: str, F=QQ):
        wavelet = Wavelet(name)

        R = LaurentPolynomialRing(F, names=("z",))
        (z,) = R.gens()

        output = ProxyOutput(f"[{name}] ", self.output)

        he, ho = split(wavelet.dec_lo, R)
        ge, go = split(wavelet.dec_hi, R)

        P_err = matrix(
            R,
            [
                [he, ge],
                [ho, go],
            ],
        )

        P_est = MatrixEstimator.solve(P_err, normalize_to=1, output=output)

        P, (delay_up, delay_down) = normalize(P_est, R, output)

        steps = Factorizer.factorize(P, output=output)

        # Sanity check
        P_rec = matrix(R, 2, 2, [[1, 0], [0, 1]])
        P_rec *= matrix(
            R,
            [
                [z ** delay_up, 0],
                [0, z ** delay_down],
            ],
        )

        R_fp64 = LaurentPolynomialRing(RR, names=("z",))
        (z_fp64,) = R_fp64.gens()
        P_rec_fp64 = matrix(R_fp64, 2, 2, [[1, 0], [0, 1]])
        P_rec_fp64 *= matrix(
            R_fp64,
            [
                [z_fp64 ** delay_up, 0],
                [0, z_fp64 ** delay_down],
            ],
        )

        for q, step in steps:
            q_fp64 = q.change_ring(RR)
            if step == LiftingStep.PREDICT:
                P_rec *= matrix(R, [[1, q], [0, 1]])
                P_rec_fp64 *= matrix(R_fp64, [[1, q_fp64], [0, 1]])
            elif step == LiftingStep.UPDATE:
                P_rec *= matrix(R, [[1, 0], [q, 1]])
                P_rec_fp64 *= matrix(R_fp64, [[1, 0], [q_fp64, 1]])
            elif step == LiftingStep.SCALE_EVEN:
                P_rec *= matrix(R, [[q, 0], [0, 1]])
                P_rec_fp64 *= matrix(R_fp64, [[q_fp64, 0], [0, 1]])
            elif step == LiftingStep.SCALE_ODD:
                P_rec *= matrix(R, [[1, 0], [0, q]])
                P_rec_fp64 *= matrix(R_fp64, [[1, 0], [0, q_fp64]])
            elif step == LiftingStep.SWAP:
                P_rec *= matrix(R, [[0, 1], [1, 0]])
                P_rec_fp64 *= matrix(R_fp64, [[0, 1], [1, 0]])
            else:
                raise ValueError(f"Unknown lifting step: {step}")

        assert P_rec == P_est, "Reconstructed matrix does not match the estimated matrix"
        output.write(f"Sanity check: OK\n")

        l2_fp64 = get_l2(P_err, P_rec_fp64)
        l2 = get_l2(P_err, P_rec)

        output.write(f"Reconstruction L2 distance in FP64 mode: {l2_fp64}\n")

        result = {
            "tap_size": len(wavelet.dec_hi),
            "l2": {
                "fp64": l2_fp64,
                "qq": l2,
            },
            "delay": {
                "even": int(delay_up),
                "odd": int(delay_down)
            },
            "steps": [],
        }

        for q, step in steps:
            step_name = {
                LiftingStep.PREDICT: "predict",
                LiftingStep.UPDATE: "update",
                LiftingStep.SCALE_EVEN: "scale-even",
                LiftingStep.SCALE_ODD: "scale-odd",
                LiftingStep.SWAP: "swap",
            }[step]

            q_vector = []
            for i in range(min_degree(q), max_degree(q)+1):
                coeff = q.coefficient(i)
                q_vector.append({
                    "numerator": int(coeff.numerator()),
                    "denominator": int(coeff.denominator())
                })

            if q_vector == [0]:
                q_vector = []

            if len(q_vector) == 0 and step in {LiftingStep.PREDICT, LiftingStep.UPDATE}:
                continue

            result["steps"].append({
                "type": step_name,
                "shift": min_degree(q),
                "coefficients": q_vector,
            })

        return result
