from sage.all import *
from lifting.matrix_estimator import MatrixEstimator
from lifting.factorizer import Factorizer, LiftingStep
from lifting.utils import max_degree, min_degree
from pywt import Wavelet
import json


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

class Processor:
    def __init__(self, output = None):
        self.output = output

    def process(self, name: str, output_path: str, F=QQ):
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

        err = (P_est - P_err).list()
        l2 = 0
        for i in err:
            l2 += sum(c**2 for c in i)

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

        for q, step in steps:
            if step == LiftingStep.PREDICT:
                P_rec *= matrix(R, [[1, q], [0, 1]])
            elif step == LiftingStep.UPDATE:
                P_rec *= matrix(R, [[1, 0], [q, 1]])
            elif step == LiftingStep.SCALE:
                P_rec *= matrix(R, [[q, 0], [0, 1/q]])

        assert P_rec == P_est, "Reconstructed matrix does not match the estimated matrix"

        output.write(f"Sanity check: OK\n")

        result = {
            "tap_size": len(wavelet.dec_hi),
            "approximation_l2": float(l2),
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
                LiftingStep.SCALE: "scale",
            }[step]

            q_vector = []
            for i in range(min_degree(q), max_degree(q)+1):
                q_vector.append(float(q.coefficient(i)))

            result["steps"].append({
                "type": step_name,
                "shift": min_degree(q),
                "coefficients": q_vector,
            })

        with open(output_path, "w") as file:
            json.dump(result, file, indent=4)
