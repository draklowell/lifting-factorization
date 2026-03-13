from sage.all import *
import math
from lifting.utils import (
    degree,
    min_degree
)
from lifting.division import (
    ldiv,
    ldivmod,
)
from enum import Enum

class LiftingStep(Enum):
    PREDICT = 1
    UPDATE = 2
    SCALE_EVEN = 3
    SCALE_ODD = 4
    SWAP = 5


class Factorizer:
    def __init__(self, output = None):
        self.output = output

    def print(self, *args):
        if self.output is not None:
            self.output.write(" ".join(map(str, args)) + "\n")

    def euclidean(self, a, b):
        result = []
        if degree(a) < degree(b):
            a, b = b, a
            result.append(0) # Swap term

        # $a_{i+1} = b_i,\quad b_{i+1} = a_i \bmod b_i, \quad q_{i+1} = a_i / b_i$
        while True:
            q, r = ldivmod(a, b)

            if degree(r) >= degree(b):
                raise ValueError("Remainder degree must be less than divisor degree")

            a = b
            b = r       
            result.append(q)
            self.print("Euclidean step computed:")
            self.print(f" Degrees: |q|={degree(q)}, |r|={degree(r)}, |a|-|b|={degree(a) - degree(b)}")
            l2q = sum(c**2 for c in q.coefficients())
            self.print(f" L2 norm: ||q||={math.sqrt(l2q)}")

            if degree(r) == -1:
                break

        return result, a

    def _factorize(self, P):
        if det(P) != 1:
            raise ValueError("Matrix must have determinant 1")

        R = P.base_ring()

        he, ho, ge, go = P[0, 0], P[1, 0], P[0, 1], P[1, 1]
        qs, a = self.euclidean(he, ho)

        self.print("Euclidean pre-factorization computed")

        self.print(f"GCD degree: {min_degree(a)}")

        #
        # $P^{(0)}_0 = \begin{bmatrix} K & 0 \\ 0 & 1/K \end{bmatrix}$
        #
        P0 = matrix(R, 2, 2, [[a, 0], [0, 1/a]])
        for q in reversed(qs):
            #
            # $P^{(0)}_{i+1} = \begin{bmatrix} q_i & 1 \\ 1 & 0 \end{bmatrix} P^{(0)}_i$
            #
            Q = matrix(R, 2, 2, [[q, 1], [1, 0]])
            P0 = Q * P0

        self.print("P0 matrix reconstructed")

        # $P^{(0)} = P^{(0)}_{n}$

        #
        # $P = P^{(0)} \begin{bmatrix} 1 & s \\ 0 & 1 \end{bmatrix}$
        #

        #
        # $\begin{bmatrix} h_e & g_e \\ h_o & g_o \end{bmatrix} = \begin{bmatrix} h_e & g_e^{(0)} \\ h_o & g_o^{(0)} \end{bmatrix} \begin{bmatrix} 1 & s \\ 0 & 1 \end{bmatrix}$
        #

        # $g_e = h_e \cdot s + g_e^{(0)}$
        # $h_e \cdot s = g_e - g_e^{(0)}$
        #
        # $s = \frac{g_e - g_e^{(0)}}{h_e}$
        #
        # And same for odd part, thay have to be the same
        d = det(P0)
        self.print("Determinant of a reconstructed P0:", d)
        b = 1/d

        ge0, go0 = P0[0][1], P0[1][1]
        se = ldiv(ge - b*ge0, he)
        so = ldiv(go - b*go0, ho)
        assert se == so
        self.print("Recovering coefficient S computed")

        return qs, a, se, b

    @staticmethod
    def factorize(P, output=None) -> list[tuple[object, LiftingStep]]:
        factorizer = Factorizer(output)
        qs, a, s, b = factorizer._factorize(P)

        steps = []
        for i, q in enumerate(qs):
            step_type = LiftingStep.PREDICT if i % 2 == 0 else LiftingStep.UPDATE
            if q == 0:
                continue

            steps.append((q, step_type))

        if len(qs) % 2 == 1:
            R = P.base_ring()
            if output is not None:
                output.write("Adding swap step to handle odd number of Euclidean steps\n")

            steps.append((R(0), LiftingStep.SWAP))

        steps.append((a*a*s/b, LiftingStep.PREDICT))
        steps.append((a, LiftingStep.SCALE_EVEN))
        steps.append((b/a, LiftingStep.SCALE_ODD))

        return steps
