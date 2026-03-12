from sage.all import *
import math
from lifting.utils import (
    degree,
    min_degree,
)
from lifting.division import (
    ldiv,
    ldivmod,
)


def euclidean(a, b):
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
        print("Euclidean step computed:")
        print(f" Degrees: |q|={degree(q)}, |r|={degree(r)}, |a|-|b|={degree(a) - degree(b)}")
        l2q = sum(c**2 for c in q.coefficients())
        print(f" L2 norm: ||q||={math.sqrt(l2q)}")

        if degree(r) == -1:
            break

    return result, a

def factorize(P):
    if det(P) != 1:
        raise ValueError("Matrix must have determinant 1")

    R = P.base_ring()

    he, ho, ge, _ = P[0, 0], P[1, 0], P[0, 1], P[1, 1]
    qs, a = euclidean(he, ho)
    print("Euclidean pre-factorization computed")

    print(f"GCD: {a}")

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
    print("P0 matrix reconstructed")

    # $P^{(0)} = P^{(0)}_{n}$

    #
    # $P = P^{(0)} \begin{bmatrix} 1 & s \\ 0 & 1 \end{bmatrix}$
    #

    #
    # $\begin{bmatrix} h_e & g_e \\ h_o & g_o \end{bmatrix} = \begin{bmatrix} h_e & g_e^{(0)} \\ h_o & g_o^{(0)} \end{bmatrix} \begin{bmatrix} 1 & s \\ 0 & 1 \end{bmatrix}$
    #

    # We need to normalize by b, because this factorization
    # is not consistent with Sweldens et al., we do not assume
    # $\det P = 1$

    #
    # $\begin{bmatrix} h_e & g_e \\ h_o & g_o \end{bmatrix} = \begin{bmatrix} h_e & g_e^{(0)} \\ h_o & g_o^{(0)} \end{bmatrix} \begin{bmatrix} 1 & s \\ 0 & b \end{bmatrix}$
    #

    # $g_e = h_e \cdot s + b \cdot g_e^{(0)}$
    # $h_e \cdot s = g_e - b\cdot g_e^{(0)}$
    #
    # $s = \frac{g_e - b\cdot g_e^{(0)}}{h_e}$
    #
    # And same for odd part, thay have to be the same

    ge0 = P0[0][1]
    s = ldiv(ge - ge0, he)
    print("s computed")

    return qs, a, s
