from sage.all import *


def min_degree(p):
    coeffs = p.monomial_coefficients()
    if not coeffs:
        return 0
    return min(coeffs.keys())


def max_degree(p):
    coeffs = p.monomial_coefficients()
    if not coeffs:
        return 0
    return max(coeffs.keys())


def degree(p):
    deg = max_degree(p) - min_degree(p)
    if deg == 0 and p == 0:
        return -1

    return deg


def laurent_toeplitz(p, rows, cols, offset: int = 0):
    T = matrix(p.base_ring(), rows, cols)
    coeffs = p.monomial_coefficients()
    for m, c in coeffs.items():
        for i in range(rows):
            j = offset + i - m
            if 0 <= j < cols:
                T[i, j] = c
    return T


def to_vector(p, offset, size):
    v = vector(p.base_ring(), size)
    for m, c in p.monomial_coefficients().items():
        i = m - offset
        if i < 0 or i >= size:
            raise ValueError(
                f"Monomial z^{m} is out of bounds for the output vector representation."
            )

        v[i] = c
    return v


def from_vector(v, offset, base_ring):
    p = base_ring(0)
    z = p.parent().gen()
    for i in range(len(v)):
        p += v[i] * z ** (i + offset)
    return p
