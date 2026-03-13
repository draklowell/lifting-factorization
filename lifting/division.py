from sage.all import *
from lifting.utils import (
    from_vector,
    laurent_toeplitz,
    min_degree,
    max_degree,
    to_vector,
    degree,
)


def try_solve(A, y):
    try:
        v_p = A.solve_right(y)
    except ValueError:
        return None

    return v_p

# Search for the longest suffix solution
# $y = Ax$
# If direction = True => first left then right
# otherwise first right then left
def find_solution(A, y, direction: bool = False):
    x = None
    for l in range(1, A.nrows()+1):
        shorter = l // 2
        longer = l - shorter

        s = longer if direction else shorter
        e = A.nrows() - (shorter if direction else longer)

        A_cur = A[:s].stack(A[e:])
        y_cur = vector(tuple(y[:s]) + tuple(y[e:]))

        x_candidate = try_solve(A_cur, y_cur)
        if x_candidate is None:
            return x
        
        x = x_candidate
    
    return x

def ldivmod(a, b):
    if degree(b) == -1:
        raise ZeroDivisionError("Cannot divide by zero polynomial.")

    if degree(a) < degree(b):
        return 0, a

    R = a.parent()

    B = laurent_toeplitz(b, degree(a) + 1, degree(a) - degree(b) + 1, min_degree(b))
    a_v = to_vector(a, min_degree(a), degree(a) + 1)

    direction = -min_degree(a) > max_degree(a)

    q_v = find_solution(B, a_v, direction)
    if q_v is None:
        raise ValueError("No solution found for division")

    r_v = a_v - B * q_v
    r = from_vector(r_v, min_degree(a), R)
    q = from_vector(q_v, min_degree(a) - min_degree(b), R)
    
    assert a == b * q + r, "Division result is incorrect: a != b*q + r"
    assert degree(r) < degree(b), "Remainder is not smaller than the divisor: degree(r) >= degree(b)"

    return q, r

def ldiv(a, b):
    q, r = ldivmod(a, b)
    if degree(r) >= 0:
        raise ValueError("Division is not exact, non-zero remainder: r != 0")

    return q