from sage.all import *
from lifting.utils import (
    from_vector,
    laurent_toeplitz,
    min_degree,
    to_vector,
    degree,
)


def try_solve(A, y):
    try:
        # print(A)
        # print(y)
        v_p = A.solve_right(y)
    except ValueError:
        return None

    return v_p

# Search for the longest suffix solution
# $y = Ax$
def find_solution(A, y):
    for i in range(A.nrows()):
        x_candidate = try_solve(A[i:], y[i:])
        if x_candidate is not None:
            return x_candidate
    
    return None
    # low = 0
    # high = A.nrows() - 1

    # while low < high:
    #     mid = (low + high) // 2
    #     x_candidate = try_solve(A[mid:], y[mid:])
    #     if x_candidate is not None:
    #         high = mid
    #     else:
    #         low = mid + 1

    # x = try_solve(A[low:], y[low:])

    # return x

def ldivmod(a, b):
    if degree(b) == -1:
        raise ZeroDivisionError("Cannot divide by zero polynomial.")

    if degree(a) < degree(b):
        return 0, a

    R = a.parent()

    B = laurent_toeplitz(b, degree(a) + 1, degree(a) - degree(b) + 1, min_degree(b))
    a_v = to_vector(a, min_degree(a), degree(a) + 1)

    q_v = find_solution(B, a_v)
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