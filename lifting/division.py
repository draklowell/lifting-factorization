from sage.all import vector

from lifting.utils import (
    degree,
    from_vector,
    laurent_toeplitz,
    max_degree,
    min_degree,
    to_vector,
)


def try_solve(A, y):
    try:
        v_p = A.solve_right(y)
    except ValueError:
        return None

    return v_p


def find_solution(A, y, direction: bool = False):
    """
    Search for the longest symmetric solution
    (i.e. a solution that satisfies Ax=y for the most
    number of points from the left and from the right)

    If direction = True => first left then right
    otherwise first right then left
    """
    x = None
    row_count = A.nrows()
    for span in range(1, row_count + 1):
        shorter = span // 2
        longer = span - shorter

        head = longer if direction else shorter
        tail = row_count - (shorter if direction else longer)

        # Fix needed: memory-heavy task copying and stacking
        # matrices
        A_cur = A[:head].stack(A[tail:])
        y_cur = vector(tuple(y[:head]) + tuple(y[tail:]))

        x_candidate = try_solve(A_cur, y_cur)
        if x_candidate is None:
            return x

        x = x_candidate

    return x


def ldivmod(a, b):
    """
    Perform polynomial division of a by b, that preserves symmetry,
    returning the quotient and remainder.
    """

    if degree(b) == -1:
        raise ZeroDivisionError("Cannot divide by zero polynomial.")

    if degree(a) < degree(b):
        return 0, a

    R = a.parent()

    B = laurent_toeplitz(
        b,
        degree(a) + 1,
        degree(a) - degree(b) + 1,
        min_degree(b),
    )
    a_v = to_vector(a, min_degree(a), degree(a) + 1)

    # Start eliminating from the
    # longer side to preserve symmetry
    direction = -min_degree(a) > max_degree(a)

    q_v = find_solution(B, a_v, direction)
    if q_v is None:
        raise ValueError("No solution found for division")

    r_v = a_v - B * q_v
    r = from_vector(r_v, min_degree(a), R)
    q = from_vector(q_v, min_degree(a) - min_degree(b), R)

    assert a == b * q + r, "Division result is incorrect: a != b*q + r"
    assert degree(r) < degree(
        b
    ), "Remainder is not smaller than the divisor: degree(r) >= degree(b)"

    return q, r


def ldivmod_candidates(a, b):
    """
    Return exact Laurent divisions obtained from every head/tail split.
    """
    if degree(b) == -1:
        raise ZeroDivisionError("Cannot divide by zero polynomial.")

    if degree(a) < degree(b):
        return [(a.parent()(0), a)]

    R = a.parent()
    row_count = degree(a) + 1
    quotient_size = degree(a) - degree(b) + 1
    B = laurent_toeplitz(
        b,
        row_count,
        quotient_size,
        min_degree(b),
    )
    a_v = to_vector(a, min_degree(a), row_count)
    candidates = []

    for head in range(quotient_size + 1):
        tail = quotient_size - head
        tail_begin = row_count - tail
        A_cur = B[:head].stack(B[tail_begin:])
        y_cur = vector(tuple(a_v[:head]) + tuple(a_v[tail_begin:]))
        q_v = try_solve(A_cur, y_cur)
        if q_v is None:
            continue

        r_v = a_v - B * q_v
        r = from_vector(r_v, min_degree(a), R)
        q = from_vector(q_v, min_degree(a) - min_degree(b), R)
        if degree(r) >= degree(b):
            continue
        if any(q == existing_q for existing_q, _ in candidates):
            continue
        assert a == b * q + r
        candidates.append((q, r))

    if not candidates:
        raise ValueError("No valid Laurent division candidate found")
    return candidates


def ldiv(a, b):
    """
    Perform polynomial division of a by b, returning the quotient.
    """
    q, r = ldivmod(a, b)
    if degree(r) >= 0:
        raise ValueError(f"Division is not exact, non-zero remainder: r = {r} != 0")

    return q
