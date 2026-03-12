from pywt import Wavelet
from sage.all import *

from lifting.matrix_estimator import MatrixEstimator
from lifting.factorizer import factorize

F = QQ
L = LaurentPolynomialRing(F, names=("z",))
(z,) = L.gens()

db4 = Wavelet("sym6")

# low pass
h = db4.rec_lo
he = 0
ho = 0
for i in range(len(h)):
    if i % 2 == 0:
        he += F(h[i]) * z ** Integer(i // 2)
    else:
        ho += F(h[i]) * z ** Integer(i // 2)

# high pass
g = db4.rec_hi
ge = 0
go = 0
for i in range(len(g)):
    if i % 2 == 0:
        ge += F(g[i]) * z ** Integer(i // 2)
    else:
        go += F(g[i]) * z ** Integer(i // 2)

# polyphase reconstruction (synthesis) matrix
H = matrix(
    L,
    [
        [he, ge],
        [ho, go],
    ],
)

H_new = MatrixEstimator.solve(H, normalize_to=-1)
d = det(H_new)
delay = d.degree()
# print(f"Determinant: {d} (delay: {delay})")

delay_up = delay // 2
delay_down = delay - delay_up

H_norm = matrix(
    L,
    [
        [z ** -delay_up, 0],
        [0, -z ** -delay_down],
    ],
) * H_new
print(f"Determinant after normalization: {det(H_norm)}")


factorize(H_norm)
