import os
import sys

import numpy as np
from pywt import Wavelet, dwt

import dtypes
from lifting import LiftingScheme

sys.set_int_max_str_digits(1000000)


def print_cmp(a, b, name):
    err = np.abs(a - b)
    err_mse = np.sum(err**2) / len(err)
    message = " " + name.upper() + ": "
    message = message + " " * (17 - len(message))

    if err_mse < 5e-6:
        message += f"\033[32m{err_mse:.2e}\033[0m"
    elif err_mse < 1e-2:
        message += f"\033[33m{err_mse:.2e}\033[0m"
    else:
        message += f"\033[31m{err_mse:.2e}\033[0m"

    print(message)


counter = 0
for filename in os.listdir("../coeffs/"):
    wavelet = filename.removesuffix(".json")
    print(f"WAVELET: {wavelet}")

    wavelet = Wavelet(wavelet)

    scheme = LiftingScheme.from_file(
        f"../coeffs/{filename}", "symmetric", dtype=dtypes.float64
    )

    test_signal = np.random.normal(scale=1, size=(1001,))

    approx_lifting, details_lifting = scheme.forward(test_signal)
    approx_dwt, details_dwt = dwt(test_signal, wavelet=wavelet, mode="symmetric")

    print_cmp(approx_dwt, approx_lifting, "approximation")
    print_cmp(details_dwt, details_lifting, "details")

    data_lifting = scheme.inverse(approx_lifting, details_lifting)
    data_lifting = data_lifting[: len(test_signal)]

    print_cmp(test_signal, data_lifting, "reconstruction")
    counter += 1

print(f"Tested {counter} wavelets")
