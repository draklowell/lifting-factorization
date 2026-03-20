import os
import sys

import numpy as np
from pywt import Wavelet, dwt

import dtypes
from lifting import DType, LiftingScheme

sys.set_int_max_str_digits(1000000)


def print_cmp(cmp, ref, name):
    err = np.abs(cmp - ref)

    l2_err = np.sqrt((np.sum(err**2)) / np.sum(ref**2))
    # linf_err = np.max(err / np.abs(ref))

    rel_err = l2_err

    message = " " + name.upper() + ": "
    message = message + " " * (17 - len(message))

    if rel_err < 1e-3:
        message += f"\033[32m{rel_err:.2e}\033[0m"
    elif rel_err < 5e-3:
        message += f"\033[33m{rel_err:.2e}\033[0m"
    else:
        message += f"\033[31m{rel_err:.2e}\033[0m"

    print(message)


counter = 0
for filename in os.listdir("../coeffs-new/"):
    if not filename.endswith("-fp64.json"):
        continue

    wavelet = filename.removesuffix("-fp64.json")
    print(f"WAVELET: {wavelet}")

    wavelet = Wavelet(wavelet)

    scheme = LiftingScheme.from_file(
        f"../coeffs-new/{filename}",
        "symmetric",
        predict_dtype=DType(dtypes.mixed_tf32xf32, dtypes.mixed_tf32xf32),
        update_dtype=DType(dtypes.mixed_tf32xf32, dtypes.mixed_tf32xf32),
        scale_dtype=DType(dtypes.mixed_tf32xf32, dtypes.mixed_tf32xf32),
    )

    test_signal = np.random.normal(scale=1, size=(1001,))

    approx_lifting, details_lifting = scheme.forward(test_signal)
    approx_dwt, details_dwt = dwt(test_signal, wavelet=wavelet, mode="symmetric")

    print_cmp(approx_lifting, approx_dwt, "approximation")
    print_cmp(details_lifting, details_dwt, "details")

    data_lifting = scheme.inverse(approx_lifting, details_lifting)
    data_lifting = data_lifting[: len(test_signal)]

    print_cmp(data_lifting, test_signal, "reconstruction")
    counter += 1

print(f"Tested {counter} wavelets")
