import json
import sys

import pywt

from processor.processor import Processor

WAVELETS = [
    "bior4.4",
    "bior2.2",
    "bior1.3",
    "coif1",
    "coif2",
    "sym8",
    "haar",
    "db20",
] + [f"db{i}" for i in range(2, 13)]


def main() -> None:
    for name in WAVELETS:
        print(f"### {name} ###")
        wavelet = pywt.Wavelet(name)

        processor = Processor(sys.stdout)
        result = processor.process(wavelet.dec_lo, wavelet.dec_hi)

        with open(f"coeffs/{name}.json", "w") as file:
            json.dump(result, file)


if __name__ == "__main__":
    main()
