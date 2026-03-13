from pywt import wavelist
from concurrent.futures import ProcessPoolExecutor

from processor import Processor
import sys

wavelet_name = "coif4"

def run_one(name):
    processor = Processor(sys.stdout)
    return processor.process(name, f"coeffs/{name}.json")

with ProcessPoolExecutor() as ex:
    futures = []
    for name in wavelist(kind="discrete"):
        if name != wavelet_name:
            continue

        # if name.startswith("db"):
        #     num = int(name[2:])

        #     # Skip large db's
        #     if num >= 20:
        #         continue

        futures.append(ex.submit(run_one, name))

    for f in futures:
        f.result()
