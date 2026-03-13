from pywt import wavelist
from concurrent.futures import ProcessPoolExecutor

from processor import Processor
import sys

#wavelets = wavelist(kind="discrete")
wavelets = [
    "bior4.4",
    "bior2.2",
    "bior1.3",
    "coif1",
    "coif2",
    "sym8",
    "haar",
    "db20",
] + [f"db{i}" for i in range(2, 13)]

def run_one(name):
    processor = Processor(sys.stdout)
    return processor.process(name, f"coeffs/{name}.json")

with ProcessPoolExecutor() as ex:
    futures = []
    for name in wavelets:
        futures.append(ex.submit(run_one, name))

    for f in futures:
        f.result()
