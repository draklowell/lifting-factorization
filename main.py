import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

from processor import Processor

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


def run_one(name: str) -> None:
    processor = Processor(sys.stdout)
    result = processor.process(name)

    with open(f"coeffs/{name}.json", "w") as file:
        json.dump(result, file)


def main() -> None:
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(run_one, name) for name in WAVELETS]
        for future in as_completed(futures):
            future.result()


if __name__ == "__main__":
    main()
