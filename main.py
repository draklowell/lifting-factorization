import argparse
import json
import sys
from pathlib import Path

from lifting.projection import ProjectionStrategy
from processor.processor import Processor


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Factor an analysis FIR filter bank into lifting steps."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="JSON file containing low_pass and high_pass coefficient arrays",
    )
    parser.add_argument("--output", type=Path, help="output lifting JSON")
    parser.add_argument(
        "--projection",
        choices=[strategy.value for strategy in ProjectionStrategy],
        default=ProjectionStrategy.FIXED_ROW.value,
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    filter_bank = json.loads(arguments.input.read_text())
    try:
        low_pass = filter_bank["low_pass"]
        high_pass = filter_bank["high_pass"]
    except KeyError as error:
        raise ValueError(f"Missing filter-bank field: {error.args[0]}") from error

    scheme = Processor(sys.stderr).process(
        low_pass,
        high_pass,
        projection_strategy=ProjectionStrategy(arguments.projection),
    )
    encoded = json.dumps(scheme, indent=2) + "\n"
    if arguments.output is None:
        sys.stdout.write(encoded)
    else:
        arguments.output.write_text(encoded)


if __name__ == "__main__":
    main()
