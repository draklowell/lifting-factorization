import re
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import cache
from pathlib import Path

ARRAY_PATTERN = re.compile(
    r"static const TYPE CAT\((?P<name>[A-Za-z0-9_]+)_, TYPE\)"
    r"\[(?P<size>\d+)\]\s*=\s*\{(?P<body>.*?)\}\s*;",
    re.DOTALL,
)
NUMBER_PATTERN = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?")
SQRT2_PATTERN = re.compile(
    r"static const TYPE CAT\(sqrt2_, TYPE\)\s*=\s*(?P<value>[^;]+);"
)


@dataclass(frozen=True)
class FilterBank:
    dec_lo: tuple[Fraction, ...]
    dec_hi: tuple[Fraction, ...]
    rec_lo: tuple[Fraction, ...]
    rec_hi: tuple[Fraction, ...]
    orthogonal: bool


def decimal_fraction(value: str | Decimal) -> Fraction:
    numerator, denominator = Decimal(value).as_integer_ratio()
    return Fraction(numerator, denominator)


@cache
def parse_coefficient_arrays(path: Path) -> dict[str, tuple[Fraction, ...]]:
    arrays = {}
    for match in ARRAY_PATTERN.finditer(Path(path).read_text()):
        values = tuple(
            decimal_fraction(value) for value in NUMBER_PATTERN.findall(match["body"])
        )
        size = int(match["size"])
        if len(values) != size:
            raise ValueError(
                f"{match['name']} declares {size} coefficients but contains "
                f"{len(values)}"
            )
        arrays[match["name"]] = values
    return arrays


@cache
def parse_sqrt2(path: Path) -> Fraction:
    match = SQRT2_PATTERN.search(Path(path).read_text())
    if match is None:
        raise ValueError("PyWavelets coefficient header does not define sqrt2")
    return decimal_fraction(match["value"].strip())


def round_to_binary64(values: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(Fraction.from_float(float(value)) for value in values)


def orthogonal_filter_bank(coefficients: tuple[Fraction, ...]) -> FilterBank:
    length = len(coefficients)
    dec_lo = tuple(reversed(coefficients))
    rec_hi = tuple(
        (-1 if index % 2 else 1) * coefficients[length - 1 - index]
        for index in range(length)
    )
    dec_hi = tuple(
        (-1 if (length - 1 - index) % 2 else 1) * coefficients[index]
        for index in range(length)
    )
    return FilterBank(dec_lo, dec_hi, coefficients, rec_hi, True)


def biorthogonal_filter_bank(
    arrays: dict[str, tuple[Fraction, ...]], order_n: int, order_m: int
) -> FilterBank:
    maximum_m = {1: 5, 2: 8, 3: 9, 4: 4, 5: 5, 6: 8}[order_n]
    length = 2 * order_m if order_n == 1 else 2 * order_m + 2
    offset = maximum_m - order_m
    reconstruction = arrays[f"bior{order_n}_0"]
    decomposition = arrays[f"bior{order_n}_{order_m}"]

    rec_lo = tuple(reconstruction[index + offset] for index in range(length))
    dec_lo = tuple(decomposition[length - 1 - index] for index in range(length))
    rec_hi = tuple(
        (-1 if index % 2 else 1) * decomposition[length - 1 - index]
        for index in range(length)
    )
    dec_hi = tuple(
        (-1 if (length - 1 - index) % 2 else 1) * reconstruction[index + offset]
        for index in range(length)
    )
    return FilterBank(dec_lo, dec_hi, rec_lo, rec_hi, False)


def load_filter_bank(
    coefficients_header: Path,
    wavelet: str,
    decimal_digits: int = 80,
    binary64: bool = False,
) -> FilterBank:
    coefficients_header = Path(coefficients_header)
    arrays = parse_coefficient_arrays(coefficients_header)
    if binary64:
        arrays = {name: round_to_binary64(values) for name, values in arrays.items()}
    name = "db1" if wavelet == "haar" else wavelet

    match = re.fullmatch(r"(db|sym|coif)(\d+)", name)
    if match is not None:
        family, order = match.groups()
        coefficients = arrays[f"{family}{order}"]
        if family == "coif":
            if binary64:
                scale = float(parse_sqrt2(coefficients_header))
                coefficients = tuple(
                    Fraction.from_float(float(value) * scale) for value in coefficients
                )
            else:
                with localcontext() as context:
                    context.prec = decimal_digits
                    scale = decimal_fraction(Decimal(2).sqrt())
                coefficients = tuple(value * scale for value in coefficients)
        return orthogonal_filter_bank(coefficients)

    match = re.fullmatch(r"(bior|rbio)(\d+)\.(\d+)", name)
    if match is not None:
        family, order_n, order_m = match.groups()
        bank = biorthogonal_filter_bank(arrays, int(order_n), int(order_m))
        if family == "bior":
            return bank
        return FilterBank(
            tuple(reversed(bank.rec_lo)),
            tuple(reversed(bank.rec_hi)),
            tuple(reversed(bank.dec_lo)),
            tuple(reversed(bank.dec_hi)),
            False,
        )

    if name == "dmey":
        return orthogonal_filter_bank(arrays[name])

    raise ValueError(f"Unsupported discrete wavelet: {wavelet}")
