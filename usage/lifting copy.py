import json
from typing import Any

import numpy as np

import dtypes


class ShiftedArray:
    shift: int
    data: np.ndarray

    def __init__(self, data: np.ndarray, shift: int):
        self.data = np.array(data)
        self.shift = shift

    @property
    def end(self) -> int:
        return self.shift + len(self.data)

    @classmethod
    def from_object(cls, obj: dict[str, Any]) -> "ShiftedArray":
        data = []
        for coeff in obj["coefficients"]:
            if isinstance(coeff, (int, float)):
                data.append(coeff)
            else:
                data.append(coeff["numerator"] / coeff["denominator"])

        return cls(data, obj["shift"])

    def __repr__(self):
        return f"ShiftedArray(data={self.data}, shift={self.shift})"

    def __neg__(self) -> "ShiftedArray":
        return ShiftedArray(-self.data, self.shift)

    def __sub__(self, other: "ShiftedArray") -> "ShiftedArray":
        return self + (-other)


class DType:
    def __init__(self, mul_dtype: dtypes.dtype, add_dtype: dtypes.dtype):
        self.mul_dtype = mul_dtype
        self.add_dtype = add_dtype

    def mul(self, a: ShiftedArray, b: int | float | ShiftedArray) -> ShiftedArray:
        if isinstance(b, (int, float)):
            return ShiftedArray(self.mul_dtype.mul(a.data, b), a.shift)

        if not isinstance(b, ShiftedArray):
            raise ValueError("Unsupported type for multiplication")

        shift_new = a.shift + b.shift + min(len(a.data), len(b.data)) - 1
        data_new = self.mul_dtype.conv(a.data, b.data)
        return ShiftedArray(data_new, shift_new)

    def add(self, a: ShiftedArray, b: ShiftedArray) -> ShiftedArray:
        if not isinstance(b, ShiftedArray):
            raise ValueError("Unsupported type for addition")

        shift_new = max(a.shift, b.shift)
        end_new = min(a.end, b.end)

        data1 = a.data[shift_new - a.shift : end_new - a.shift]
        data2 = b.data[shift_new - b.shift : end_new - b.shift]
        data_new = self.add_dtype.add(data1, data2)

        return ShiftedArray(data_new, shift_new)


class LiftingStep:
    def forward(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return even, odd

    def inverse(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return even, odd


class LiftingStepPredict(LiftingStep):
    kernel: ShiftedArray
    dtype: DType

    def __init__(self, kernel: ShiftedArray, dtype: DType):
        self.kernel = kernel
        self.dtype = dtype

    def forward(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        pred = self.dtype.mul(even, self.kernel)
        return even, self.dtype.add(odd, pred)

    def inverse(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        pred = self.dtype.mul(even, self.kernel)
        return even, self.dtype.add(odd, -pred)


class LiftingStepUpdate(LiftingStep):
    kernel: ShiftedArray
    dtype: DType

    def __init__(self, kernel: ShiftedArray, dtype: DType):
        self.kernel = kernel
        self.dtype = dtype

    def forward(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        upd = self.dtype.mul(self.kernel, odd)
        return self.dtype.add(even, upd), odd

    def inverse(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        upd = self.dtype.mul(self.kernel, odd)
        return self.dtype.add(even, -upd), odd


class LiftingStepScaleEven(LiftingStep):
    factor: float
    dtype: DType

    def __init__(self, factor: float, dtype: DType):
        self.factor = factor
        self.dtype = dtype

    def forward(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return self.dtype.mul(even, self.factor), odd

    def inverse(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return self.dtype.mul(even, 1 / self.factor), odd


class LiftingStepScaleOdd(LiftingStep):
    factor: float
    dtype: DType

    def __init__(self, factor: float, dtype: DType):
        self.factor = factor
        self.dtype = dtype

    def forward(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return even, self.dtype.mul(odd, self.factor)

    def inverse(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return even, self.dtype.mul(odd, 1 / self.factor)


class LiftingStepSwap(LiftingStep):
    def forward(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return odd, even

    def inverse(
        self, even: ShiftedArray, odd: ShiftedArray
    ) -> tuple[ShiftedArray, ShiftedArray]:
        return odd, even


class LiftingScheme:
    def __init__(
        self,
        mode: str,
        tap_size: int,
        steps: list[LiftingStep],
        delays: tuple[int, int],
    ):
        self.mode = mode
        self.tap_size = tap_size
        self.steps = steps
        self.delays = delays

    @classmethod
    def from_object(
        cls,
        obj: dict[str, Any],
        mode: str,
        predict_dtype: DType,
        update_dtype: DType,
        scale_dtype: DType,
    ):
        steps = []

        for step in obj["steps"]:
            kernel = ShiftedArray.from_object(step)

            if step["type"] == "predict":
                steps.append(LiftingStepPredict(kernel, predict_dtype))
            elif step["type"] == "update":
                steps.append(LiftingStepUpdate(kernel, update_dtype))
            elif step["type"] == "scale-even":
                assert kernel.shift == 0, "Scale steps must have zero shift"
                assert (
                    len(kernel.data) == 1
                ), "Scale steps must have a single coefficient"
                steps.append(LiftingStepScaleEven(kernel.data[0], scale_dtype))
            elif step["type"] == "scale-odd":
                assert kernel.shift == 0, "Scale steps must have zero shift"
                assert (
                    len(kernel.data) == 1
                ), "Scale steps must have a single coefficient"
                steps.append(LiftingStepScaleOdd(kernel.data[0], scale_dtype))
            elif step["type"] == "swap":
                steps.append(LiftingStepSwap())
            else:
                raise ValueError(f"Unsupported lifting step type: {step['type']}")

        return cls(
            mode=mode,
            tap_size=obj["tap_size"],
            steps=steps,
            delays=(
                obj["delay"]["even"],
                obj["delay"]["odd"],
            ),
        )

    @classmethod
    def from_file(
        cls,
        path: str,
        mode: str,
        predict_dtype: DType,
        update_dtype: DType,
        scale_dtype: DType,
    ):
        with open(path, "r") as f:
            obj = json.load(f)

        return cls.from_object(
            obj,
            mode=mode,
            predict_dtype=predict_dtype,
            update_dtype=update_dtype,
            scale_dtype=scale_dtype,
        )

    def pad(self, data: np.ndarray, shape: int) -> np.ndarray:
        if self.mode == "symmetric":
            return np.pad(data, shape, mode="symmetric")

        if self.mode == "periodic":
            return np.pad(data, shape, mode="wrap")

        if self.mode == "zero":
            return np.pad(data, shape, mode="constant", constant_values=0)

        if self.mode == "constant":
            return np.pad(data, shape, mode="edge")

        raise ValueError(f"Unsupported padding mode: {self.mode}")

    def forward(self, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        L = self.tap_size
        length = (len(data) + L - 1) // 2
        direct_shift = L // 2

        data = self.pad(data, (L - 1, L - 1))

        even = data[::2]
        odd = data[1::2]

        # Shift when multiplying by kernel without
        # polyphase split
        even = ShiftedArray(even, self.delays[0])
        odd = ShiftedArray(odd, self.delays[1])

        for step in self.steps:
            even, odd = step.forward(even, odd)

        even_shift = even.shift - direct_shift
        odd_shift = odd.shift - direct_shift

        if even_shift > 0 or odd_shift > 0:
            raise ValueError("Coefficients have positive shift after forward lifting")

        even = even.data[-even_shift : -even_shift + length]
        odd = odd.data[-odd_shift : -odd_shift + length]

        return even, odd

    def inverse(self, even: np.ndarray, odd: np.ndarray) -> np.ndarray:
        L = self.tap_size
        length = 2 * len(even) - L + 2
        direct_shift = L // 2

        even = ShiftedArray(even, 0)
        odd = ShiftedArray(odd, 0)

        for step in reversed(self.steps):
            even, odd = step.inverse(even, odd)

        even_shift = even.shift + self.delays[1] - direct_shift
        odd_shift = odd.shift + self.delays[0] - direct_shift

        # Odd is shifted
        odd_shift += 1

        even = even.data[-even_shift : -even_shift + length // 2]
        odd = odd.data[-odd_shift : -odd_shift + length // 2]

        data = np.empty(length, dtype=np.float64)
        data[1::2] = even
        data[::2] = odd

        return data
