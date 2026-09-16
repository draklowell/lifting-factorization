from dataclasses import dataclass, field

import numpy as np
import pywt

from processor.candidate_search import ScoredFactorization, select_factorization
from usage import dtypes
from usage.lifting import LiftingScheme


def signal_cases(length: int, seed: int) -> tuple[np.ndarray, ...]:
    rng = np.random.default_rng(seed + length)
    index = np.arange(length, dtype=np.float64)
    impulse = np.zeros(length, dtype=np.float64)
    impulse[length // 2] = 1.0
    boundary_impulse = np.zeros(length, dtype=np.float64)
    boundary_impulse[0] = 1.0
    alternating = np.where(index % 2 == 0, 1.0, -1.0)
    cancellation = alternating * (1.0 + index * np.finfo(np.float32).eps)
    return (
        rng.normal(size=length),
        rng.uniform(-1.0, 1.0, size=length),
        impulse,
        boundary_impulse,
        alternating,
        np.ones(length, dtype=np.float64),
        np.linspace(-1.0, 1.0, length, dtype=np.float64),
        cancellation,
    )


def relative_l2(reference: np.ndarray, actual: np.ndarray) -> float:
    reference = np.asarray(reference, dtype=np.float64)
    actual = np.asarray(actual, dtype=np.float64)
    if reference.shape != actual.shape:
        raise ValueError(
            f"Shape mismatch: reference {reference.shape}, actual {actual.shape}"
        )
    denominator = max(float(np.linalg.norm(reference)), np.finfo(np.float64).tiny)
    return float(np.linalg.norm(actual - reference)) / denominator


def paired(approximation: np.ndarray, detail: np.ndarray) -> np.ndarray:
    return np.concatenate(
        [
            np.asarray(approximation, dtype=np.float64),
            np.asarray(detail, dtype=np.float64),
        ]
    )


def score_scheme(
    scheme: dict,
    wavelet: pywt.Wavelet,
    lengths: tuple[int, ...] = (32, 101, 1009),
    seed: int = 1729,
) -> tuple[float, ...]:
    lifting = LiftingScheme.from_object(
        scheme,
        mode="symmetric",
        dtype=dtypes.float32,
    )
    forward_error = 0.0
    inverse_error = 0.0
    round_trip_error = 0.0
    peak_amplification = 0.0
    finite = True

    for length in lengths:
        for signal in signal_cases(length, seed):
            approximation_ref, detail_ref = pywt.dwt(
                signal,
                wavelet,
                mode="symmetric",
            )
            inverse_ref = pywt.idwt(
                approximation_ref,
                detail_ref,
                wavelet,
                mode="symmetric",
            )
            approximation, detail, forward_peak = lifting.forward_with_amplification(
                signal.astype(np.float32)
            )
            inverse, _ = lifting.inverse_with_amplification(
                approximation_ref.astype(np.float32),
                detail_ref.astype(np.float32),
            )
            round_trip, inverse_peak = lifting.inverse_with_amplification(
                approximation,
                detail,
            )

            finite = finite and all(
                np.isfinite(array).all()
                for array in (approximation, detail, inverse, round_trip)
            )
            forward_error = max(
                forward_error,
                relative_l2(
                    paired(approximation_ref, detail_ref),
                    paired(approximation, detail),
                ),
            )
            inverse_error = max(inverse_error, relative_l2(inverse_ref, inverse))
            round_trip_error = max(
                round_trip_error,
                relative_l2(signal, round_trip[:length]),
            )
            peak_amplification = max(
                peak_amplification,
                forward_peak,
                inverse_peak,
            )

    return (
        0.0 if finite else 1.0,
        max(forward_error, inverse_error),
        round_trip_error,
        peak_amplification,
    )


@dataclass
class CandidateSelector:
    wavelet: pywt.Wavelet
    tap_size: int
    lengths: tuple[int, ...] = (32, 101, 1009)
    seed: int = 1729
    selected: ScoredFactorization | None = field(default=None, init=False)
    candidate_count: int = field(default=0, init=False)

    def __call__(self, candidates):
        self.candidate_count = len(candidates)
        self.selected = select_factorization(
            candidates,
            self.tap_size,
            lambda scheme: score_scheme(
                scheme,
                self.wavelet,
                lengths=self.lengths,
                seed=self.seed,
            ),
        )
        return self.selected.steps, self.selected.delays
