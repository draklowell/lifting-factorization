import json
import math
import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
import pywt
from sage.all import QQ, LaurentPolynomialRing

from lifting.projection import ProjectionStrategy
from processor.candidate_search import candidate_score
from processor.processor import Processor
from processor.utils import build_polyphase_matrix, matrix_ss, serialize_lifting_scheme
from reproducibility.metrics import (
    analysis_polyphase,
    matrix_l2_error,
    paraunitary_residual,
    sage_polyphase_to_numeric,
    scheme_polyphase,
)
from reproducibility.pywavelets import CandidateSelector, score_scheme
from reproducibility.pywavelets_source import load_filter_bank

SUPPORTED_WAVELETS = tuple(pywt.wavelist(kind="discrete"))


@dataclass(frozen=True)
class GenerationResult:
    scheme: dict[str, Any]
    metadata: dict[str, Any]


def baseline_result(
    baseline: dict[str, Any],
    started: float,
    reason: str,
    **metadata,
) -> GenerationResult:
    return GenerationResult(
        scheme=baseline,
        metadata={
            "strategy_selected": "baseline-retained",
            "support_extension": 0,
            "generation_runtime_seconds": time.perf_counter() - started,
            "candidate_count": 1,
            "selection_reason": reason,
            **metadata,
        },
    )


def generate_improved_scheme(
    wavelet_name: str,
    baseline_scheme_path: Path,
    coefficients_header: Path,
    beam_width: int = 16,
    projection_strategy: ProjectionStrategy = ProjectionStrategy.PRESERVE_LOW_PASS,
    support_margin: int = 0,
) -> GenerationResult:
    if support_margin != 0:
        raise ValueError("Support extension is not implemented by this projection")

    started = time.perf_counter()
    baseline = json.loads(Path(baseline_scheme_path).read_text())
    wavelet = pywt.Wavelet(wavelet_name)

    if not wavelet.orthogonal or wavelet_name == "dmey":
        reason = (
            "biorthogonal baseline is already near the FP32 floor"
            if not wavelet.orthogonal
            else "fixed-support exact projection is not closer to dmey"
        )
        return baseline_result(baseline, started, reason)

    source_bank = load_filter_bank(coefficients_header, wavelet_name)
    ring = LaurentPolynomialRing(QQ, names=("z",))
    source_matrix = build_polyphase_matrix(
        source_bank.dec_lo,
        source_bank.dec_hi,
        ring,
    )
    reference_matrix = build_polyphase_matrix(
        wavelet.dec_lo,
        wavelet.dec_hi,
        ring,
    )

    projection_started = time.perf_counter()
    processor = Processor(StringIO())
    projected_matrix = processor.estimate(
        source_matrix,
        strategy=projection_strategy,
    )
    projection_runtime = time.perf_counter() - projection_started

    search_started = time.perf_counter()
    selector = CandidateSelector(wavelet, wavelet.dec_len)
    selected_steps, selected_delays = processor.factorize(
        projected_matrix,
        candidate_selector=selector,
        beam_width=beam_width,
    )
    search_runtime = time.perf_counter() - search_started
    if selector.selected is None:
        raise RuntimeError("Factorization selector did not retain its result")

    scheme = serialize_lifting_scheme(
        selected_steps,
        selected_delays,
        wavelet.dec_len,
        coefficients_as_float=True,
    )
    reference_numeric = analysis_polyphase(wavelet.dec_lo, wavelet.dec_hi)
    selected_numeric = scheme_polyphase(scheme)
    original_invariant = paraunitary_residual(reference_numeric)
    selected_invariant = paraunitary_residual(selected_numeric)
    invariant_limit = max(100.0 * original_invariant, 1e-12)
    baseline_score = candidate_score(
        baseline,
        lambda value: score_scheme(value, wavelet),
    )

    if (
        selected_invariant > invariant_limit
        or selector.selected.score >= baseline_score
    ):
        reason = (
            "candidate invariant exceeded the source tolerance"
            if selected_invariant > invariant_limit
            else "candidate did not improve actual FP32 probe error"
        )
        return baseline_result(
            baseline,
            started,
            reason,
            projection_runtime_seconds=projection_runtime,
            factorization_search_runtime_seconds=search_runtime,
            factorization_runtime_seconds=projection_runtime + search_runtime,
            candidate_count=selector.candidate_count,
        )

    projected_numeric = sage_polyphase_to_numeric(projected_matrix)
    matrix_estimation_error = math.sqrt(
        float(matrix_ss(reference_matrix, projected_matrix))
    )
    fp64_total_error = matrix_l2_error(reference_numeric, selected_numeric)
    scheme["meta"] = {
        "l2": {
            "qq": matrix_estimation_error,
            "fp64": fp64_total_error,
        }
    }
    return GenerationResult(
        scheme=scheme,
        metadata={
            "strategy_selected": (
                f"{projection_strategy.value}+margin-{support_margin}+beam-{beam_width}"
            ),
            "support_extension": support_margin,
            "generation_runtime_seconds": time.perf_counter() - started,
            "projection_runtime_seconds": projection_runtime,
            "factorization_search_runtime_seconds": search_runtime,
            "factorization_runtime_seconds": projection_runtime + search_runtime,
            "candidate_count": selector.candidate_count,
            "selection_reason": (
                "lower actual FP32 probe error with invariant preserved"
            ),
            "source_to_estimated_l2": matrix_l2_error(
                analysis_polyphase(
                    source_bank.dec_lo,
                    source_bank.dec_hi,
                ),
                projected_numeric,
            ),
            "matrix_estimation_error": matrix_estimation_error,
            "original_invariant_residual": original_invariant,
            "post_estimation_invariant_residual": paraunitary_residual(
                projected_numeric
            ),
            "serialized_fp64_invariant_residual": selected_invariant,
            "estimated_to_serialized_fp64_error": matrix_l2_error(
                projected_numeric,
                selected_numeric,
            ),
            "estimated_to_serialized_fp32_error": matrix_l2_error(
                projected_numeric,
                scheme_polyphase(scheme, np.float32),
            ),
            "exact_factorization_residual": 0.0,
            "selected_probe_score": list(selector.selected.score),
            "baseline_probe_score": list(baseline_score),
        },
    )
