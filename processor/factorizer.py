from lifting.factorization import (
    canonicalize_scale_delays,
    euclidean,
    factorization_candidates,
    normalize_det,
    pack,
    reconstruct,
    recover,
)
from processor.utils import reconstruct as reconstruct_steps


class FactorizerStage:
    NONE = 0
    NORMALIZED = 1
    EUCLIDEAN = 2
    RECONSTRUCT = 3
    RECOVER = 4


class Factorizer:
    stage: FactorizerStage
    P = None
    delays: tuple[int, int] | None = None
    qs = None
    a = None
    P0 = None
    s = None
    b = None

    def __init__(self, P):
        self.stage = FactorizerStage.NONE
        self.original = P
        self.P = P

    def normalize(self):
        if self.stage != FactorizerStage.NONE:
            raise ValueError("Invalid stage for normalization step")

        self.P, self.delays = normalize_det(self.P)
        self.stage = FactorizerStage.NORMALIZED

        return self.P, self.delays

    def euclidean(self):
        if self.stage not in {FactorizerStage.NONE, FactorizerStage.NORMALIZED}:
            raise ValueError("Invalid stage for euclidean step")

        self.qs, self.a = euclidean(self.P[0, 0], self.P[1, 0])
        self.stage = FactorizerStage.EUCLIDEAN

        return self.qs, self.a

    def build_candidates(self, beam_width=16):
        if self.stage != FactorizerStage.NORMALIZED:
            raise ValueError("Invalid stage for candidate factorization step")

        candidates = []
        for steps in factorization_candidates(self.P, beam_width=beam_width):
            runtime_steps, runtime_delays = canonicalize_scale_delays(
                steps,
                self.delays,
            )
            assert (
                reconstruct_steps(
                    runtime_steps,
                    runtime_delays,
                    self.original.base_ring(),
                )
                == self.original
            )
            candidates.append((runtime_steps, runtime_delays))
        return candidates

    def reconstruct(self):
        if self.stage != FactorizerStage.EUCLIDEAN:
            raise ValueError("Invalid stage for reconstruct step")

        self.P0 = reconstruct(self.qs, self.a)
        self.stage = FactorizerStage.RECONSTRUCT

        return self.P0

    def recover(self):
        if self.stage != FactorizerStage.RECONSTRUCT:
            raise ValueError("Invalid stage for recover step")

        self.s, self.b = recover(self.P, self.P0)
        self.stage = FactorizerStage.RECOVER

        return self.s, self.b

    def pack(self):
        if self.stage != FactorizerStage.RECOVER:
            raise ValueError("Invalid stage for pack step")

        return (
            pack(self.qs, self.a, self.s, self.b, self.P),
            self.delays,
        )
