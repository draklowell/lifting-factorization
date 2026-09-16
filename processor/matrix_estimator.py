from lifting.matrix_estimation import (
    EstimationSpecification,
    build_solution_space,
    build_specification,
    build_system,
    normalize,
    reconstruct,
    select_target_monomial,
    solve,
)


class MatrixEstimationStage:
    NONE = 0
    SYSTEM_BUILT = 1
    TARGET_SELECTED = 2
    SOLUTION_SPACE_BUILT = 3
    SOLVED = 4
    NORMALIZED = 5


class MatrixEstimator:
    stage: MatrixEstimationStage
    esp: float = None
    H = None
    spec: EstimationSpecification = None
    T = None
    x_initial = None
    y = None
    x_particular = None
    x_basis = None
    x = None

    def __init__(self, H, eps=None):
        self.stage = MatrixEstimationStage.NONE
        self.esp = eps
        self.H = H
        self.spec = build_specification(H)

    def build_system(self):
        if self.stage != MatrixEstimationStage.NONE:
            raise ValueError("Invalid stage for system building step")

        self.T, self.x_initial = build_system(self.H, self.spec)
        self.stage = MatrixEstimationStage.SYSTEM_BUILT

        return self.T, self.x_initial

    def select_target_monomial(self, return_eps=False, normalize_to=None):
        if self.stage != MatrixEstimationStage.SYSTEM_BUILT:
            raise ValueError("Invalid stage for target selection step")

        self.y, eps = select_target_monomial(
            self.H,
            self.spec,
            eps=self.esp,
            return_eps=True,
            normalize_to=normalize_to,
        )
        self.stage = MatrixEstimationStage.TARGET_SELECTED

        if return_eps:
            return self.y, eps

        return self.y

    def build_solution_space(self):
        if self.stage != MatrixEstimationStage.TARGET_SELECTED:
            raise ValueError("Invalid stage for solution space building step")

        self.x_particular, self.x_basis = build_solution_space(self.T, self.y)
        self.stage = MatrixEstimationStage.SOLUTION_SPACE_BUILT

        return self.x_particular, self.x_basis

    def solve(self):
        if self.stage != MatrixEstimationStage.SOLUTION_SPACE_BUILT:
            raise ValueError("Invalid stage for solving step")

        self.x = solve(self.x_initial, self.x_particular, self.x_basis)
        self.stage = MatrixEstimationStage.SOLVED

        return self.x

    def normalize(self, normalize_to=1):
        if self.stage != MatrixEstimationStage.SOLVED:
            raise ValueError("Invalid stage for normalization step")

        self.x = normalize(self.x, self.y, normalize_to)
        self.stage = MatrixEstimationStage.NORMALIZED

        return self.x

    def reconstruct(self):
        if self.stage not in {
            MatrixEstimationStage.SOLVED,
            MatrixEstimationStage.NORMALIZED,
        }:
            raise ValueError("Invalid stage for reconstruction step")

        return reconstruct(self.H, self.x, self.spec)
