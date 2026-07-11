"""Matter-wave phase actuators and physical authority checks."""

from .phase import (
    AchromaticPhaseResponse,
    CoplanarSquareLoopArray,
    ElectrostaticPlateGeometry,
    ElectrostaticPhaseResponse,
    ElectrostaticValidityReport,
    IdealPhasePlate,
    PhaseAuthorityError,
    PhaseAuthorityReport,
    resolve_phase_response,
)
from .electrostatic_3d import (
    analyze_electrostatic_plate_3d,
    Electrostatic3DAnalysis,
    ElectrostaticFieldMap,
    ElectrostaticMultislicePropagator,
    LaplaceFringeFieldModel,
    MultisliceComparison,
)
from .fem_import import (
    FEMFieldImportReport,
    FEMFieldImportResult,
    FEMGridSpec,
    load_fem_csv,
    load_fem_npz,
    save_fem_npz,
)

__all__ = [
    "AchromaticPhaseResponse",
    "CoplanarSquareLoopArray",
    "ElectrostaticPlateGeometry",
    "ElectrostaticPhaseResponse",
    "ElectrostaticValidityReport",
    "IdealPhasePlate",
    "PhaseAuthorityError",
    "PhaseAuthorityReport",
    "resolve_phase_response",
    "ElectrostaticFieldMap",
    "Electrostatic3DAnalysis",
    "ElectrostaticMultislicePropagator",
    "LaplaceFringeFieldModel",
    "MultisliceComparison",
    "analyze_electrostatic_plate_3d",
    "FEMFieldImportReport",
    "FEMFieldImportResult",
    "FEMGridSpec",
    "load_fem_csv",
    "load_fem_npz",
    "save_fem_npz",
]
