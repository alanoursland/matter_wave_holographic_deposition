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
]
