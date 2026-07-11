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
]
