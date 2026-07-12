"""Layer-by-layer stochastic surface deposition models."""

from .electrical import (
    ContactArrayElectricalResult,
    ContactElectricalStack,
    SOIIsolationResult,
    SOIIsolationStack,
    extract_contact_array_electrical,
    extract_soi_isolation,
)
from .surface import (
    DepositionMaterial,
    DepositionResult,
    SurfaceState,
    deposit_layer,
)

__all__ = [
    "ContactArrayElectricalResult",
    "ContactElectricalStack",
    "DepositionMaterial",
    "DepositionResult",
    "SurfaceState",
    "SOIIsolationResult",
    "SOIIsolationStack",
    "deposit_layer",
    "extract_contact_array_electrical",
    "extract_soi_isolation",
]
