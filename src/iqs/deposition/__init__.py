"""Layer-by-layer stochastic surface deposition models."""

from .electrical import (
    ContactArrayElectricalResult,
    ContactElectricalStack,
    extract_contact_array_electrical,
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
    "deposit_layer",
    "extract_contact_array_electrical",
]
