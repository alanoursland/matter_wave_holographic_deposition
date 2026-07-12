"""Layer-by-layer stochastic surface deposition models."""

from .surface import (
    DepositionMaterial,
    DepositionResult,
    SurfaceState,
    deposit_layer,
)

__all__ = [
    "DepositionMaterial",
    "DepositionResult",
    "SurfaceState",
    "deposit_layer",
]
