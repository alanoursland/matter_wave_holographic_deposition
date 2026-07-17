"""Programmable device models."""

from .squid_array import SQUIDArray
from .segmented_phase_plate import (
    APERTURE_RADIUS_M,
    ARRAY_SHAPE,
    DOMAIN_EXTENT_M,
    DOMAIN_SHAPE,
    PITCH_M,
    aperture_averaged_profiles,
    aperture_centers,
    segmented_phase_plate_model,
)

__all__ = [
    "APERTURE_RADIUS_M",
    "ARRAY_SHAPE",
    "DOMAIN_EXTENT_M",
    "DOMAIN_SHAPE",
    "PITCH_M",
    "SQUIDArray",
    "aperture_averaged_profiles",
    "aperture_centers",
    "segmented_phase_plate_model",
]
