"""Reference geometry for the 3x3 segmented electrostatic phase plate."""

from __future__ import annotations

import numpy as np

from iqs.actuators import ElectrostaticDomain, segmented_three_plate_aperture_array


UM = 1e-6
ARRAY_SHAPE = (3, 3)
PITCH_M = 8 * UM
APERTURE_RADIUS_M = 2.5 * UM
DOMAIN_EXTENT_M = (
    (-16 * UM, 16 * UM),
    (-16 * UM, 16 * UM),
    (-15 * UM, 15 * UM),
)
DOMAIN_SHAPE = (65, 65, 97)


def aperture_centers():
    """Return the phase-readout centers in T31's established ordering."""
    x = (np.arange(ARRAY_SHAPE[0]) - 1) * PITCH_M
    y = (np.arange(ARRAY_SHAPE[1]) - 1) * PITCH_M
    return tuple((float(xi), float(yi)) for xi in x for yi in y)


def segmented_phase_plate_model(voltages):
    """Build the reference three-plate, nine-segment electrostatic model."""
    controls = np.asarray(voltages, dtype=float)
    if controls.size != int(np.prod(ARRAY_SHAPE)):
        raise ValueError("segmented phase plate requires nine control voltages")
    domain = ElectrostaticDomain(
        extent=DOMAIN_EXTENT_M,
        shape=DOMAIN_SHAPE,
        boundary_policy="grounded_box",
    )
    return segmented_three_plate_aperture_array(
        domain,
        plate_z_m=(-6.25 * UM, 0.0, 6.25 * UM),
        plate_thickness_m=1.25 * UM,
        center_voltages_V=controls.reshape(ARRAY_SHAPE),
        aperture_radius_m=APERTURE_RADIUS_M,
        pitch_m=PITCH_M,
        array_shape=ARRAY_SHAPE,
        segment_gap_m=1.0 * UM,
        outer_plate_half_width_m=(14 * UM, 14 * UM),
        metadata={"device": "reference 3x3 segmented phase plate"},
    )


def aperture_averaged_profiles(field_map):
    """Average the potential over each circular aperture at every z plane."""
    xx, yy = np.meshgrid(field_map.x_m, field_map.y_m, indexing="xy")
    profiles = []
    for cx, cy in aperture_centers():
        mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= APERTURE_RADIUS_M ** 2
        profiles.append(np.asarray(field_map.potential_V[:, mask].mean(axis=1)))
    return np.asarray(profiles)


__all__ = [
    "APERTURE_RADIUS_M",
    "ARRAY_SHAPE",
    "DOMAIN_EXTENT_M",
    "DOMAIN_SHAPE",
    "PITCH_M",
    "UM",
    "aperture_averaged_profiles",
    "aperture_centers",
    "segmented_phase_plate_model",
]
