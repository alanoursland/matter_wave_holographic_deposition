"""Electrical extraction for patterned metal/semiconductor contacts."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ContactElectricalStack:
    """Measured stack parameters used by the geometry extractor.

    Resistivities use SI units. ``specific_contact_resistivity_ohm_m2`` is
    the interface quantity conventionally reported in ohm cm^2.
    """

    specific_contact_resistivity_ohm_m2: float
    semiconductor_resistivity_ohm_m: float
    test_bias_v: float = 1.0

    def __post_init__(self):
        values = (
            self.specific_contact_resistivity_ohm_m2,
            self.semiconductor_resistivity_ohm_m,
            self.test_bias_v,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("electrical stack values must be finite")
        if min(values) <= 0:
            raise ValueError("electrical stack values must be positive")


@dataclass(frozen=True)
class ContactArrayElectricalResult:
    effective_area_m2: np.ndarray
    interface_resistance_ohm: np.ndarray
    spreading_resistance_ohm: np.ndarray
    pair_resistance_ohm: np.ndarray
    pair_leakage_a: np.ndarray

    @property
    def worst_contact_resistance_ohm(self) -> float:
        total = self.interface_resistance_ohm + self.spreading_resistance_ohm
        return float(np.max(total))

    @property
    def worst_pair_leakage_a(self) -> float:
        finite = self.pair_leakage_a[np.isfinite(self.pair_leakage_a)]
        return float(np.max(finite)) if finite.size else 0.0


@dataclass(frozen=True)
class SOIIsolationStack:
    """Buried-oxide and exposed-surface isolation for etched SOI mesas."""

    box_resistivity_ohm_m: float
    box_thickness_m: float
    surface_sheet_resistance_ohm_sq: float
    test_bias_v: float = 1.0

    def __post_init__(self):
        values = (
            self.box_resistivity_ohm_m,
            self.box_thickness_m,
            self.surface_sheet_resistance_ohm_sq,
            self.test_bias_v,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("SOI isolation values must be finite")
        if min(values) <= 0:
            raise ValueError("SOI isolation values must be positive")


@dataclass(frozen=True)
class SOIIsolationResult:
    mesa_area_m2: np.ndarray
    pair_gap_m: np.ndarray
    box_pair_resistance_ohm: np.ndarray
    surface_pair_resistance_ohm: np.ndarray
    pair_leakage_a: np.ndarray

    @property
    def worst_pair_leakage_a(self) -> float:
        return float(np.max(self.pair_leakage_a))


@dataclass(frozen=True)
class TwoTerminalSOIStack:
    """Finite SOI channel and neighboring-mesa surface leakage parameters."""

    semiconductor_resistivity_ohm_m: float
    device_layer_thickness_m: float
    channel_length_m: float
    channel_width_m: float
    surface_sheet_resistance_ohm_sq: float
    monitor_geometry_factor: float
    test_bias_v: float = 1.0

    def __post_init__(self):
        values = (
            self.semiconductor_resistivity_ohm_m,
            self.device_layer_thickness_m,
            self.channel_length_m,
            self.channel_width_m,
            self.surface_sheet_resistance_ohm_sq,
            self.monitor_geometry_factor,
            self.test_bias_v,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("two-terminal SOI stack values must be finite")
        if min(values) <= 0:
            raise ValueError("two-terminal SOI stack values must be positive")


@dataclass(frozen=True)
class TwoTerminalSOIResult:
    contact_resistance_ohm: np.ndarray
    channel_resistance_ohm: float
    total_resistance_ohm: float
    device_current_a: float
    monitor_leakage_a: float
    device_to_leakage_ratio: float
    contact_voltage_drop_v: np.ndarray
    channel_voltage_drop_v: float
    contact_current_density_a_m2: np.ndarray
    channel_current_density_a_m2: float
    channel_power_w: float


def extract_two_terminal_soi_device(
    contact_resistance_ohm,
    effective_contact_area_m2,
    *,
    stack: TwoTerminalSOIStack,
) -> TwoTerminalSOIResult:
    """Extract a two-contact resistor and worst-case grounded-monitor leakage."""
    contact_r = np.asarray(contact_resistance_ohm, dtype=float)
    areas = np.asarray(effective_contact_area_m2, dtype=float)
    if contact_r.shape != (2,) or areas.shape != (2,):
        raise ValueError("exactly two contact resistances and areas are required")
    if (np.any(np.isnan(contact_r)) or np.any(contact_r <= 0)
            or np.any(~np.isfinite(areas)) or np.any(areas < 0)):
        raise ValueError("contact resistances and areas must be nonnegative and valid")

    channel_r = (
        stack.semiconductor_resistivity_ohm_m * stack.channel_length_m
        / (stack.channel_width_m * stack.device_layer_thickness_m)
    )
    total_r = float(np.sum(contact_r) + channel_r)
    current = float(stack.test_bias_v / total_r)
    monitor_leakage = float(
        stack.test_bias_v * stack.monitor_geometry_factor
        / stack.surface_sheet_resistance_ohm_sq
    )
    current_density = np.zeros(2, dtype=float)
    np.divide(current, areas, out=current_density, where=areas > 0)
    contact_drop = np.full(2, np.nan, dtype=float)
    np.multiply(current, contact_r, out=contact_drop, where=np.isfinite(contact_r))
    ratio = current / monitor_leakage if monitor_leakage > 0 else float("inf")
    return TwoTerminalSOIResult(
        contact_resistance_ohm=contact_r.copy(),
        channel_resistance_ohm=float(channel_r),
        total_resistance_ohm=total_r,
        device_current_a=current,
        monitor_leakage_a=monitor_leakage,
        device_to_leakage_ratio=float(ratio),
        contact_voltage_drop_v=contact_drop,
        channel_voltage_drop_v=float(current * channel_r),
        contact_current_density_a_m2=current_density,
        channel_current_density_a_m2=float(
            current / (stack.channel_width_m * stack.device_layer_thickness_m)
        ),
        channel_power_w=float(current ** 2 * channel_r),
    )


def extract_contact_array_electrical(
    metal_thickness_m,
    contact_masks,
    *,
    pixel_pitch_m: float,
    continuity_threshold_m: float,
    stack: ContactElectricalStack,
) -> ContactArrayElectricalResult:
    """Extract contact resistance and pair leakage from realized morphology.

    Each conducting pixel in a nominal contact contributes interface area.
    A circular equal-area contact gives the classical half-space spreading
    resistance. Pair leakage includes the first-order mutual spreading term
    for two surface contacts on one continuous semiconductor half-space.

    This is an ohmic DC extraction. It does not model Schottky barriers,
    quantum confinement, grain-boundary resistance, or dielectric isolation.
    """
    thickness = np.asarray(metal_thickness_m, dtype=float)
    masks = np.asarray(contact_masks, dtype=bool)
    if thickness.ndim != 2 or masks.ndim != 3 or masks.shape[1:] != thickness.shape:
        raise ValueError("contact masks must have shape (contacts, *metal.shape)")
    if masks.shape[0] < 2:
        raise ValueError("at least two contacts are required")
    if not np.all(np.isfinite(thickness)) or np.any(thickness < 0):
        raise ValueError("metal thickness must be finite and nonnegative")
    if not np.isfinite(pixel_pitch_m) or pixel_pitch_m <= 0:
        raise ValueError("pixel pitch must be positive")
    if not np.isfinite(continuity_threshold_m) or continuity_threshold_m <= 0:
        raise ValueError("continuity threshold must be positive")

    conducting = thickness >= continuity_threshold_m
    pixel_area = pixel_pitch_m ** 2
    areas = np.array([
        np.count_nonzero(conducting & mask) * pixel_area for mask in masks
    ])
    interface_r = np.full(masks.shape[0], np.inf)
    spreading_r = np.full(masks.shape[0], np.inf)
    present = areas > 0
    interface_r[present] = stack.specific_contact_resistivity_ohm_m2 / areas[present]
    radii = np.sqrt(areas[present] / np.pi)
    spreading_r[present] = stack.semiconductor_resistivity_ohm_m / (4.0 * radii)

    coordinates = np.indices(thickness.shape, dtype=float)
    centers = np.empty((masks.shape[0], 2), dtype=float)
    for index, mask in enumerate(masks):
        if not np.any(mask):
            raise ValueError("contact masks must not be empty")
        centers[index] = [np.mean(axis[mask]) for axis in coordinates]
    centers *= pixel_pitch_m

    count = masks.shape[0]
    pair_r = np.full((count, count), np.inf)
    pair_i = np.zeros((count, count))
    for left in range(count):
        for right in range(left + 1, count):
            if not (present[left] and present[right]):
                continue
            distance = float(np.linalg.norm(centers[left] - centers[right]))
            mutual = stack.semiconductor_resistivity_ohm_m / (np.pi * distance)
            substrate = max(
                spreading_r[left] + spreading_r[right] - mutual,
                0.0,
            )
            resistance = interface_r[left] + interface_r[right] + substrate
            leakage = stack.test_bias_v / resistance
            pair_r[left, right] = pair_r[right, left] = resistance
            pair_i[left, right] = pair_i[right, left] = leakage

    return ContactArrayElectricalResult(
        effective_area_m2=areas,
        interface_resistance_ohm=interface_r,
        spreading_resistance_ohm=spreading_r,
        pair_resistance_ohm=pair_r,
        pair_leakage_a=pair_i,
    )


def extract_soi_isolation(
    mesa_masks,
    *,
    pixel_pitch_m: float,
    stack: SOIIsolationStack,
) -> SOIIsolationResult:
    """Extract BOX and exposed-surface leakage between isolated SOI mesas.

    The BOX branch sends current vertically through each mesa footprint to an
    ideally conducting handle wafer, which is a conservative upper-current
    bound because handle spreading resistance is omitted. The surface branch
    is a sheet-resistance path across the edge-to-edge trench gap.
    """
    masks = np.asarray(mesa_masks, dtype=bool)
    if masks.ndim != 3 or masks.shape[0] < 2:
        raise ValueError("mesa masks must have shape (mesas, rows, columns)")
    if not np.isfinite(pixel_pitch_m) or pixel_pitch_m <= 0:
        raise ValueError("pixel pitch must be positive")
    if np.any(np.sum(masks, axis=(1, 2)) == 0):
        raise ValueError("mesa masks must not be empty")

    pixel_area = pixel_pitch_m ** 2
    areas = np.sum(masks, axis=(1, 2), dtype=float) * pixel_area
    widths = np.sqrt(areas)
    coordinates = np.indices(masks.shape[1:], dtype=float)
    centers = np.array([
        [np.mean(axis[mask]) for axis in coordinates] for mask in masks
    ]) * pixel_pitch_m

    count = masks.shape[0]
    gaps = np.full((count, count), np.inf)
    box_r = np.full((count, count), np.inf)
    surface_r = np.full((count, count), np.inf)
    leakage = np.zeros((count, count))
    vertical_r = stack.box_resistivity_ohm_m * stack.box_thickness_m / areas
    for left in range(count):
        for right in range(left + 1, count):
            distance = float(np.linalg.norm(centers[left] - centers[right]))
            gap = distance - 0.5 * (widths[left] + widths[right])
            if gap <= 0:
                raise ValueError("mesa masks overlap or touch")
            facing_width = min(widths[left], widths[right])
            pair_box = vertical_r[left] + vertical_r[right]
            pair_surface = (
                stack.surface_sheet_resistance_ohm_sq * gap / facing_width
            )
            current = stack.test_bias_v * (
                1.0 / pair_box + 1.0 / pair_surface
            )
            for matrix, value in (
                (gaps, gap),
                (box_r, pair_box),
                (surface_r, pair_surface),
                (leakage, current),
            ):
                matrix[left, right] = matrix[right, left] = value

    return SOIIsolationResult(
        mesa_area_m2=areas,
        pair_gap_m=gaps,
        box_pair_resistance_ohm=box_r,
        surface_pair_resistance_ohm=surface_r,
        pair_leakage_a=leakage,
    )
