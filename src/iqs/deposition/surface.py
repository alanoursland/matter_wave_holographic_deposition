"""Minimal stochastic, multi-material surface evolution.

The model is intentionally 2.5D: every material is represented by a local
thickness field over one advancing surface. It captures finite dose, sticking,
registration, projection blur, and post-landing diffusion without claiming
atomistic crystallization or interface chemistry.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter, shift


@dataclass(frozen=True)
class DepositionMaterial:
    name: str
    atomic_volume_m3: float

    def __post_init__(self):
        if not self.name:
            raise ValueError("material name must not be empty")
        if (not np.isfinite(self.atomic_volume_m3)
                or self.atomic_volume_m3 <= 0):
            raise ValueError("atomic volume must be positive")


@dataclass
class SurfaceState:
    shape: tuple[int, int]
    pixel_pitch_m: float
    thickness_m: dict[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self):
        if len(self.shape) != 2 or min(self.shape) < 2:
            raise ValueError("surface shape must be two-dimensional")
        if not np.isfinite(self.pixel_pitch_m) or self.pixel_pitch_m <= 0:
            raise ValueError("pixel pitch must be positive")
        for name, values in self.thickness_m.items():
            array = np.asarray(values, dtype=float)
            if array.shape != self.shape or not np.all(np.isfinite(array)):
                raise ValueError(f"invalid thickness field for {name!r}")
            if np.any(array < 0):
                raise ValueError("material thickness cannot be negative")
            self.thickness_m[name] = array.copy()

    @property
    def total_height_m(self):
        if not self.thickness_m:
            return np.zeros(self.shape, dtype=float)
        return np.sum(list(self.thickness_m.values()), axis=0)

    def material_thickness(self, name):
        return self.thickness_m.get(name, np.zeros(self.shape, dtype=float))


@dataclass(frozen=True)
class DepositionResult:
    deposited_thickness_m: np.ndarray
    registered_target_m: np.ndarray
    attempted_atoms: float
    retained_atoms: int
    mean_retention_probability: float


def deposit_layer(
    surface: SurfaceState,
    material: DepositionMaterial,
    target_thickness_m,
    *,
    nominal_retention_probability: float,
    actual_retention_probability,
    projection_sigma_m: float = 0.0,
    diffusion_sigma_m: float = 0.0,
    registration_offset_m: tuple[float, float] = (0.0, 0.0),
    rng: np.random.Generator | None = None,
) -> DepositionResult:
    """Deposit one calibrated stochastic layer onto ``surface``.

    Commanded dose is calibrated for ``nominal_retention_probability``.
    Spatially varying actual retention represents material/interface sticking
    and post-landing survival. Gaussian projection blur acts on attempted dose;
    diffusion acts on retained material and conserves deposited volume.
    """
    target = np.asarray(target_thickness_m, dtype=float)
    if target.shape != surface.shape or not np.all(np.isfinite(target)):
        raise ValueError("target thickness must be finite and match the surface")
    if np.any(target < 0):
        raise ValueError("target thickness cannot be negative")
    if (not np.isfinite(nominal_retention_probability)
            or not 0 < nominal_retention_probability <= 1):
        raise ValueError("nominal retention probability must lie in (0, 1]")
    retention = np.asarray(actual_retention_probability, dtype=float)
    try:
        retention = np.broadcast_to(retention, surface.shape).copy()
    except ValueError as exc:
        raise ValueError("actual retention is not broadcastable to surface") from exc
    if (not np.all(np.isfinite(retention)) or np.any(retention < 0)
            or np.any(retention > 1)):
        raise ValueError("actual retention probabilities must lie in [0, 1]")
    if projection_sigma_m < 0 or diffusion_sigma_m < 0:
        raise ValueError("blur and diffusion sigmas cannot be negative")
    if len(registration_offset_m) != 2 or not np.all(
            np.isfinite(registration_offset_m)):
        raise ValueError("registration offset must be a finite pair")

    rng = np.random.default_rng() if rng is None else rng
    shift_pixels = tuple(
        value / surface.pixel_pitch_m for value in registration_offset_m
    )
    registered = shift(
        target,
        shift=shift_pixels,
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    )
    sigma_projection_pixels = projection_sigma_m / surface.pixel_pitch_m
    if sigma_projection_pixels > 0:
        registered = gaussian_filter(
            registered, sigma_projection_pixels, mode="constant")

    pixel_area = surface.pixel_pitch_m ** 2
    commanded_atoms = (
        registered * pixel_area
        / material.atomic_volume_m3
        / nominal_retention_probability
    )
    retained_mean = commanded_atoms * retention
    retained_counts = rng.poisson(np.maximum(retained_mean, 0.0))
    deposited = retained_counts * material.atomic_volume_m3 / pixel_area
    sigma_diffusion_pixels = diffusion_sigma_m / surface.pixel_pitch_m
    if sigma_diffusion_pixels > 0:
        deposited = gaussian_filter(
            deposited, sigma_diffusion_pixels, mode="constant")

    existing = surface.material_thickness(material.name)
    surface.thickness_m[material.name] = existing + deposited
    attempted = float(commanded_atoms.sum())
    weighted_retention = float(
        np.sum(commanded_atoms * retention) / (attempted + 1e-30)
    )
    return DepositionResult(
        deposited_thickness_m=deposited,
        registered_target_m=registered,
        attempted_atoms=attempted,
        retained_atoms=int(retained_counts.sum()),
        mean_retention_probability=weighted_retention,
    )
