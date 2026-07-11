"""Physical models for programmable matter-wave phase actuators.

The holography solver needs a local transmission ``exp(i * phi(x, y))``.
That is a valid mathematical actuator contract, but magnetic flux assigned
to a loop does not automatically provide that transmission.  This module
keeps the ideal contract separate from an electromagnetic loop model.

For a filamentary loop lying in the x-y plane, every current element
``dl'`` has zero z component.  In Coulomb gauge,

    A(r) = mu0 I / (4 pi) integral(dl' / |r-r'|),

so ``A_z`` is identically zero.  A straight normally incident ray has
``dl = dz z_hat`` and therefore accumulates no axial thin-screen phase.
Observable Aharonov-Bohm phases can still occur in geometries whose paths
form a loop around confined flux, but that is not the local pixel model.
"""

from dataclasses import dataclass

import numpy as np

from iqs.constants import e_C, hbar, mu_0


class AchromaticPhaseResponse:
    """Velocity-independent phase response for an ideal or AB plate."""

    name = "achromatic"
    periodic_controls = True

    @staticmethod
    def phase_scales(wavelengths, lambda0):
        wavelengths = _validated_wavelengths(wavelengths, lambda0)
        return np.ones_like(wavelengths)


class ElectrostaticPhaseResponse:
    """Thin static-potential phase response in the straight-ray limit.

    For a fixed potential integral, ``phi = -q integral(V dl)/(hbar v)``.
    Nonrelativistically ``lambda = h/(m v)``, hence a control phase defined
    at ``lambda0`` scales as ``phi(lambda) = phi0 * lambda/lambda0``.
    """

    name = "electrostatic"
    periodic_controls = False

    @staticmethod
    def phase_scales(wavelengths, lambda0):
        wavelengths = _validated_wavelengths(wavelengths, lambda0)
        return wavelengths / float(lambda0)


@dataclass(frozen=True)
class ElectrostaticValidityReport:
    """Thin-screen validity and hardware requirements for one phase map."""

    ok: bool
    max_voltage_V: float
    rms_voltage_V: float
    max_energy_ratio: float
    max_deflection_rad: float
    rms_deflection_rad: float
    max_walkoff_m: float
    max_walkoff_pitch_ratio: float
    max_transverse_field_V_m: float
    phase_span_rad: float
    violations: tuple[str, ...]

    def require(self) -> "ElectrostaticValidityReport":
        if not self.ok:
            raise PhaseAuthorityError(
                "electrostatic thin-screen gate failed: "
                + "; ".join(self.violations)
            )
        return self


@dataclass(frozen=True)
class ElectrostaticPlateGeometry:
    """Physical geometry used to validate an electrostatic phase plate.

    ``pixel_pitch_m`` and ``interaction_length_m`` belong to the physical
    control plane, not a demagnified image-side simulation grid.
    """

    pixel_pitch_m: float
    interaction_length_m: float
    kinetic_energy_eV: float
    particle_mass_kg: float
    charge_C: float = e_C
    max_energy_ratio: float = 0.1
    max_walkoff_pitch_ratio: float = 0.1
    max_deflection_rad: float = 0.1
    voltage_limit_V: float | None = None
    transverse_field_limit_V_m: float | None = None

    def __post_init__(self):
        positive = (
            self.pixel_pitch_m, self.interaction_length_m,
            self.kinetic_energy_eV, self.particle_mass_kg,
            abs(self.charge_C), self.max_energy_ratio,
            self.max_walkoff_pitch_ratio, self.max_deflection_rad,
        )
        if not all(np.isfinite(value) and value > 0 for value in positive):
            raise ValueError("electrostatic plate parameters must be positive")
        optional = (self.voltage_limit_V, self.transverse_field_limit_V_m)
        if not all(value is None or (np.isfinite(value) and value > 0)
                   for value in optional):
            raise ValueError("optional hardware limits must be positive")

    def evaluate(self, phase_controls) -> ElectrostaticValidityReport:
        """Evaluate voltage, kick, and finite-thickness walkoff.

        A spatially uniform phase is removed because it is globally
        unobservable and can be supplied by the reference electrode. The
        remaining voltage map is the minimum-peak realization obtained by
        centering the phase range.
        """
        phase = np.asarray(phase_controls, dtype=float)
        if phase.ndim != 2 or min(phase.shape) < 2:
            raise ValueError("phase_controls must be a 2D array at least 2x2")
        if not np.all(np.isfinite(phase)):
            raise ValueError("phase_controls must be finite")

        energy_J = self.kinetic_energy_eV * e_C
        velocity = np.sqrt(2 * energy_J / self.particle_mass_kg)
        k0 = self.particle_mass_kg * velocity / hbar

        offset = 0.5 * (float(phase.max()) + float(phase.min()))
        phase_centered = phase - offset
        voltage = (
            -phase_centered * hbar * velocity
            / (self.charge_C * self.interaction_length_m)
        )

        edge_order = 2 if min(phase.shape) >= 3 else 1
        grad_y, grad_x = np.gradient(
            phase_centered, self.pixel_pitch_m,
            self.pixel_pitch_m, edge_order=edge_order)
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)
        deflection = grad_mag / k0
        walkoff = deflection * self.interaction_length_m

        field_y, field_x = np.gradient(
            voltage, self.pixel_pitch_m,
            self.pixel_pitch_m, edge_order=edge_order)
        field_mag = np.sqrt(field_x**2 + field_y**2)
        energy_ratio = np.abs(self.charge_C * voltage) / energy_J

        max_voltage = float(np.abs(voltage).max())
        max_field = float(field_mag.max())
        max_energy = float(energy_ratio.max())
        max_theta = float(deflection.max())
        max_walkoff = float(walkoff.max())
        walkoff_ratio = max_walkoff / self.pixel_pitch_m
        violations = []
        if max_energy > self.max_energy_ratio:
            violations.append(
                f"|qV|/E={max_energy:.3g} > {self.max_energy_ratio:.3g}")
        if max_theta > self.max_deflection_rad:
            violations.append(
                f"theta={max_theta:.3g} rad > "
                f"{self.max_deflection_rad:.3g} rad")
        if walkoff_ratio > self.max_walkoff_pitch_ratio:
            violations.append(
                f"walkoff/pitch={walkoff_ratio:.3g} > "
                f"{self.max_walkoff_pitch_ratio:.3g}")
        if (self.voltage_limit_V is not None
                and max_voltage > self.voltage_limit_V):
            violations.append(
                f"|V|={max_voltage:.3g} V > "
                f"{self.voltage_limit_V:.3g} V")
        if (self.transverse_field_limit_V_m is not None
                and max_field > self.transverse_field_limit_V_m):
            violations.append(
                f"|E_perp|={max_field:.3g} V/m > "
                f"{self.transverse_field_limit_V_m:.3g} V/m")

        return ElectrostaticValidityReport(
            ok=not violations,
            max_voltage_V=max_voltage,
            rms_voltage_V=float(np.sqrt(np.mean(voltage**2))),
            max_energy_ratio=max_energy,
            max_deflection_rad=max_theta,
            rms_deflection_rad=float(np.sqrt(np.mean(deflection**2))),
            max_walkoff_m=max_walkoff,
            max_walkoff_pitch_ratio=float(walkoff_ratio),
            max_transverse_field_V_m=max_field,
            phase_span_rad=float(np.ptp(phase)),
            violations=tuple(violations),
        )


def _validated_wavelengths(wavelengths, lambda0):
    wavelengths = np.asarray(wavelengths, dtype=float).reshape(-1)
    lambda0 = float(lambda0)
    if wavelengths.size == 0 or not np.all(np.isfinite(wavelengths)):
        raise ValueError("wavelengths must be a non-empty finite array")
    if np.any(wavelengths <= 0) or not np.isfinite(lambda0) or lambda0 <= 0:
        raise ValueError("wavelengths and lambda0 must be positive")
    return wavelengths


def resolve_phase_response(response):
    """Resolve a response name or validate a response object."""
    if isinstance(response, str):
        key = response.strip().lower()
        if key in {"ideal", "ab", "achromatic"}:
            return AchromaticPhaseResponse()
        if key in {"electrostatic", "electric"}:
            return ElectrostaticPhaseResponse()
        raise ValueError(
            "phase response must be 'achromatic' or 'electrostatic'"
        )
    if not hasattr(response, "name") or not callable(
            getattr(response, "phase_scales", None)):
        raise TypeError("phase response must define name and phase_scales()")
    if not hasattr(response, "periodic_controls"):
        raise TypeError("phase response must define periodic_controls")
    return response


class PhaseAuthorityError(RuntimeError):
    """Raised when an actuator cannot supply the requested phase span."""


@dataclass(frozen=True)
class PhaseAuthorityReport:
    """Result of mapping an actuator geometry to a local phase screen."""

    model: str
    phase_span_rad: float
    phase_rms_rad: float
    required_span_rad: float
    local_thin_screen: bool
    field_free_verified: bool | None
    reason: str

    @property
    def ok(self) -> bool:
        return self.local_thin_screen and self.phase_span_rad >= self.required_span_rad

    def require(self) -> "PhaseAuthorityReport":
        """Return this report or raise when the authority gate fails."""
        if not self.ok:
            raise PhaseAuthorityError(
                f"{self.model} has phase span {self.phase_span_rad:.3e} rad; "
                f"{self.required_span_rad:.3e} rad is required. {self.reason}"
            )
        return self


class IdealPhasePlate:
    """Explicit ideal phase-only actuator used by the inverse solver."""

    def __init__(self, phase_screen):
        phase = np.asarray(phase_screen, dtype=float)
        if phase.ndim != 2:
            raise ValueError("phase_screen must be a two-dimensional array")
        if not np.all(np.isfinite(phase)):
            raise ValueError("phase_screen must contain only finite values")
        self.phase_screen = phase

    def transmission(self) -> np.ndarray:
        return np.exp(1j * self.phase_screen)

    def authority_report(self, required_span_rad=np.pi) -> PhaseAuthorityReport:
        centered = self.phase_screen - self.phase_screen.mean()
        return PhaseAuthorityReport(
            model="ideal-phase-plate",
            phase_span_rad=float(np.ptp(self.phase_screen)),
            phase_rms_rad=float(np.sqrt(np.mean(centered**2))),
            required_span_rad=float(required_span_rad),
            local_thin_screen=True,
            field_free_verified=None,
            reason=(
                "This is an ideal actuator contract; electromagnetic and "
                "fabrication feasibility are intentionally unspecified."
            ),
        )


class CoplanarSquareLoopArray:
    """Filamentary square current loops in one x-y plane.

    The vector-potential calculation is useful for arbitrary trajectories.
    The axial thin-screen result is analytic: it is zero for normally
    incident straight rays because this geometry has ``A_z = 0`` in the
    Coulomb-gauge Biot-Savart representation.
    """

    def __init__(self, centers_xy, side_length, z_plane=0.0,
                 n_segments_per_side=24, wire_radius=None):
        centers = np.asarray(centers_xy, dtype=float)
        if centers.ndim != 2 or centers.shape[1] != 2:
            raise ValueError("centers_xy must have shape (n_loops, 2)")
        if side_length <= 0:
            raise ValueError("side_length must be positive")
        if n_segments_per_side < 2:
            raise ValueError("n_segments_per_side must be at least 2")

        self.centers_xy = centers
        self.side_length = float(side_length)
        self.z_plane = float(z_plane)
        self.n_segments_per_side = int(n_segments_per_side)
        self.wire_radius = float(
            wire_radius if wire_radius is not None else side_length * 1e-3
        )
        if self.wire_radius <= 0:
            raise ValueError("wire_radius must be positive")

        self._local_midpoints, self._local_dl = self._square_segments()

    def _square_segments(self):
        half = self.side_length / 2
        corners = np.array([
            [-half, -half, 0.0],
            [half, -half, 0.0],
            [half, half, 0.0],
            [-half, half, 0.0],
            [-half, -half, 0.0],
        ])
        midpoints = []
        elements = []
        for start, stop in zip(corners[:-1], corners[1:]):
            dl = (stop - start) / self.n_segments_per_side
            t = (np.arange(self.n_segments_per_side) + 0.5)
            midpoints.append(start + t[:, None] * dl)
            elements.append(np.repeat(dl[None, :], self.n_segments_per_side, axis=0))
        return np.vstack(midpoints), np.vstack(elements)

    def vector_potential(self, points_xyz, currents_A) -> np.ndarray:
        """Return Coulomb-gauge vector potential at arbitrary points.

        The filament singularity is regularized by ``wire_radius``.  Open
        path phases computed from this gauge are not independently
        observable; closed-path differences are gauge invariant.
        """
        points = np.asarray(points_xyz, dtype=float)
        original_shape = points.shape
        if points.ndim < 1 or original_shape[-1] != 3:
            raise ValueError("points_xyz must have final dimension 3")
        flat_points = points.reshape(-1, 3)

        currents = np.asarray(currents_A, dtype=float).reshape(-1)
        if currents.size != self.centers_xy.shape[0]:
            raise ValueError("one current is required per loop center")

        result = np.zeros_like(flat_points)
        prefactor = mu_0 / (4 * np.pi)
        for center, current in zip(self.centers_xy, currents):
            if current == 0:
                continue
            shift = np.array([center[0], center[1], self.z_plane])
            wire_points = self._local_midpoints + shift
            separation = flat_points[:, None, :] - wire_points[None, :, :]
            distance = np.sqrt(
                np.sum(separation**2, axis=-1) + self.wire_radius**2
            )
            result += prefactor * current * np.sum(
                self._local_dl[None, :, :] / distance[:, :, None], axis=1
            )
        return result.reshape(original_shape)

    def coulomb_gauge_path_phase(self, paths_xyz, currents_A,
                                 charge_C=e_C) -> np.ndarray:
        """Integrate ``q A dot dl / hbar`` along one or more paths.

        ``paths_xyz`` has shape ``(..., n_points, 3)``.  Use this method to
        construct closed contours or path differences; an isolated open-path
        value is gauge dependent.
        """
        paths = np.asarray(paths_xyz, dtype=float)
        if paths.ndim < 2 or paths.shape[-1] != 3 or paths.shape[-2] < 2:
            raise ValueError("paths_xyz must have shape (..., n_points, 3)")
        dl = np.diff(paths, axis=-2)
        midpoints = 0.5 * (paths[..., 1:, :] + paths[..., :-1, :])
        potential = self.vector_potential(midpoints, currents_A)
        integral = np.sum(potential * dl, axis=(-2, -1))
        return charge_C * integral / hbar

    def axial_phase_screen(self, x, y, z_start, z_stop, currents_A,
                           charge_C=e_C) -> np.ndarray:
        """Return the phase for straight z-directed rays.

        The arguments are retained to make the physical contract explicit.
        For coplanar filamentary loops the result is identically zero; no
        expensive field sampling is necessary.
        """
        x_arr, y_arr = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        )
        currents = np.asarray(currents_A, dtype=float).reshape(-1)
        if currents.size != self.centers_xy.shape[0]:
            raise ValueError("one current is required per loop center")
        if not np.isfinite(z_start) or not np.isfinite(z_stop):
            raise ValueError("z_start and z_stop must be finite")
        if not np.isfinite(charge_C):
            raise ValueError("charge_C must be finite")
        return np.zeros_like(x_arr, dtype=float)

    def authority_report(self, currents_A, required_span_rad=np.pi
                         ) -> PhaseAuthorityReport:
        currents = np.asarray(currents_A, dtype=float).reshape(-1)
        if currents.size != self.centers_xy.shape[0]:
            raise ValueError("one current is required per loop center")
        return PhaseAuthorityReport(
            model="coplanar-square-loop-array/normal-incidence",
            phase_span_rad=0.0,
            phase_rms_rad=0.0,
            required_span_rad=float(required_span_rad),
            local_thin_screen=False,
            field_free_verified=False,
            reason=(
                "All loop current elements lie in the x-y plane, so A_z is "
                "zero and straight z-directed rays acquire no local phase. "
                "The loop field is also not confined away from the beam."
            ),
        )
