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
