"""Measured cross-spectral phase-noise falsification for electrode arrays."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from iqs.constants import e_C, hbar, m_He
from iqs.experiments.phase_stability import particle_velocity
from iqs.experiments.thermal_phase_noise import pairwise_differential_phase_rms


@dataclass(frozen=True)
class SpectralPhaseNoiseConfig:
    energy_eV: float = 30_000.0
    phase_budget_rad: float = 0.05
    mass_kg: float = m_He
    charge_e: float = 1.0
    minimum_bandwidth_factor: float = 10.0
    maximum_gap_factor: float = 0.2
    maximum_instrument_fraction: float = 0.2

    def __post_init__(self):
        positive = (
            self.energy_eV,
            self.phase_budget_rad,
            self.mass_kg,
            self.charge_e,
            self.minimum_bandwidth_factor,
            self.maximum_gap_factor,
            self.maximum_instrument_fraction,
        )
        if not all(np.isfinite(value) and value > 0 for value in positive):
            raise ValueError("spectral phase-noise settings must be positive")


@dataclass(frozen=True)
class SpectralPhaseNoiseResult:
    decision: str
    reasons: tuple[str, ...]
    velocity_m_s: float
    interaction_time_s: float
    frequency_count: int
    maximum_frequency_Hz: float
    required_maximum_frequency_Hz: float
    maximum_frequency_gap_Hz: float
    allowed_maximum_frequency_gap_Hz: float
    coverage_sufficient: bool
    instrument_resolved: bool
    total_worst_pair: tuple[int, int]
    total_worst_pair_rms_rad: float
    lower_bound_worst_pair: tuple[int, int]
    lower_bound_worst_pair_rms_rad: float
    instrument_worst_pair_rms_rad: float | None
    phase_budget_rad: float

    def to_dict(self):
        payload = asdict(self)
        payload["reasons"] = list(self.reasons)
        payload["total_worst_pair"] = list(self.total_worst_pair)
        payload["lower_bound_worst_pair"] = list(self.lower_bound_worst_pair)
        return payload


def _trapezoid_weights(values):
    values = np.asarray(values, dtype=float)
    weights = np.empty_like(values)
    weights[0] = 0.5 * (values[1] - values[0])
    weights[-1] = 0.5 * (values[-1] - values[-2])
    weights[1:-1] = 0.5 * (values[2:] - values[:-2])
    return weights


def _validated_spectrum(frequency_Hz, voltage_csd, electrode_count):
    frequency = np.asarray(frequency_Hz, dtype=float)
    csd = np.asarray(voltage_csd, dtype=complex)
    expected = (frequency.size, electrode_count, electrode_count)
    if frequency.ndim != 1 or frequency.size < 3:
        raise ValueError("frequency_Hz must be a 1D array with at least 3 samples")
    if csd.shape != expected:
        raise ValueError(f"voltage CSD must have shape {expected}")
    if not np.all(np.isfinite(frequency)) or not np.all(np.isfinite(csd)):
        raise ValueError("frequency and voltage CSD must be finite")
    if frequency[0] < 0 or np.any(np.diff(frequency) <= 0):
        raise ValueError("frequency_Hz must be nonnegative and strictly increasing")
    hermitian_error = np.max(np.abs(csd - np.swapaxes(csd.conj(), 1, 2)))
    scale = max(float(np.max(np.abs(csd))), np.finfo(float).tiny)
    if hermitian_error > 1e-9 * scale:
        raise ValueError("voltage CSD must be Hermitian at every frequency")
    eigenvalues = np.linalg.eigvalsh(csd)
    if eigenvalues.min() < -1e-9 * scale:
        raise ValueError("voltage CSD must be positive semidefinite")
    return frequency, csd


def phase_transfer_matrix(
    frequency_Hz,
    z_m,
    response_per_V,
    energy_eV=30_000.0,
    mass_kg=m_He,
    charge_e=1.0,
):
    """Return phase/electrode transfer ``H[f, aperture, electrode]``."""
    frequency = np.asarray(frequency_Hz, dtype=float)
    z = np.asarray(z_m, dtype=float)
    response = np.asarray(response_per_V, dtype=float)
    if frequency.ndim != 1 or not np.all(np.isfinite(frequency)):
        raise ValueError("frequency_Hz must be finite and one dimensional")
    if z.ndim != 1 or z.size < 2 or np.any(np.diff(z) <= 0):
        raise ValueError("z_m must be a strictly increasing 1D array")
    if response.ndim != 3 or response.shape[2] != z.size:
        raise ValueError("response_per_V must have shape (aperture, electrode, z)")
    if not np.all(np.isfinite(response)):
        raise ValueError("response_per_V must be finite")
    velocity = particle_velocity(energy_eV, mass_kg)
    weights = _trapezoid_weights(z)
    time_s = (z - z[0]) / velocity
    Fourier_kernel = np.exp(2j * np.pi * frequency[:, None] * time_s[None, :])
    integrated_response = np.einsum(
        "fz,oez,z->foe", Fourier_kernel, response, weights, optimize=True
    )
    return -charge_e * e_C / (hbar * velocity) * integrated_response


def phase_covariance_spectral_density(phase_transfer, voltage_csd):
    """Contract the electrode CSD into an aperture-phase CSD."""
    transfer = np.asarray(phase_transfer, dtype=complex)
    csd = np.asarray(voltage_csd, dtype=complex)
    if transfer.ndim != 3:
        raise ValueError("phase transfer must have shape (frequency, aperture, electrode)")
    if csd.shape != (
        transfer.shape[0], transfer.shape[2], transfer.shape[2]
    ):
        raise ValueError("voltage CSD and phase transfer shapes do not match")
    return np.einsum(
        "foe,fed,fpd->fop", transfer, csd, transfer.conj(), optimize=True
    )


def integrate_phase_covariance(frequency_Hz, phase_csd):
    """Integrate a one-sided phase CSD into a real covariance matrix."""
    frequency = np.asarray(frequency_Hz, dtype=float)
    spectrum = np.asarray(phase_csd, dtype=complex)
    if spectrum.ndim != 3 or spectrum.shape[0] != frequency.size:
        raise ValueError("phase CSD and frequency shapes do not match")
    covariance = np.trapezoid(spectrum, frequency, axis=0).real
    return 0.5 * (covariance + covariance.T)


def spectral_phase_covariance(
    frequency_Hz,
    voltage_csd,
    z_m,
    response_per_V,
    config: SpectralPhaseNoiseConfig,
):
    """Contract and integrate a one-sided electrode voltage CSD."""
    response = np.asarray(response_per_V)
    frequency, csd = _validated_spectrum(
        frequency_Hz, voltage_csd, response.shape[1]
    )
    transfer = phase_transfer_matrix(
        frequency,
        z_m,
        response,
        config.energy_eV,
        config.mass_kg,
        config.charge_e,
    )
    phase_csd = phase_covariance_spectral_density(transfer, csd)
    return integrate_phase_covariance(frequency, phase_csd)


def _worst_pair(pair_rms):
    upper = np.triu_indices(pair_rms.shape[0], k=1)
    flat = int(np.argmax(pair_rms[upper]))
    pair = (int(upper[0][flat]), int(upper[1][flat]))
    return pair, float(pair_rms[pair])


def analyze_spectral_phase_noise(
    frequency_Hz,
    measured_voltage_csd,
    z_m,
    response_per_V,
    config: SpectralPhaseNoiseConfig,
    instrument_voltage_csd=None,
    phase_transfer=None,
):
    """Return a conservative three-outcome verdict for a measured CSD."""
    response = np.asarray(response_per_V)
    frequency, measured = _validated_spectrum(
        frequency_Hz, measured_voltage_csd, response.shape[1]
    )
    instrument = None
    if instrument_voltage_csd is not None:
        _, instrument = _validated_spectrum(
            frequency, instrument_voltage_csd, response.shape[1]
        )
    if phase_transfer is None:
        transfer = phase_transfer_matrix(
            frequency,
            z_m,
            response,
            config.energy_eV,
            config.mass_kg,
            config.charge_e,
        )
    else:
        transfer = np.asarray(phase_transfer, dtype=complex)
        expected = (frequency.size, response.shape[0], response.shape[1])
        if transfer.shape != expected or not np.all(np.isfinite(transfer)):
            raise ValueError(f"precomputed phase transfer must have shape {expected}")
    total_covariance = integrate_phase_covariance(
        frequency,
        phase_covariance_spectral_density(transfer, measured),
    )
    total_variance = pairwise_differential_phase_rms(total_covariance) ** 2
    total_pair, total_rms = _worst_pair(np.sqrt(total_variance))

    instrument_rms = None
    if instrument is None:
        lower_variance = np.zeros_like(total_variance)
    else:
        instrument_covariance = integrate_phase_covariance(
            frequency,
            phase_covariance_spectral_density(transfer, instrument),
        )
        instrument_variance = (
            pairwise_differential_phase_rms(instrument_covariance) ** 2
        )
        instrument_rms = _worst_pair(np.sqrt(instrument_variance))[1]
        lower_variance = np.maximum(total_variance - instrument_variance, 0.0)
    lower_pair, lower_rms = _worst_pair(np.sqrt(lower_variance))

    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    interaction_time = (float(np.max(z_m)) - float(np.min(z_m))) / velocity
    required_maximum = config.minimum_bandwidth_factor / interaction_time
    allowed_gap = config.maximum_gap_factor / interaction_time
    maximum_gap = float(np.max(np.diff(frequency)))
    coverage_sufficient = bool(
        np.isclose(frequency[0], 0.0, rtol=0, atol=np.finfo(float).eps)
        and frequency[-1] >= required_maximum
        and maximum_gap <= allowed_gap
    )
    instrument_resolved = bool(
        instrument_rms is not None
        and instrument_rms
        <= config.maximum_instrument_fraction * config.phase_budget_rad
    )

    reasons = []
    if lower_rms > config.phase_budget_rad:
        decision = "falsified"
        reasons.append("instrument-subtracted lower bound exceeds the phase budget")
    elif coverage_sufficient and instrument_resolved and total_rms <= config.phase_budget_rad:
        decision = "not_falsified"
        reasons.append("resolved total spectrum remains within the phase budget")
    else:
        decision = "inconclusive"
        if not coverage_sufficient:
            reasons.append("frequency coverage cannot support a pass")
        if not instrument_resolved:
            reasons.append("instrument spectrum does not resolve the phase budget")
        if total_rms > config.phase_budget_rad and lower_rms <= config.phase_budget_rad:
            reasons.append("result overlaps the instrument-limited decision boundary")

    return SpectralPhaseNoiseResult(
        decision=decision,
        reasons=tuple(reasons),
        velocity_m_s=velocity,
        interaction_time_s=interaction_time,
        frequency_count=int(frequency.size),
        maximum_frequency_Hz=float(frequency[-1]),
        required_maximum_frequency_Hz=required_maximum,
        maximum_frequency_gap_Hz=maximum_gap,
        allowed_maximum_frequency_gap_Hz=allowed_gap,
        coverage_sufficient=coverage_sufficient,
        instrument_resolved=instrument_resolved,
        total_worst_pair=total_pair,
        total_worst_pair_rms_rad=total_rms,
        lower_bound_worst_pair=lower_pair,
        lower_bound_worst_pair_rms_rad=lower_rms,
        instrument_worst_pair_rms_rad=instrument_rms,
        phase_budget_rad=config.phase_budget_rad,
    )


__all__ = [
    "SpectralPhaseNoiseConfig",
    "SpectralPhaseNoiseResult",
    "analyze_spectral_phase_noise",
    "integrate_phase_covariance",
    "phase_covariance_spectral_density",
    "phase_transfer_matrix",
    "spectral_phase_covariance",
]
