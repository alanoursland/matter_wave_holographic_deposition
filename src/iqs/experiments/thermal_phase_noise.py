"""Thermal RC voltage-noise floor for electrostatic matter-wave phase plates.

For a resistor ``R`` at temperature ``T`` driving an electrode capacitance
``C``, the equilibrium voltage autocovariance is

    <V(t)V(t+u)> = (k_B T / C) exp(-|u| / RC).

An ion integrates that voltage during its phase-critical flight time ``tau``.
The exact variance of the integrated voltage is

    2 (k_B T / C) [tau RC - (RC)^2 (1 - exp(-tau/RC))].

This is the time-domain equivalent of filtering the Johnson spectrum by the
rectangular transit window. A field-resolved variant also accepts the actual
dimensionless potential response ``g(z)`` of an electrode. Two equal path
channels with correlation ``rho`` multiply the single-channel variance by
``2(1-rho)``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from iqs.constants import e_C, h, hbar, k_B, m_He
from iqs.experiments.phase_stability import particle_velocity


@dataclass(frozen=True)
class ThermalPhaseNoiseConfig:
    """Circuit and beam settings for the thermal phase-noise gate."""

    energy_eV: float = 30_000.0
    path_length_m: float = 30e-6
    temperature_K: float = 4.0
    resistance_ohm: float = 50.0
    capacitance_F: float = 100e-15
    path_correlation: float = 0.0
    phase_budget_rad: float = 0.05
    mass_kg: float = m_He
    charge_e: float = 1.0

    def __post_init__(self):
        positive = (
            self.energy_eV,
            self.path_length_m,
            self.temperature_K,
            self.resistance_ohm,
            self.capacitance_F,
            self.phase_budget_rad,
            self.mass_kg,
            self.charge_e,
        )
        if not all(np.isfinite(value) and value > 0 for value in positive):
            raise ValueError("thermal phase-noise settings must be positive")
        if not np.isfinite(self.path_correlation) or not -1 <= self.path_correlation <= 1:
            raise ValueError("path_correlation must lie in [-1, 1]")


@dataclass(frozen=True)
class ThermalPhaseNoiseResult:
    decision: str
    coupling_model: str
    velocity_m_s: float
    interaction_length_m: float
    flight_time_s: float
    dc_phase_gain_rad_per_V: float
    rc_time_s: float
    transit_to_rc_ratio: float
    thermal_phase_rms_rad: float
    phase_budget_rad: float
    budget_ratio: float
    maximum_resistance_ohm: float
    required_path_correlation: float
    maximum_temperature_K: float
    quantum_crossover_temperature_K: float
    classical_temperature_ratio: float

    def to_dict(self):
        return asdict(self)


def _stable_rc_bracket(transit_time_s, rc_time_s):
    """Return ``tau*RC - RC^2*(1-exp(-tau/RC))`` stably."""
    tau = np.asarray(transit_time_s, dtype=float)
    rc = np.asarray(rc_time_s, dtype=float)
    if np.any(tau <= 0) or np.any(rc <= 0):
        raise ValueError("transit and RC times must be positive")
    x = tau / rc
    direct = x + np.expm1(-x)
    series = x**2 / 2 - x**3 / 6 + x**4 / 24 - x**5 / 120
    dimensionless = np.where(x < 1e-3, series, direct)
    return rc**2 * dimensionless


def integrated_voltage_variance(
    transit_time_s,
    temperature_K,
    resistance_ohm,
    capacitance_F,
):
    """Exact single-electrode variance of ``integral V(t) dt``."""
    tau = np.asarray(transit_time_s, dtype=float)
    temperature = np.asarray(temperature_K, dtype=float)
    resistance = np.asarray(resistance_ohm, dtype=float)
    capacitance = np.asarray(capacitance_F, dtype=float)
    if (np.any(tau <= 0) or np.any(temperature < 0)
            or np.any(resistance < 0) or np.any(capacitance <= 0)):
        raise ValueError("invalid thermal RC parameters")
    rc = resistance * capacitance
    broadcast = np.broadcast_arrays(tau, temperature, resistance, capacitance)
    tau, temperature, resistance, capacitance = broadcast
    rc = resistance * capacitance
    result = np.zeros_like(rc, dtype=float)
    active = (temperature > 0) & (resistance > 0)
    if np.any(active):
        bracket = _stable_rc_bracket(tau[active], rc[active])
        result[active] = 2 * k_B * temperature[active] / capacitance[active] * bracket
    return float(result) if result.ndim == 0 else result


def _validated_response_profile(z_m, response_per_V):
    z = np.asarray(z_m, dtype=float)
    response = np.asarray(response_per_V, dtype=float)
    if z.ndim != 1 or response.ndim != 1 or z.shape != response.shape:
        raise ValueError("z_m and response_per_V must be equal-length 1D arrays")
    if z.size < 2 or not np.all(np.isfinite(z)) or not np.all(np.isfinite(response)):
        raise ValueError("field response profile must contain finite samples")
    if np.any(np.diff(z) <= 0):
        raise ValueError("z_m must be strictly increasing")
    return z, response


def weighted_integrated_voltage_variance(
    z_m,
    response_per_V,
    velocity_m_s,
    temperature_K,
    resistance_ohm,
    capacitance_F,
):
    """Variance of ``integral g(z(t)) V(t) dt`` for an RC noise source.

    The response is treated as piecewise constant at each spatial cell's
    endpoint average. The exponential covariance is integrated exactly over
    every pair of cells, retaining both the white-noise and quasistatic limits
    without requiring a time grid finer than ``R*C``.
    """
    z, response = _validated_response_profile(z_m, response_per_V)
    values = (velocity_m_s, temperature_K, resistance_ohm, capacitance_F)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("thermal response settings must be finite")
    if velocity_m_s <= 0 or temperature_K < 0 or resistance_ohm < 0 or capacitance_F <= 0:
        raise ValueError("invalid thermal response settings")
    if temperature_K == 0 or resistance_ohm == 0 or not np.any(response):
        return 0.0

    cell_response, covariance_integrals = _response_cells_and_rc_covariance(
        z, response, velocity_m_s, resistance_ohm, capacitance_F
    )
    variance = (
        k_B * temperature_K / capacitance_F
        * cell_response @ covariance_integrals @ cell_response
    )
    return float(max(variance, 0.0))


def _response_cells_and_rc_covariance(
    z_m,
    response_per_V,
    velocity_m_s,
    resistance_ohm,
    capacitance_F,
):
    """Return piecewise response values and exactly integrated RC kernel."""
    edges = (z_m - z_m[0]) / velocity_m_s
    widths = np.diff(edges)
    cell_response = 0.5 * (response_per_V[:-1] + response_per_V[1:])
    rc = resistance_ohm * capacitance_F
    scaled_widths = widths / rc
    covariance_integrals = np.zeros((widths.size, widths.size), dtype=float)
    diagonal = 2 * rc**2 * (scaled_widths + np.expm1(-scaled_widths))
    np.fill_diagonal(covariance_integrals, diagonal)
    edge_factors = -np.expm1(-scaled_widths)
    for left in range(widths.size - 1):
        right = np.arange(left + 1, widths.size)
        gap = edges[right] - edges[left + 1]
        values = (
            rc**2
            * edge_factors[left]
            * edge_factors[right]
            * np.exp(-gap / rc)
        )
        covariance_integrals[left, right] = values
        covariance_integrals[right, left] = values
    return cell_response, covariance_integrals


def multi_electrode_phase_covariance(
    config: ThermalPhaseNoiseConfig,
    z_m,
    response_per_V,
    electrode_correlation=None,
):
    """Return the aperture-phase covariance from multiple electrode channels.

    ``response_per_V`` has shape ``(readout, electrode, z)``. Every electrode
    uses the same R, C, and T in this baseline model. ``electrode_correlation``
    is a frequency-independent channel correlation matrix; identity is used
    by default. The output covariance is in radian squared.
    """
    z = np.asarray(z_m, dtype=float)
    response = np.asarray(response_per_V, dtype=float)
    if response.ndim != 3 or response.shape[2:] != z.shape:
        raise ValueError("response_per_V must have shape (readout, electrode, z)")
    _validated_response_profile(z, response[0, 0])
    if not np.all(np.isfinite(response)):
        raise ValueError("multi-electrode response must be finite")
    electrode_count = response.shape[1]
    if electrode_correlation is None:
        correlation = np.eye(electrode_count)
    else:
        correlation = np.asarray(electrode_correlation, dtype=float)
    if correlation.shape != (electrode_count, electrode_count):
        raise ValueError("electrode correlation matrix has the wrong shape")
    if not np.all(np.isfinite(correlation)) or not np.allclose(
        correlation, correlation.T, rtol=1e-10, atol=1e-12
    ):
        raise ValueError("electrode correlation matrix must be finite and symmetric")
    if not np.allclose(np.diag(correlation), 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError("electrode correlation matrix must have a unit diagonal")
    if np.linalg.eigvalsh(correlation).min() < -1e-10:
        raise ValueError("electrode correlation matrix must be positive semidefinite")
    if config.temperature_K == 0 or config.resistance_ohm == 0:
        return np.zeros((response.shape[0], response.shape[0]))

    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    first_cells, rc_kernel = _response_cells_and_rc_covariance(
        z,
        response[0, 0],
        velocity,
        config.resistance_ohm,
        config.capacitance_F,
    )
    del first_cells
    cells = 0.5 * (response[:, :, :-1] + response[:, :, 1:])
    voltage_covariance = (
        k_B * config.temperature_K / config.capacitance_F
        * np.einsum(
            "oea,ef,ab,pfb->op",
            cells,
            correlation,
            rc_kernel,
            cells,
            optimize=True,
        )
    )
    phase_factor = config.charge_e * e_C / hbar
    phase_covariance = phase_factor**2 * voltage_covariance
    return 0.5 * (phase_covariance + phase_covariance.T)


def pairwise_differential_phase_rms(phase_covariance):
    """Return the RMS matrix for all pairwise phase differences."""
    covariance = np.asarray(phase_covariance, dtype=float)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("phase covariance must be square")
    variances = (
        np.diag(covariance)[:, None]
        + np.diag(covariance)[None, :]
        - 2 * covariance
    )
    return np.sqrt(np.maximum(variances, 0.0))


def thermal_phase_rms(
    config: ThermalPhaseNoiseConfig,
    z_m=None,
    response_per_V=None,
) -> float:
    """RMS differential phase caused by classical equilibrium RC noise."""
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    if (z_m is None) != (response_per_V is None):
        raise ValueError("z_m and response_per_V must be supplied together")
    if z_m is None:
        tau = config.path_length_m / velocity
        single_variance = integrated_voltage_variance(
            tau,
            config.temperature_K,
            config.resistance_ohm,
            config.capacitance_F,
        )
    else:
        single_variance = weighted_integrated_voltage_variance(
            z_m,
            response_per_V,
            velocity,
            config.temperature_K,
            config.resistance_ohm,
            config.capacitance_F,
        )
    differential_factor = 2.0 * (1.0 - config.path_correlation)
    charge = config.charge_e * e_C
    return float(charge / hbar * np.sqrt(differential_factor * single_variance))


def maximum_resistance_for_budget(
    config: ThermalPhaseNoiseConfig,
    z_m=None,
    response_per_V=None,
) -> float:
    """Largest Thevenin resistance consistent with the phase budget."""
    if config.path_correlation == 1:
        return float("inf")

    def sigma_at(log_resistance):
        candidate = ThermalPhaseNoiseConfig(
            **{**asdict(config), "resistance_ohm": 10.0**log_resistance}
        )
        return thermal_phase_rms(candidate, z_m, response_per_V)

    low, high = -18.0, 18.0
    if sigma_at(high) <= config.phase_budget_rad:
        return float("inf")
    for _ in range(100):
        middle = 0.5 * (low + high)
        if sigma_at(middle) <= config.phase_budget_rad:
            low = middle
        else:
            high = middle
    return float(10.0**low)


def required_path_correlation(
    config: ThermalPhaseNoiseConfig,
    z_m=None,
    response_per_V=None,
) -> float:
    """Minimum equal-channel correlation needed at the configured R, C, T."""
    uncorrelated = ThermalPhaseNoiseConfig(
        **{**asdict(config), "path_correlation": 0.0}
    )
    sigma_zero = thermal_phase_rms(uncorrelated, z_m, response_per_V)
    if sigma_zero <= config.phase_budget_rad:
        return 0.0
    required = 1.0 - (config.phase_budget_rad / sigma_zero) ** 2
    return float(np.clip(required, 0.0, 1.0))


def maximum_temperature_for_budget(
    config: ThermalPhaseNoiseConfig,
    z_m=None,
    response_per_V=None,
) -> float:
    """Classical temperature ceiling at fixed circuit and correlation."""
    sigma = thermal_phase_rms(config, z_m, response_per_V)
    if sigma == 0:
        return float("inf")
    return float(config.temperature_K * (config.phase_budget_rad / sigma) ** 2)


def analyze_thermal_phase_noise(
    config: ThermalPhaseNoiseConfig,
    z_m=None,
    response_per_V=None,
) -> ThermalPhaseNoiseResult:
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    if (z_m is None) != (response_per_V is None):
        raise ValueError("z_m and response_per_V must be supplied together")
    if z_m is None:
        interaction_length = config.path_length_m
        coupling_model = "rectangular full-coupling approximation"
        dc_gain = -config.charge_e * e_C / hbar * config.path_length_m / velocity
    else:
        z, response = _validated_response_profile(z_m, response_per_V)
        interaction_length = float(z[-1] - z[0])
        coupling_model = "field-weighted aperture response"
        dc_gain = float(
            -config.charge_e * e_C / (hbar * velocity)
            * np.trapezoid(response, z)
        )
    tau = interaction_length / velocity
    rc_time = config.resistance_ohm * config.capacitance_F
    sigma = thermal_phase_rms(config, z_m, response_per_V)
    ratio = sigma / config.phase_budget_rad
    crossover = h / (k_B * tau)
    return ThermalPhaseNoiseResult(
        decision="falsified" if ratio > 1 else "not_falsified",
        coupling_model=coupling_model,
        velocity_m_s=velocity,
        interaction_length_m=interaction_length,
        flight_time_s=tau,
        dc_phase_gain_rad_per_V=dc_gain,
        rc_time_s=rc_time,
        transit_to_rc_ratio=tau / rc_time,
        thermal_phase_rms_rad=sigma,
        phase_budget_rad=config.phase_budget_rad,
        budget_ratio=ratio,
        maximum_resistance_ohm=maximum_resistance_for_budget(
            config, z_m, response_per_V),
        required_path_correlation=required_path_correlation(
            config, z_m, response_per_V),
        maximum_temperature_K=maximum_temperature_for_budget(
            config, z_m, response_per_V),
        quantum_crossover_temperature_K=crossover,
        classical_temperature_ratio=config.temperature_K / crossover,
    )


__all__ = [
    "ThermalPhaseNoiseConfig",
    "ThermalPhaseNoiseResult",
    "analyze_thermal_phase_noise",
    "integrated_voltage_variance",
    "maximum_resistance_for_budget",
    "maximum_temperature_for_budget",
    "multi_electrode_phase_covariance",
    "pairwise_differential_phase_rms",
    "required_path_correlation",
    "thermal_phase_rms",
    "weighted_integrated_voltage_variance",
]
