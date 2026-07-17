"""Falsification analysis for electrostatic matter-wave phase stability.

The programmable phase is calibrated at the start of an exposure cycle.
Only subsequent differential-potential drift is uncorrected.  For a slowly
varying effective path differential ``delta_V`` over a phase-critical leg,

    delta_phi = q * delta_V * (L / v) / hbar.

This module converts a measured differential-voltage time series into that
phase error and applies a conservative decision rule.  It deliberately
returns ``inconclusive`` when the record is too short or the instrument floor
cannot resolve the required voltage stability.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import NormalDist

import numpy as np

from iqs.constants import e_C, hbar, m_He


@dataclass(frozen=True)
class PhaseStabilityConfig:
    """Physical and statistical settings for a phase-stability gate."""

    energy_eV: float = 30_000.0
    path_length_m: float = 0.01
    recalibration_interval_s: float = 60.0
    phase_budget_rad: float = 0.05
    mass_kg: float = m_He
    charge_e: float = 1.0
    min_cycles: int = 20
    instrument_noise_rms_v: float | None = None
    instrument_peak_excursion_v: float | None = None
    confidence_sigma: float = 3.0

    def __post_init__(self):
        positive = (
            self.energy_eV,
            self.path_length_m,
            self.recalibration_interval_s,
            self.phase_budget_rad,
            self.mass_kg,
            self.charge_e,
            self.confidence_sigma,
        )
        if not all(np.isfinite(value) and value > 0 for value in positive):
            raise ValueError("physical settings must be finite and positive")
        if self.min_cycles < 2:
            raise ValueError("min_cycles must be at least two")
        if (self.instrument_noise_rms_v is not None
                and (not np.isfinite(self.instrument_noise_rms_v)
                     or self.instrument_noise_rms_v < 0)):
            raise ValueError("instrument noise must be finite and non-negative")
        if (self.instrument_peak_excursion_v is not None
                and (not np.isfinite(self.instrument_peak_excursion_v)
                     or self.instrument_peak_excursion_v < 0)):
            raise ValueError(
                "instrument peak excursion must be finite and non-negative"
            )


@dataclass(frozen=True)
class PhaseStabilityResult:
    """Summary returned by :func:`analyze_phase_stability`."""

    decision: str
    reasons: tuple[str, ...]
    duration_s: float
    sample_count: int
    complete_cycles: int
    velocity_m_s: float
    flight_time_s: float
    phase_sensitivity_rad_per_v: float
    voltage_budget_v: float
    instrument_noise_rms_v: float | None
    instrument_peak_excursion_v: float | None
    noise_excursion_bound_v: float | None
    linear_drift_v_s: float
    differential_rms_v: float
    p95_peak_excursion_v: float
    worst_peak_excursion_v: float
    p95_rms_excursion_v: float
    p95_peak_phase_rad: float
    worst_peak_phase_rad: float
    allan_tau_s: np.ndarray
    allan_deviation_v: np.ndarray

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["reasons"] = list(self.reasons)
        payload["allan_tau_s"] = self.allan_tau_s.tolist()
        payload["allan_deviation_v"] = self.allan_deviation_v.tolist()
        return payload


def particle_velocity(energy_eV: float, mass_kg: float = m_He) -> float:
    """Nonrelativistic particle velocity at kinetic energy ``energy_eV``."""
    energy_eV = float(energy_eV)
    mass_kg = float(mass_kg)
    if not np.isfinite(energy_eV) or energy_eV <= 0:
        raise ValueError("energy_eV must be finite and positive")
    if not np.isfinite(mass_kg) or mass_kg <= 0:
        raise ValueError("mass_kg must be finite and positive")
    return float(np.sqrt(2.0 * energy_eV * e_C / mass_kg))


def phase_sensitivity_rad_per_v(config: PhaseStabilityConfig) -> float:
    """Phase drift per volt of effective differential-potential drift."""
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    return float(
        config.charge_e * e_C * config.path_length_m / (hbar * velocity)
    )


def voltage_budget(config: PhaseStabilityConfig) -> float:
    """Maximum differential-voltage excursion for the allocated phase budget."""
    return float(config.phase_budget_rad / phase_sensitivity_rad_per_v(config))


def _validated_series(time_s, differential_voltage_v):
    time_s = np.asarray(time_s, dtype=float).reshape(-1)
    voltage_v = np.asarray(differential_voltage_v, dtype=float).reshape(-1)
    if time_s.size != voltage_v.size or time_s.size < 3:
        raise ValueError("time and voltage must have the same length >= 3")
    if not np.all(np.isfinite(time_s)) or not np.all(np.isfinite(voltage_v)):
        raise ValueError("time and voltage must be finite")
    if np.any(np.diff(time_s) <= 0):
        raise ValueError("time must be strictly increasing")
    return time_s, voltage_v


def calibration_window_excursions(
    time_s,
    differential_voltage_v,
    interval_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Peak and RMS drift after calibration for overlapping complete cycles.

    Every sample that has a full ``interval_s`` of subsequent data is treated
    as a possible calibration instant. Static differential voltage cancels;
    excursions are measured relative to the value at that instant.
    """
    time_s, voltage_v = _validated_series(time_s, differential_voltage_v)
    interval_s = float(interval_s)
    if not np.isfinite(interval_s) or interval_s <= 0:
        raise ValueError("interval_s must be finite and positive")

    last_start = time_s[-1] - interval_s
    starts = np.flatnonzero(time_s <= last_start)
    peaks = np.empty(starts.size, dtype=float)
    rms = np.empty(starts.size, dtype=float)
    for output_index, start in enumerate(starts):
        stop = np.searchsorted(
            time_s, time_s[start] + interval_s, side="right")
        delta = voltage_v[start:stop] - voltage_v[start]
        peaks[output_index] = np.max(np.abs(delta))
        rms[output_index] = np.sqrt(np.mean(delta**2))
    if peaks.size == 0:
        raise ValueError("record is shorter than one recalibration interval")
    return peaks, rms


def overlapping_allan_deviation(
    time_s,
    values,
    averaging_times_s=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return overlapping Allan deviation on an interpolated uniform grid."""
    time_s, values = _validated_series(time_s, values)
    dt = float(np.median(np.diff(time_s)))
    uniform_time = np.arange(time_s[0], time_s[-1] + 0.5 * dt, dt)
    uniform_values = np.interp(uniform_time, time_s, values)
    count = uniform_values.size

    if averaging_times_s is None:
        max_m = max(1, count // 4)
        m_values = np.unique(np.geomspace(1, max_m, 24).astype(int))
    else:
        requested = np.asarray(averaging_times_s, dtype=float).reshape(-1)
        if (requested.size == 0 or not np.all(np.isfinite(requested))
                or np.any(requested <= 0)):
            raise ValueError("averaging times must be finite and positive")
        m_values = np.unique(np.maximum(1, np.rint(requested / dt).astype(int)))

    taus = []
    deviations = []
    cumulative = np.concatenate(([0.0], np.cumsum(uniform_values)))
    for m in m_values:
        if 2 * m >= count:
            continue
        means = (cumulative[m:] - cumulative[:-m]) / m
        adjacent_difference = means[m:] - means[:-m]
        taus.append(m * dt)
        deviations.append(np.sqrt(0.5 * np.mean(adjacent_difference**2)))
    return np.asarray(taus), np.asarray(deviations)


def analyze_phase_stability(
    time_s,
    differential_voltage_v,
    config: PhaseStabilityConfig,
) -> PhaseStabilityResult:
    """Analyze a measurement and return a conservative architecture verdict."""
    time_s, voltage_v = _validated_series(time_s, differential_voltage_v)
    peaks, rms = calibration_window_excursions(
        time_s, voltage_v, config.recalibration_interval_s)
    allan_tau, allan_deviation = overlapping_allan_deviation(time_s, voltage_v)

    duration = float(time_s[-1] - time_s[0])
    complete_cycles = int(np.floor(duration / config.recalibration_interval_s))
    velocity = particle_velocity(config.energy_eV, config.mass_kg)
    flight_time = config.path_length_m / velocity
    sensitivity = phase_sensitivity_rad_per_v(config)
    budget_v = config.phase_budget_rad / sensitivity
    p95_peak = float(np.percentile(peaks, 95))
    worst_peak = float(np.max(peaks))
    p95_rms = float(np.percentile(rms, 95))
    drift_rate = float(np.polyfit(time_s - time_s[0], voltage_v, 1)[0])

    reasons = []
    decision = "inconclusive"
    noise_bound = config.instrument_peak_excursion_v
    if noise_bound is None and config.instrument_noise_rms_v is not None:
        # Bound the largest absolute difference across all samples in a
        # calibration window, not just one sample pair. The configured sigma
        # defines the family-wise Gaussian tail probability; Bonferroni then
        # allocates it across comparisons. A measured reference-window peak
        # is preferred because real instruments often contain non-Gaussian
        # and correlated noise.
        dt = float(np.median(np.diff(time_s)))
        comparisons = max(1, int(np.ceil(config.recalibration_interval_s / dt)))
        normal = NormalDist()
        one_sided_tail = 1.0 - normal.cdf(config.confidence_sigma)
        per_comparison_tail = one_sided_tail / comparisons
        family_z = normal.inv_cdf(1.0 - per_comparison_tail)
        noise_bound = (
            family_z * np.sqrt(2.0) * config.instrument_noise_rms_v
        )

    if complete_cycles < config.min_cycles:
        reasons.append(
            f"record contains {complete_cycles} complete cycles; "
            f"at least {config.min_cycles} are required"
        )
    elif noise_bound is None:
        reasons.append("instrument noise floor or reference excursion was not supplied")
    elif max(0.0, p95_peak - noise_bound) > budget_v:
        decision = "falsified"
        reasons.append(
            "95th-percentile calibration-window excursion exceeds the "
            "voltage budget after a conservative instrument-noise bound"
        )
    elif noise_bound > budget_v / 3.0:
        reasons.append(
            "instrument excursion bound exceeds one third of the voltage budget"
        )
    elif worst_peak + noise_bound <= budget_v:
        decision = "not_falsified"
        reasons.append(
            "every complete calibration window remains inside the voltage "
            "budget including the instrument-noise bound"
        )
    else:
        reasons.append(
            "measurement overlaps the decision boundary after accounting "
            "for the instrument-noise bound"
        )

    return PhaseStabilityResult(
        decision=decision,
        reasons=tuple(reasons),
        duration_s=duration,
        sample_count=int(time_s.size),
        complete_cycles=complete_cycles,
        velocity_m_s=velocity,
        flight_time_s=flight_time,
        phase_sensitivity_rad_per_v=sensitivity,
        voltage_budget_v=budget_v,
        instrument_noise_rms_v=config.instrument_noise_rms_v,
        instrument_peak_excursion_v=config.instrument_peak_excursion_v,
        noise_excursion_bound_v=(None if noise_bound is None else float(noise_bound)),
        linear_drift_v_s=drift_rate,
        differential_rms_v=float(np.std(voltage_v)),
        p95_peak_excursion_v=p95_peak,
        worst_peak_excursion_v=worst_peak,
        p95_rms_excursion_v=p95_rms,
        p95_peak_phase_rad=p95_peak * sensitivity,
        worst_peak_phase_rad=worst_peak * sensitivity,
        allan_tau_s=allan_tau,
        allan_deviation_v=allan_deviation,
    )


def load_voltage_csv(path: str | Path):
    """Load ``time_s`` plus differential or two-path voltage columns.

    Supported schemas are:

    - ``time_s,differential_voltage_v``
    - ``time_s,path_a_voltage_v,path_b_voltage_v``
    """
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    if data.dtype.names is None:
        raise ValueError("CSV must contain a header row")
    names = set(data.dtype.names)
    if "time_s" not in names:
        raise ValueError("CSV requires a time_s column")
    if "differential_voltage_v" in names:
        differential = data["differential_voltage_v"]
    elif {"path_a_voltage_v", "path_b_voltage_v"}.issubset(names):
        differential = data["path_a_voltage_v"] - data["path_b_voltage_v"]
    else:
        raise ValueError(
            "CSV requires differential_voltage_v or both path voltage columns"
        )
    return _validated_series(data["time_s"], differential)


def synthetic_voltage_channels(
    config: PhaseStabilityConfig,
    mode: str = "pass",
    sample_period_s: float = 0.5,
    cycles: int = 30,
    seed: int = 42,
):
    """Generate two channels with large common mode and controlled drift.

    ``pass`` stays comfortably below the architecture budget. ``fail`` has a
    calibration-window drift of roughly three times the budget. The returned
    instrument noise is the RMS noise of the differential channel.
    """
    if mode not in {"pass", "fail"}:
        raise ValueError("mode must be 'pass' or 'fail'")
    if sample_period_s <= 0 or cycles < 2:
        raise ValueError("sample period must be positive and cycles >= 2")

    budget_v = voltage_budget(config)
    duration = cycles * config.recalibration_interval_s
    time_s = np.arange(0.0, duration + 0.5 * sample_period_s, sample_period_s)
    rng = np.random.default_rng(seed)

    drift_per_cycle = (0.20 if mode == "pass" else 3.0) * budget_v
    drift_rate = drift_per_cycle / config.recalibration_interval_s
    differential_true = drift_rate * time_s
    differential_true += 0.03 * budget_v * np.sin(
        2 * np.pi * time_s / (3.7 * config.recalibration_interval_s)
    )

    differential_noise_rms = 0.02 * budget_v
    channel_noise_rms = differential_noise_rms / np.sqrt(2.0)
    common_mode = (
        1e-3
        + 2e-5 * np.sin(2 * np.pi * time_s / 17.0)
        + 1e-7 * time_s
    )
    path_a = (
        common_mode + 0.5 * differential_true
        + rng.normal(0.0, channel_noise_rms, time_s.size)
    )
    path_b = (
        common_mode - 0.5 * differential_true
        + rng.normal(0.0, channel_noise_rms, time_s.size)
    )
    return time_s, path_a, path_b, differential_noise_rms


__all__ = [
    "PhaseStabilityConfig",
    "PhaseStabilityResult",
    "analyze_phase_stability",
    "calibration_window_excursions",
    "load_voltage_csv",
    "overlapping_allan_deviation",
    "particle_velocity",
    "phase_sensitivity_rad_per_v",
    "synthetic_voltage_channels",
    "voltage_budget",
]
