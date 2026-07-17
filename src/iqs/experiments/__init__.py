"""Reusable experiment and sweep entry points."""

from .ablation import ablation_study
from .decoherence_sweep import decoherence_sweep
from .disorder_robustness import disorder_robustness
from .dose_fidelity import dose_fidelity_curve, dose_to_ssim
from .multi_species import multi_species_deposition
from .column_aberration import (
    aberrated_intensity,
    transverse_k_squared,
    wavelength_from_energy,
)
from .electrostatic_column import (
    ElectrostaticRayTracer,
    ImagePlane,
    RayTrace,
    TransferMatrixTrace,
    aperture_clearances,
    image_planes,
)
from .quantum_dephasing import (
    array_factor_intensity,
    coherence_matrix_from_pairwise_phase_rms,
    dephasing_metrics,
    observable_phase_covariance,
)
from .phase_stability import (
    PhaseStabilityConfig,
    PhaseStabilityResult,
    analyze_phase_stability,
    load_voltage_csv,
)
from .thermal_phase_noise import (
    ThermalPhaseNoiseConfig,
    ThermalPhaseNoiseResult,
    analyze_thermal_phase_noise,
    multi_electrode_phase_covariance,
    pairwise_differential_phase_rms,
    weighted_integrated_voltage_variance,
)
from .pressure_sweep import pressure_sweep
from .sideband_selectivity import sideband_selectivity_sweep
from .spectral_phase_noise import (
    SpectralPhaseNoiseConfig,
    SpectralPhaseNoiseResult,
    analyze_spectral_phase_noise,
    quantum_rc_voltage_csd,
)

__all__ = [
    "ablation_study",
    "decoherence_sweep",
    "disorder_robustness",
    "dose_fidelity_curve",
    "dose_to_ssim",
    "multi_species_deposition",
    "aberrated_intensity",
    "transverse_k_squared",
    "wavelength_from_energy",
    "ElectrostaticRayTracer",
    "ImagePlane",
    "RayTrace",
    "TransferMatrixTrace",
    "aperture_clearances",
    "image_planes",
    "array_factor_intensity",
    "coherence_matrix_from_pairwise_phase_rms",
    "dephasing_metrics",
    "observable_phase_covariance",
    "PhaseStabilityConfig",
    "PhaseStabilityResult",
    "analyze_phase_stability",
    "load_voltage_csv",
    "ThermalPhaseNoiseConfig",
    "ThermalPhaseNoiseResult",
    "analyze_thermal_phase_noise",
    "multi_electrode_phase_covariance",
    "pairwise_differential_phase_rms",
    "weighted_integrated_voltage_variance",
    "pressure_sweep",
    "sideband_selectivity_sweep",
    "SpectralPhaseNoiseConfig",
    "SpectralPhaseNoiseResult",
    "analyze_spectral_phase_noise",
    "quantum_rc_voltage_csd",
]
