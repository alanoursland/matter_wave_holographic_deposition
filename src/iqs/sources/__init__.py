"""Beam-source interface (fable5 T12).

Downstream code consumes only the physical source parameters
(Δλ/λ, ξ⊥, current, σ_θ) via `SourceParams`; where those numbers come
from is a provider's business.  `DirectSource` states them outright;
`KuramotoPatentSource` derives them from the patent's AB-Kuramoto model.
"""

from iqs.sources.interface import (
    SourceParams, DirectSource, KuramotoPatentSource,
)
from iqs.sources.coherence import gaussian_wavelength_samples, fractional_rms

__all__ = [
    'SourceParams', 'DirectSource', 'KuramotoPatentSource',
    'gaussian_wavelength_samples', 'fractional_rms',
    'CavityGeometry', 'CoherentMatterwaveBeam', 'SPECIES',
    'add_species', 'coherence_sigma_theta', 'sample_phase_noise',
]


def __getattr__(name):
    if name in {
        'CavityGeometry', 'CoherentMatterwaveBeam', 'SPECIES',
        'add_species', 'coherence_sigma_theta', 'sample_phase_noise',
    }:
        from . import coherent_matterwave

        return getattr(coherent_matterwave, name)
    raise AttributeError(name)
