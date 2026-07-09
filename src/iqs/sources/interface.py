"""Source parameterization decoupled from the Kuramoto model (fable5 T12/M2).

The holographic pipeline never needs the Kuramoto abstraction: everything
downstream consumes a small set of measurable beam properties —
monochromaticity Δλ/λ, transverse coherence length ξ⊥, beam current, and
the transverse phase-noise RMS σ_θ (equivalently a coherence factor
r = exp(−σ_θ²/2)).  This module defines that interface plus two providers:

- `DirectSource` — you state the (measured or specified) parameters.
- `KuramotoPatentSource` — the AB-Kuramoto synchronization model from
  US Patent 9,502,202 B2, clearly labeled as one *model* of where the
  parameters come from rather than the foundation of the pipeline.
"""

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SourceParams:
    """The physical beam-source parameters downstream code consumes.

    Attributes
    ----------
    dlam_frac : float
        Monochromaticity Δλ/λ (= Δv/v = ΔE/E / 2).
    xi_perp : float
        Transverse coherence length of the source [m]; sets the
        correlation length of the phase-noise screens.
    current_A : float
        Beam current [A]; sets the ion arrival rate (T11) and the
        space-charge occupancy (T10).
    sigma_theta : float
        Transverse phase-noise RMS [rad].  For Gaussian phase disorder
        this is equivalent to a coherence factor r = exp(−σ_θ²/2).
    provider : str
        Which provider produced these numbers (for logging/plots).
    meta : dict
        Provider-specific extras (e.g. the Kuramoto sync result).
    """
    dlam_frac: float
    xi_perp: float
    current_A: float
    sigma_theta: float
    provider: str = 'direct'
    meta: dict = field(default_factory=dict)

    @property
    def r_equivalent(self):
        """Coherence factor r = ⟨e^{iθ}⟩ = exp(−σ_θ²/2)."""
        return float(np.exp(-self.sigma_theta**2 / 2))

    @staticmethod
    def sigma_theta_from_r(r):
        """Inverse mapping σ_θ = √(−2 ln r) (fable5 E2/T6)."""
        r = float(np.clip(r, 1e-12, 1.0))
        return float(np.sqrt(max(0.0, -2.0 * np.log(r))))


class DirectSource:
    """Source parameterized directly by measurable quantities.

    Use this when the beam properties are known (spec sheet, measurement,
    or a design study sweeping them) — no synchronization model involved.
    """

    def __init__(self, dlam_frac, xi_perp, current_A, sigma_theta=0.0):
        self.dlam_frac = dlam_frac
        self.xi_perp = xi_perp
        self.current_A = current_A
        self.sigma_theta = sigma_theta

    def source_params(self, seed=None, verbose=False):
        """Return (SourceParams, sync_like) — sync_like mimics the shape
        of a Kuramoto sync result for downstream compatibility."""
        params = SourceParams(
            dlam_frac=self.dlam_frac, xi_perp=self.xi_perp,
            current_A=self.current_A, sigma_theta=self.sigma_theta,
            provider='direct',
        )
        sync_like = {'r_final': params.r_equivalent, 'mode': 'direct'}
        if verbose:
            print(f"    Source (direct): σ_θ = {self.sigma_theta:.4f} rad "
                  f"(r_eq = {params.r_equivalent:.4f}), "
                  f"Δλ/λ = {self.dlam_frac:.2%}, "
                  f"ξ⊥ = {self.xi_perp*1e9:.0f} nm, "
                  f"I = {self.current_A:.3e} A")
        return params, sync_like


class KuramotoPatentSource:
    """The patent's AB-Kuramoto synchronization model as a provider.

    Wraps a `CoherentMatterwaveBeam` (US 9,502,202 B2 model): runs the
    Kuramoto synchronization and maps its order parameter to the physical
    source parameters via σ_θ = √(−2 ln r) and Δλ/λ = dE_frac/2.  This is
    the most speculative element of the chain (fable5 M2: a common
    external A-field cannot synchronize relative phases) — kept as one
    clearly-labeled option, not the foundation.

    Note: imports the top-level `coherent_matterwave_beam` module lazily;
    the full move of the CMWB model into iqs/sources is reorg milestone 3.
    """

    def __init__(self, cmwb, xi_perp):
        self.cmwb = cmwb
        self.xi_perp = xi_perp

    def source_params(self, seed=None, verbose=False, mode='auto'):
        sync = self.cmwb.synchronize(seed=seed, verbose=verbose, mode=mode)
        params = SourceParams(
            dlam_frac=self.cmwb.dE_frac / 2,
            xi_perp=self.xi_perp,
            current_A=self.cmwb.beam_current_A,
            sigma_theta=SourceParams.sigma_theta_from_r(sync['r_final']),
            provider='kuramoto-patent',
            meta={'sync': sync},
        )
        return params, sync
