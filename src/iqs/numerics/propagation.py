"""
Angular-spectrum free-space propagator for matter waves and scalar fields.

Reference implementation used by both the patterned-substrate pipeline
(sim_v9 / PatternedSubstratePipeline) and the holographic-caging pipeline
(inverse_holography / HolographicCagingPipeline).
"""

import numpy as np
import torch
import torch.nn.functional as F

from iqs.numerics.device import get_device


class AngularSpectrumPropagator:
    """
    Precomputed angular-spectrum propagator for a fixed geometry.

    Transfer function:
        H(kx, ky) = exp(±i kz z),   kz = sqrt(k0² - kx² - ky²)
    Evanescent components (kz² ≤ 0) are zeroed.  Note this makes
    ``backward`` the inverse of ``forward`` only on the propagating
    subspace: power carried by evanescent modes is lost, not restored.

    The field is zero-padded by ``pad_factor`` before the FFT (fable5 M4):
    an unpadded FFT imposes periodic boundary conditions, so power
    diffracted out of the frame wraps around and re-enters — at z of a few
    aperture widths this contaminates the target-plane intensity and the
    diffraction-efficiency metric.  Padding lets that power leave the frame.

    In addition, the transfer function is band-limited per the
    Matsushima–Shimobaba criterion (``band_limit=True``): transverse
    frequencies whose walk-off exceeds the padded frame under-sample
    exp(i kz z) and alias back in; they are zeroed.  With this limit the
    result is converged with respect to ``pad_factor`` (verified: SSIM of
    a GD holography solve is pad-independent for pad ≥ 2), which unpadded
    or naively padded ASM is not.

    A norm gate (fable5 T4) monitors the fraction of input power that falls
    in the evanescent band and prints a warning once per instance when it
    exceeds ``evanescent_warn_frac``.  A large fraction almost always means
    the input field is unphysical (e.g. a transverse phase ramp at ~k0,
    fable5 E1), not that the physics intends near-field structure.

    Parameters
    ----------
    N : int
        Grid size (N × N pixels).
    L : float
        Physical side length of the square grid [m].
    k0 : float
        Wavenumber of the matter wave [rad/m].
    z : float
        Propagation distance [m].
    device : torch.device, optional
        Computation device.  Defaults to GPU if available, else CPU.
    pad_factor : int
        Zero-padding factor for the FFT grid (1 = legacy unpadded
        behaviour; default 4: field embedded in a (4N)² frame, which
        keeps residual band-edge wraparound small — see band_limit).
    evanescent_warn_frac : float
        Warn (once) when more than this fraction of the input power is
        evanescent-zeroed at the input plane.

    Examples
    --------
    >>> prop = AngularSpectrumPropagator(N=256, L=400e-9, k0=k0, z=z)
    >>> psi_out  = prop.forward(psi_in_tensor)
    >>> psi_back = prop.backward(psi_out)
    """

    def __init__(self, N: int, L: float, k0: float, z: float,
                 device: torch.device | None = None,
                 pad_factor: int = 4,
                 band_limit: bool = True,
                 evanescent_warn_frac: float = 0.02):
        if device is None:
            device = get_device()

        self.N      = N
        self.L      = L
        self.dx     = L / N
        self.k0     = k0
        self.z      = z
        self.device = device

        self.pad_factor = max(int(pad_factor), 1)
        self.N_pad      = self.pad_factor * N
        self._pad       = (self.N_pad - N) // 2

        self.evanescent_warn_frac = evanescent_warn_frac
        self.last_evanescent_frac = 0.0
        self._warned = False

        # Spatial-frequency axes on the padded grid (same dx, finer dk)
        kx = torch.fft.fftfreq(self.N_pad, self.dx, device=device,
                                dtype=torch.float64) * 2 * np.pi
        ky = torch.fft.fftfreq(self.N_pad, self.dx, device=device,
                                dtype=torch.float64) * 2 * np.pi
        KX, KY = torch.meshgrid(kx, ky, indexing='ij')

        # Axial wavenumber — evanescent modes are zeroed
        kz_sq = k0 ** 2 - KX ** 2 - KY ** 2
        valid = kz_sq > 0

        # Band-limited ASM (fable5 T3).  Two cutoffs, both required:
        #
        # 1. Matsushima–Shimobaba anti-aliasing limit: at large z the
        #    transfer-function phase kz*z is under-sampled at high
        #    transverse frequency, and the DFT wraps that power back into
        #    the frame:  f_M = 1 / (lambda sqrt((2 df z)^2 + 1)),
        #    df = 1/(N_pad dx).
        # 2. Geometric aperture limit, INDEPENDENT of pad_factor: for a
        #    source window and observation window both of size L at
        #    distance z, only angles with tan(theta) <= L/z can carry
        #    power from source to observation window.  Steeper components
        #    physically exit sideways and never return; keeping the cutoff
        #    pad-independent is what makes results converge as pad_factor
        #    grows (the pure Matsushima limit rises with pad and re-admits
        #    near-grazing modes that alias at any finite pad).
        if band_limit:
            lam0 = 2 * np.pi / k0
            df = 1.0 / (self.N_pad * self.dx)
            f_M = 1.0 / (lam0 * np.sqrt((2 * df * z) ** 2 + 1))
            tan_g = L / z
            k_geo = k0 * tan_g / np.sqrt(1.0 + tan_g ** 2)
            k_lim = min(2 * np.pi * f_M, k_geo)
            self.k_limit = k_lim
            valid = valid & (KX.abs() <= k_lim) & (KY.abs() <= k_lim)
        else:
            self.k_limit = k0

        kz = torch.where(valid,
                         torch.sqrt(torch.clamp(kz_sq, min=0)),
                         torch.zeros_like(kz_sq))

        self._valid = valid
        self._H_fwd = torch.where(
            valid,
            torch.exp(1j * kz * z),
            torch.zeros_like(kz, dtype=torch.complex128),
        )
        self._H_bwd = torch.where(
            valid,
            torch.exp(-1j * kz * z),
            torch.zeros_like(kz, dtype=torch.complex128),
        )

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _embed(self, psi_t: torch.Tensor) -> torch.Tensor:
        """Zero-pad (N, N) -> (N_pad, N_pad), field centred."""
        if self.pad_factor == 1:
            return psi_t
        p = self._pad
        return F.pad(psi_t, (p, p, p, p))

    def _crop(self, psi_t: torch.Tensor) -> torch.Tensor:
        """Crop (N_pad, N_pad) -> central (N, N)."""
        if self.pad_factor == 1:
            return psi_t
        p = self._pad
        return psi_t[p:p + self.N, p:p + self.N]

    def _check_evanescent(self, spec_t: torch.Tensor) -> None:
        """Norm gate: warn once if the input spectrum is heavily evanescent."""
        with torch.no_grad():
            power = torch.abs(spec_t) ** 2
            total = power.sum()
            if total > 0:
                frac = float(1.0 - power[self._valid].sum() / total)
            else:
                frac = 0.0
        self.last_evanescent_frac = frac
        if frac > self.evanescent_warn_frac and not self._warned:
            self._warned = True
            print(f"  [AngularSpectrumPropagator] WARNING: {frac:.1%} of the "
                  f"input power lies outside the propagating band "
                  f"(evanescent |k_t| > k0, or beyond the band-limited "
                  f"aperture k_t > k_limit) and will be zeroed.  Check the "
                  f"input field for unphysical transverse phase structure "
                  f"(fable5 E1/T4).")

    def _apply(self, psi_t: torch.Tensor, H: torch.Tensor) -> torch.Tensor:
        spec = torch.fft.fft2(self._embed(psi_t))
        self._check_evanescent(spec)
        return self._crop(torch.fft.ifft2(spec * H))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def forward(self, psi_t: torch.Tensor) -> torch.Tensor:
        """
        Propagate *psi_t* forward by z.

        Parameters
        ----------
        psi_t : torch.Tensor, dtype complex128, shape (N, N)

        Returns
        -------
        torch.Tensor, dtype complex128, shape (N, N)
        """
        return self._apply(psi_t, self._H_fwd)

    def backward(self, psi_t: torch.Tensor) -> torch.Tensor:
        """
        Propagate *psi_t* backward by z (conjugate transfer function).

        Inverse of :meth:`forward` on the propagating subspace only —
        evanescent content zeroed by ``forward`` is not recoverable.

        Parameters
        ----------
        psi_t : torch.Tensor, dtype complex128, shape (N, N)

        Returns
        -------
        torch.Tensor, dtype complex128, shape (N, N)
        """
        return self._apply(psi_t, self._H_bwd)
