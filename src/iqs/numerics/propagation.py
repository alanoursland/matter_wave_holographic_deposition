"""
Angular-spectrum free-space propagator for matter waves and scalar fields.

Reference implementation used by both the patterned-substrate pipeline
(sim_v9 / PatternedSubstratePipeline) and the holographic-caging pipeline
(inverse_holography / HolographicCagingPipeline).
"""

import numpy as np
import torch

from iqs.numerics.device import get_device


class AngularSpectrumPropagator:
    """
    Precomputed angular-spectrum propagator for a fixed geometry.

    Transfer function:
        H(kx, ky) = exp(±i kz z),   kz = sqrt(k0² - kx² - ky²)
    Evanescent components (kz² ≤ 0) are zeroed.

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

    Examples
    --------
    >>> prop = AngularSpectrumPropagator(N=256, L=400e-9, k0=k0, z=z)
    >>> psi_out  = prop.forward(psi_in_tensor)
    >>> psi_back = prop.backward(psi_out)
    """

    def __init__(self, N: int, L: float, k0: float, z: float,
                 device: torch.device | None = None):
        if device is None:
            device = get_device()

        self.N      = N
        self.L      = L
        self.dx     = L / N
        self.k0     = k0
        self.z      = z
        self.device = device

        # Spatial-frequency axes
        kx = torch.fft.fftfreq(N, self.dx, device=device,
                                dtype=torch.float64) * 2 * np.pi
        ky = torch.fft.fftfreq(N, self.dx, device=device,
                                dtype=torch.float64) * 2 * np.pi
        KX, KY = torch.meshgrid(kx, ky, indexing='ij')

        # Axial wavenumber — evanescent modes are zeroed
        kz_sq = k0 ** 2 - KX ** 2 - KY ** 2
        valid = kz_sq > 0
        kz = torch.where(valid,
                         torch.sqrt(torch.clamp(kz_sq, min=0)),
                         torch.zeros_like(kz_sq))

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
        return torch.fft.ifft2(torch.fft.fft2(psi_t) * self._H_fwd)

    def backward(self, psi_t: torch.Tensor) -> torch.Tensor:
        """
        Propagate *psi_t* backward by z (conjugate transfer function).

        Parameters
        ----------
        psi_t : torch.Tensor, dtype complex128, shape (N, N)

        Returns
        -------
        torch.Tensor, dtype complex128, shape (N, N)
        """
        return torch.fft.ifft2(torch.fft.fft2(psi_t) * self._H_bwd)
