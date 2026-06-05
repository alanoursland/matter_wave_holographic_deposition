"""
Wavefunction normalization helpers.

Normalization convention: integral of |psi|^2 dA = 1
on a uniform 2D grid with pixel area dx^2.
"""

import numpy as np
import torch


def normalize_np(psi: np.ndarray, dx: float) -> np.ndarray:
    """
    Normalize *psi* so that sum(|psi|^2) * dx^2 = 1.

    Parameters
    ----------
    psi : np.ndarray, dtype complex
        Wavefunction on an N×N grid.
    dx : float
        Grid spacing [m].

    Returns
    -------
    np.ndarray
        Normalized copy of *psi*.
    """
    return psi / np.sqrt(np.sum(np.abs(psi) ** 2) * dx ** 2)


def normalize_torch(psi_t: torch.Tensor, dx: float) -> torch.Tensor:
    """
    Normalize *psi_t* so that sum(|psi|^2) * dx^2 = 1.

    Parameters
    ----------
    psi_t : torch.Tensor, dtype complex
        Wavefunction on an N×N grid.
    dx : float
        Grid spacing [m].

    Returns
    -------
    torch.Tensor
        Normalized tensor (same device and dtype as *psi_t*).
    """
    return psi_t / torch.sqrt(torch.sum(torch.abs(psi_t) ** 2) * dx ** 2)
