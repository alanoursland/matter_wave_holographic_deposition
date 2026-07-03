"""
diamond_caging — backward-compatibility wrapper.

The implementation has moved to:
    iqs.lattices.diamond         (DiamondNetwork)
    iqs.lattices.density_mapping (DensityPeakMapper)

This module re-exports DiamondNetworkSimulator as a thin facade so that
existing code that does ``from diamond_caging import DiamondNetworkSimulator``
continues to work without modification.
"""

import numpy as np

from iqs.lattices.diamond import DiamondNetwork
from iqs.lattices.density_mapping import DensityPeakMapper

__all__ = ['DiamondNetworkSimulator']


class DiamondNetworkSimulator(DiamondNetwork):
    """
    Backward-compatible wrapper around DiamondNetwork + DensityPeakMapper.

    .. deprecated::
        Use ``iqs.lattices.diamond.DiamondNetwork`` and
        ``iqs.lattices.density_mapping.DensityPeakMapper`` directly.
    """

    def __init__(self, Lx: int, Ly: int, N_grid: int = 256):
        super().__init__(Lx, Ly)
        self.N_grid = N_grid

    def evolve_caging(self, density_2d: np.ndarray,
                      phi_cage: float = np.pi,
                      J_hop: float = 1.0,
                      T_evolve: float = 40.0,
                      dt: float = 0.05,
                      n_peaks: int = 16,
                      disorder_W: float = 0.0,
                      disorder_frac: float = 1.0,
                      verbose: bool = True) -> dict:
        """Map density to A-sites then run time evolution."""
        mapper = DensityPeakMapper(N_grid=self.N_grid, n_peaks=n_peaks)
        psi0   = mapper.map_to_a_sites(density_2d, self)
        return self.evolve(
            psi0,
            phi           = phi_cage,
            J             = J_hop,
            T             = T_evolve,
            dt            = dt,
            disorder_W    = disorder_W,
            disorder_frac = disorder_frac,
            n_peaks       = mapper.last_n_peaks,
            verbose       = verbose,
        )
