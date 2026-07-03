"""
Density-to-lattice state mapper for the diamond-network caging stage.

Converts a 2D deposition density image (on the simulation grid) into an
initial quantum state on the diamond lattice by:

  1. Downsampling the density to lattice dimensions.
  2. Smoothing with a Gaussian filter.
  3. Finding local intensity maxima.
  4. Loading the top-N peaks onto A-sites (hub sites).
  5. Normalising the resulting state vector.
"""

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter, zoom


class DensityPeakMapper:
    """
    Map a 2D deposition density to an initial state on a DiamondNetwork.

    Parameters
    ----------
    N_grid : int
        Side length of the input density grid (pixels).
    n_peaks : int
        Maximum number of peaks to load onto A-sites.

    Attributes
    ----------
    last_n_peaks : int
        Actual number of peaks found and loaded in the most recent call
        to :meth:`map_to_a_sites`.  Useful for recording in result dicts.
    """

    def __init__(self, N_grid: int = 256, n_peaks: int = 16):
        self.N_grid      = N_grid
        self.n_peaks     = n_peaks
        self.last_n_peaks = 0

    def map_to_a_sites(self, density_2d: np.ndarray,
                       lattice) -> np.ndarray:
        """
        Build a normalised initial state from *density_2d*.

        Parameters
        ----------
        density_2d : np.ndarray, shape (N_grid, N_grid)
            2D intensity / probability-density image.
        lattice : DiamondNetwork
            Target lattice.  Must expose ``Lx``, ``Ly``, ``dim``,
            and ``_site(ix, iy, t)``.

        Returns
        -------
        psi0 : np.ndarray, shape (lattice.dim,), dtype complex
            Normalised initial state with amplitude loaded onto A-sites.
        """
        Lx, Ly = lattice.Lx, lattice.Ly
        dim    = lattice.dim

        # ── Downsample density to lattice dimensions ──────────────────
        scale_x   = Lx / self.N_grid
        scale_y   = Ly / self.N_grid
        density_ds = np.maximum(
            zoom(density_2d, (scale_x, scale_y), order=1), 0)
        d_smooth  = gaussian_filter(density_ds, sigma=1)

        # ── Find local maxima ─────────────────────────────────────────
        local_max  = maximum_filter(d_smooth, size=3) == d_smooth
        peak_mask  = local_max & (d_smooth > 0.1 * d_smooth.max())
        peak_idx   = np.argwhere(peak_mask)
        peak_vals  = d_smooth[peak_mask]

        # Take the top n_peaks by intensity
        order_     = np.argsort(peak_vals)[::-1]
        peak_idx   = peak_idx[order_[:self.n_peaks]]
        self.last_n_peaks = len(peak_idx)

        # ── Load peaks onto A-sites ───────────────────────────────────
        psi_np = np.zeros(dim, dtype=complex)
        for px, py in peak_idx:
            ix_ = min(px, Lx - 1)
            iy_ = min(py, Ly - 1)
            psi_np[lattice._site(ix_, iy_, 0)] += np.sqrt(d_smooth[px, py])

        # ── Normalise ─────────────────────────────────────────────────
        norm = np.linalg.norm(psi_np)
        if norm < 1e-12:
            # Fallback: uniform A-site loading
            psi_np[::5] = 1.0 / np.sqrt(Lx * Ly)
        else:
            psi_np /= norm

        return psi_np
