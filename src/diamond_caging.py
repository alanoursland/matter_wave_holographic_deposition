import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter, zoom
from typing import Tuple, Dict, Any

class DiamondNetworkSimulator:
    """
    Simulator for 2D Diamond Network dynamics.
    Replaces the Lieb lattice to provide true, flux-dependent Aharonov-Bohm caging.
    The unit cell has 5 sites: A (hub), B1/B2 (horizontal diamond), C1/C2 (vertical diamond).
    """
    
    def __init__(self, Lx: int, Ly: int, N_grid: int = 256):
        self.Lx = Lx
        self.Ly = Ly
        self.N_grid = N_grid
        self.dim = 5 * Lx * Ly

    def _site(self, ix: int, iy: int, t: int) -> int:
        # t=0: A, t=1: B1, t=2: B2, t=3: C1, t=4: C2
        return 5 * (ix * self.Ly + iy) + t

    def build_hamiltonian(self, phi: float = np.pi, J: float = 1.0, 
                          disorder_W: float = 0.0, disorder_frac: float = 1.0) -> np.ndarray:
        """
        Builds the Hamiltonian for the 2D Diamond network.
        A flux of Φ is threaded through each diamond plaquette.
        """
        H = np.zeros((self.dim, self.dim), dtype=complex)
        
        # Split the flux symmetrically across the upper/lower paths
        phase = np.exp(1j * phi / 2)
        phase_conj = np.exp(-1j * phi / 2)

        for ix in range(self.Lx):
            for iy in range(self.Ly):
                a = self._site(ix, iy, 0)
                b1 = self._site(ix, iy, 1) # Horizontal top path
                b2 = self._site(ix, iy, 2) # Horizontal bottom path
                c1 = self._site(ix, iy, 3) # Vertical right path
                c2 = self._site(ix, iy, 4) # Vertical left path

                a_right = self._site((ix + 1) % self.Lx, iy, 0)
                a_up = self._site(ix, (iy + 1) % self.Ly, 0)

                # Horizontal Diamond (A -> B1/B2 -> A_right)
                H[a, b1] = J
                H[b1, a] = J
                H[b1, a_right] = J * phase
                H[a_right, b1] = J * np.conj(phase)

                H[a, b2] = J
                H[b2, a] = J
                H[b2, a_right] = J * phase_conj
                H[a_right, b2] = J * np.conj(phase_conj)

                # Vertical Diamond (A -> C1/C2 -> A_up)
                H[a, c1] = J
                H[c1, a] = J
                H[c1, a_up] = J * phase
                H[a_up, c1] = J * np.conj(phase)

                H[a, c2] = J
                H[c2, a] = J
                H[c2, a_up] = J * phase_conj
                H[a_up, c2] = J * np.conj(phase_conj)

        # Apply random on-site disorder
        if disorder_W > 0:
            n_disorder = int(disorder_frac * self.dim)
            disorder_sites = np.random.choice(self.dim, n_disorder, replace=False)
            disorder_vals = np.random.uniform(-disorder_W / 2, disorder_W / 2, n_disorder)
            for s, dv in zip(disorder_sites, disorder_vals):
                H[s, s] += dv

        return H

    def validate_spectrum(self, phi: float = np.pi, J: float = 1.0, 
                          abort_on_fail: bool = True) -> Tuple[bool, np.ndarray, np.ndarray]:
        """
        Validates the emergence of AB caging. At Φ=π, all bands in the 
        Diamond network become perfectly flat, yielding exactly 3 unique eigenvalues.
        """
        H = self.build_hamiltonian(phi=phi, J=J, disorder_W=0.0)
        evals = np.linalg.eigvalsh(H)
        
        bands = np.unique(np.round(evals, 2))
        n_uniq = len(bands)
        
        # True AB caging yields 3 unique macroscopic eigenvalues
        ok = n_uniq <= 3 

        print(f"\n  DIAMOND SPECTRUM CHECK  ({self.Lx}×{self.Ly}, Φ={phi/np.pi:.1f}π)")
        print(f"  Unique eigenvalues: {n_uniq}  ({'OK (AB Caging)' if ok else 'FAIL'})")
            
        if not ok and abort_on_fail:
            raise RuntimeError(f"Diamond AB caging check failed ({n_uniq} unique eigenvalues, expected ≤3 at Φ=π).")
            
        return ok, H, evals

    def evolve_caging(self, density_2d: np.ndarray, phi_cage: float = np.pi, 
                      J_hop: float = 1.0, T_evolve: float = 40.0, dt: float = 0.05, 
                      n_peaks: int = 16, disorder_W: float = 0.0, 
                      disorder_frac: float = 1.0, verbose: bool = True) -> Dict[str, Any]:
        """
        Executes time evolution, correctly loading the peaks onto the A (hub) sites.
        """
        if verbose:
            print(f"\n  STAGE 4: 2D Diamond Caging  [Φ={phi_cage/np.pi:.2f}π, "
                  f"{self.Lx}×{self.Ly}" + (f", W={disorder_W:.2f}J" if disorder_W > 0 else "") + "]")

        try:
            flat_ok, H, evals_H = self.validate_spectrum(phi=phi_cage, J=J_hop, abort_on_fail=False)
        except RuntimeError:
            flat_ok = False
            H = self.build_hamiltonian(phi_cage, J_hop, disorder_W, disorder_frac)

        scale_x = self.Lx / self.N_grid
        scale_y = self.Ly / self.N_grid
        density_ds = np.maximum(zoom(density_2d, (scale_x, scale_y), order=1), 0)
        d_smooth = gaussian_filter(density_ds, sigma=1)
        
        local_max = maximum_filter(d_smooth, size=3) == d_smooth
        peak_mask = local_max & (d_smooth > 0.1 * d_smooth.max())
        peak_idx = np.argwhere(peak_mask)
        peak_vals = d_smooth[peak_mask]
        
        order_ = np.argsort(peak_vals)[::-1]
        peak_idx = peak_idx[order_[:n_peaks]]

        psi_diamond = np.zeros(self.dim, dtype=complex)
        for px, py in peak_idx:
            ix_ = min(px, self.Lx - 1)
            iy_ = min(py, self.Ly - 1)
            # Load directly onto the A-sites (t=0)
            psi_diamond[self._site(ix_, iy_, 0)] += np.sqrt(d_smooth[px, py])

        norm = np.linalg.norm(psi_diamond)
        if norm < 1e-12:
            psi_diamond[::5] = 1.0 / np.sqrt(self.Lx * self.Ly)
        else:
            psi_diamond /= norm

        psi0_diamond = psi_diamond.copy()
        
        evals_H, evecs_H = np.linalg.eigh(H)
        times = np.arange(0, T_evolve, dt)
        fid_h, spr_h = [], []

        # Map sites to approximate spatial coordinates for spread calculation
        site_x = np.array([ix + (0.5 if t in [1,2] else 0.0) for ix in range(self.Lx) for iy in range(self.Ly) for t in range(5)])
        site_y = np.array([iy + (0.5 if t in [3,4] else 0.0) for ix in range(self.Lx) for iy in range(self.Ly) for t in range(5)])

        coeff = evecs_H.conj().T @ psi_diamond

        for t in times:
            psi_t = evecs_H @ (np.exp(-1j * evals_H * t) * coeff)
            probs = np.abs(psi_t)**2
            
            fid_h.append(np.abs(np.dot(psi0_diamond.conj(), psi_t))**2)
            
            mu_x = np.dot(site_x, probs)
            mu_y = np.dot(site_y, probs)
            spr_h.append(np.sqrt(np.dot((site_x - mu_x)**2 + (site_y - mu_y)**2, probs)))

        return {
            'times': times,
            'fidelity': np.array(fid_h),
            'spread': np.array(spr_h),
            'phi': phi_cage,
            'n_peaks': len(peak_idx),
            'flat_ok': flat_ok,
        }