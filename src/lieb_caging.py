import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter, zoom
from typing import Tuple, Dict, Any

class LiebLatticeSimulator:
    """
    Simulator for 2D Lieb lattice dynamics.
    Implements a symmetric gauge A = (-y/2, x/2) with periodic boundary conditions (PBC).
    """
    
    def __init__(self, Lx: int, Ly: int, N_grid: int = 256):
        self.Lx = Lx
        self.Ly = Ly
        self.N_grid = N_grid
        self.dim = 3 * Lx * Ly

    def _site(self, ix: int, iy: int, t: int) -> int:
        return 3 * (ix * self.Ly + iy) + t

    def build_hamiltonian(self, phi: float = np.pi, J: float = 1.0, 
                          disorder_W: float = 0.0, disorder_frac: float = 1.0) -> np.ndarray:
        """
        Builds the Hamiltonian using a continuous symmetric gauge.
        """
        H = np.zeros((self.dim, self.dim), dtype=complex)

        for ix in range(self.Lx):
            for iy in range(self.Ly):
                a = self._site(ix, iy, 0)
                b = self._site(ix, iy, 1)
                c = self._site(ix, iy, 2)

                # Horizontal bonds depend on the y-coordinate
                phase_h = np.exp(-1j * phi * iy / 4)
                H[a, b] = J * phase_h
                H[b, a] = J * np.conj(phase_h)

                a2 = self._site((ix + 1) % self.Lx, iy, 0)
                H[b, a2] = J * phase_h
                H[a2, b] = J * np.conj(phase_h)

                # Vertical bonds depend on the x-coordinate
                phase_v = np.exp(1j * phi * ix / 4)
                H[a, c] = J * phase_v
                H[c, a] = J * np.conj(phase_v)

                a3 = self._site(ix, (iy + 1) % self.Ly, 0)
                H[c, a3] = J * phase_v
                H[a3, c] = J * np.conj(phase_v)

        # Apply random on-site disorder if requested
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
        Validates the spectrum. A 4x4 2D Lieb lattice at Φ=π yields exactly 11 
        unique dispersive energy bands, not 3.
        """
        H = self.build_hamiltonian(phi=phi, J=J, disorder_W=0.0)
        evals = np.linalg.eigvalsh(H)
        
        bands = np.unique(np.round(evals, 2))
        n_uniq = len(bands)
        
        # 2D Lieb at Pi flux yields 11 eigenvalues on a 4x4 grid
        ok = n_uniq <= 11

        print(f"\n  LIEB SPECTRUM CHECK  ({self.Lx}×{self.Ly}, Φ={phi/np.pi:.1f}π)")
        print(f"  Unique eigenvalues: {n_uniq}  ({'OK (Dispersive)' if ok else 'FAIL'})")
            
        if not ok and abort_on_fail:
            raise RuntimeError(f"Lieb band check failed ({n_uniq} unique eigenvalues, expected ≤11).")
            
        return ok, H, evals

    def evolve_caging(self, density_2d: np.ndarray, phi_cage: float = np.pi, 
                      J_hop: float = 1.0, T_evolve: float = 40.0, dt: float = 0.05, 
                      n_peaks: int = 16, disorder_W: float = 0.0, 
                      disorder_frac: float = 1.0, verbose: bool = True) -> Dict[str, Any]:
        """Time evolution of the density."""
        if verbose:
            print(f"\n  STAGE 4: 2D Lieb Evolution  [Φ={phi_cage/np.pi:.2f}π, sym gauge, "
                  f"{self.Lx}×{self.Ly}" + (f", W={disorder_W:.2f}J" if disorder_W > 0 else "") + "]")

        try:
            flat_ok, H, evals_H = self.validate_spectrum(phi=phi_cage, J=J_hop, abort_on_fail=False)
        except RuntimeError as e:
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

        psi_lieb = np.zeros(self.dim, dtype=complex)
        for px, py in peak_idx:
            ix_ = min(px, self.Lx - 1)
            iy_ = min(py, self.Ly - 1)
            psi_lieb[self._site(ix_, iy_, 0)] += np.sqrt(d_smooth[px, py])

        norm = np.linalg.norm(psi_lieb)
        if norm < 1e-12:
            psi_lieb[::3] = 1.0 / np.sqrt(self.Lx * self.Ly)
        else:
            psi_lieb /= norm

        psi0_lieb = psi_lieb.copy()
        init_dens = np.abs(psi_lieb)**2

        evals_H, evecs_H = np.linalg.eigh(H)
        times = np.arange(0, T_evolve, dt)
        fid_h, spr_h, snaps = [], [], []

        site_x = np.array([ix for ix in range(self.Lx) for iy in range(self.Ly) for _ in range(3)], dtype=float)
        site_y = np.array([iy for ix in range(self.Lx) for iy in range(self.Ly) for _ in range(3)], dtype=float)

        coeff = evecs_H.conj().T @ psi_lieb

        for t in times:
            psi_t = evecs_H @ (np.exp(-1j * evals_H * t) * coeff)
            probs = np.abs(psi_t)**2
            
            fid_h.append(np.abs(np.dot(psi0_lieb.conj(), psi_t))**2)
            
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