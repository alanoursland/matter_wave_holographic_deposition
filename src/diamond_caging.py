import numpy as np
import torch
from scipy.ndimage import gaussian_filter, maximum_filter, zoom
from typing import Tuple, Dict, Any

# Select GPU if available, else CPU
_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

class DiamondNetworkSimulator:
    """
    Simulator for 2D Diamond Network dynamics.
    Replaces the Lieb lattice to provide true, flux-dependent Aharonov-Bohm caging.
    The unit cell has 5 sites: A (hub), B1/B2 (horizontal diamond), C1/C2 (vertical diamond).
    GPU-accelerated via PyTorch for diagonalization and time evolution.
    """

    def __init__(self, Lx: int, Ly: int, N_grid: int = 256):
        self.Lx = Lx
        self.Ly = Ly
        self.N_grid = N_grid
        self.dim = 5 * Lx * Ly

    def _site(self, ix: int, iy: int, t: int) -> int:
        return 5 * (ix * self.Ly + iy) + t

    def build_hamiltonian(self, phi: float = np.pi, J: float = 1.0,
                          disorder_W: float = 0.0, disorder_frac: float = 1.0) -> np.ndarray:
        """
        Builds the Hamiltonian for the 2D Diamond network using vectorized index arrays.
        Returns numpy ndarray for compatibility.
        """
        Lx, Ly = self.Lx, self.Ly
        dim = self.dim

        # Build all (ix, iy) pairs
        ix_all = np.arange(Lx)
        iy_all = np.arange(Ly)
        IX, IY = np.meshgrid(ix_all, iy_all, indexing='ij')  # (Lx, Ly)
        IX_f = IX.ravel()  # (Lx*Ly,)
        IY_f = IY.ravel()

        # Site indices for all cells — shape (Lx*Ly,)
        a  = 5 * (IX_f * Ly + IY_f)       # A sites
        b1 = a + 1                          # B1 (horizontal top)
        b2 = a + 2                          # B2 (horizontal bottom)
        c1 = a + 3                          # C1 (vertical right)
        c2 = a + 4                          # C2 (vertical left)

        # Neighbor A sites (PBC)
        a_right = 5 * (((IX_f + 1) % Lx) * Ly + IY_f)
        a_up    = 5 * (IX_f * Ly + ((IY_f + 1) % Ly))

        phase_p = J * np.exp(1j * phi / 2)
        phase_m = J * np.exp(-1j * phi / 2)

        # Build sparse-style row/col/val arrays, then assign all at once
        H = np.zeros((dim, dim), dtype=complex)

        # Horizontal diamond: A <-> B1 <-> A_right, A <-> B2 <-> A_right
        # A <-> B1 (no phase)
        H[a, b1] = J;  H[b1, a] = J
        # B1 <-> A_right (phase +)
        H[b1, a_right] = phase_p;  H[a_right, b1] = np.conj(phase_p)
        # A <-> B2 (no phase)
        H[a, b2] = J;  H[b2, a] = J
        # B2 <-> A_right (phase -)
        H[b2, a_right] = phase_m;  H[a_right, b2] = np.conj(phase_m)

        # Vertical diamond: A <-> C1 <-> A_up, A <-> C2 <-> A_up
        H[a, c1] = J;  H[c1, a] = J
        H[c1, a_up] = phase_p;  H[a_up, c1] = np.conj(phase_p)
        H[a, c2] = J;  H[c2, a] = J
        H[c2, a_up] = phase_m;  H[a_up, c2] = np.conj(phase_m)

        # On-site disorder
        if disorder_W > 0:
            n_disorder = int(disorder_frac * dim)
            sites = np.random.choice(dim, n_disorder, replace=False)
            vals = np.random.uniform(-disorder_W / 2, disorder_W / 2, n_disorder)
            H[sites, sites] += vals

        return H

    def validate_spectrum(self, phi: float = np.pi, J: float = 1.0,
                          abort_on_fail: bool = True) -> Tuple[bool, np.ndarray, np.ndarray]:
        """
        Validates AB caging. At phi=pi, exactly 3 unique eigenvalues expected.
        Uses GPU for diagonalization.
        """
        H = self.build_hamiltonian(phi=phi, J=J, disorder_W=0.0)
        H_t = torch.tensor(H, dtype=torch.complex128, device=_device)
        evals_t = torch.linalg.eigvalsh(H_t)
        evals = evals_t.cpu().numpy()

        bands = np.unique(np.round(evals, 2))
        n_uniq = len(bands)
        ok = n_uniq <= 3

        print(f"\n  DIAMOND SPECTRUM CHECK  ({self.Lx}x{self.Ly}, phi={phi/np.pi:.1f}pi)")
        print(f"  Unique eigenvalues: {n_uniq}  ({'OK (AB Caging)' if ok else 'FAIL'})")

        if not ok and abort_on_fail:
            raise RuntimeError(f"Diamond AB caging check failed ({n_uniq} unique eigenvalues, expected <=3 at phi=pi).")

        return ok, H, evals

    def evolve_caging(self, density_2d: np.ndarray, phi_cage: float = np.pi,
                      J_hop: float = 1.0, T_evolve: float = 40.0, dt: float = 0.05,
                      n_peaks: int = 16, disorder_W: float = 0.0,
                      disorder_frac: float = 1.0, verbose: bool = True) -> Dict[str, Any]:
        """
        Time evolution with GPU-accelerated diagonalization and vectorized time stepping.
        Peak finding stays on CPU (scipy). Returns numpy arrays.
        """
        if verbose:
            print(f"\n  STAGE 4: 2D Diamond Caging  [phi={phi_cage/np.pi:.2f}pi, "
                  f"{self.Lx}x{self.Ly}" + (f", W={disorder_W:.2f}J" if disorder_W > 0 else "") + "]")

        Lx, Ly = self.Lx, self.Ly
        dim = self.dim

        # Spectrum check (clean phi=pi only)
        flat_ok = True
        if abs(phi_cage - np.pi) < 0.01 and disorder_W == 0:
            flat_ok, _, _ = self.validate_spectrum(phi=phi_cage, J=J_hop, abort_on_fail=False)

        # Build Hamiltonian (numpy), move to GPU for diag
        H_np = self.build_hamiltonian(phi=phi_cage, J=J_hop,
                                       disorder_W=disorder_W,
                                       disorder_frac=disorder_frac)
        H_t = torch.tensor(H_np, dtype=torch.complex128, device=_device)

        # --- Peak finding on CPU (scipy) ---
        scale_x = Lx / self.N_grid
        scale_y = Ly / self.N_grid
        density_ds = np.maximum(zoom(density_2d, (scale_x, scale_y), order=1), 0)
        d_smooth = gaussian_filter(density_ds, sigma=1)

        local_max = maximum_filter(d_smooth, size=3) == d_smooth
        peak_mask = local_max & (d_smooth > 0.1 * d_smooth.max())
        peak_idx = np.argwhere(peak_mask)
        peak_vals = d_smooth[peak_mask]

        order_ = np.argsort(peak_vals)[::-1]
        peak_idx = peak_idx[order_[:n_peaks]]

        # Load peaks onto A-sites
        psi_np = np.zeros(dim, dtype=complex)
        for px, py in peak_idx:
            ix_ = min(px, Lx - 1)
            iy_ = min(py, Ly - 1)
            psi_np[self._site(ix_, iy_, 0)] += np.sqrt(d_smooth[px, py])

        norm = np.linalg.norm(psi_np)
        if norm < 1e-12:
            psi_np[::5] = 1.0 / np.sqrt(Lx * Ly)
        else:
            psi_np /= norm

        psi0_np = psi_np.copy()
        init_dens = np.abs(psi_np)**2

        # --- GPU diagonalization ---
        evals_t, evecs_t = torch.linalg.eigh(H_t)

        # --- Build site coordinate and mask arrays (vectorized) ---
        ix_all = np.arange(Lx)
        iy_all = np.arange(Ly)
        IX, IY = np.meshgrid(ix_all, iy_all, indexing='ij')
        # Expand to all 5 sublattices
        site_x_np = np.repeat(IX.ravel(), 5).astype(float)
        site_y_np = np.repeat(IY.ravel(), 5).astype(float)

        # Localization mask: sites within +/-1 cell of any peak
        local_mask_np = np.zeros(dim, dtype=bool)
        for px, py in peak_idx:
            ix_ = min(px, Lx - 1)
            iy_ = min(py, Ly - 1)
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    base = 5 * (((ix_ + dx) % Lx) * Ly + ((iy_ + dy) % Ly))
                    local_mask_np[base:base+5] = True

        # Move vectors to GPU
        psi_t = torch.tensor(psi_np, dtype=torch.complex128, device=_device)
        psi0_t = torch.tensor(psi0_np, dtype=torch.complex128, device=_device)
        site_x_t = torch.tensor(site_x_np, dtype=torch.float64, device=_device)
        site_y_t = torch.tensor(site_y_np, dtype=torch.float64, device=_device)
        local_mask_t = torch.tensor(local_mask_np, device=_device)

        # Expansion coefficients: coeff = V^H @ psi, shape (dim,)
        coeff_t = evecs_t.conj().T @ psi_t

        # --- Vectorized time evolution ---
        times_np = np.arange(0, T_evolve, dt)
        N_times = len(times_np)
        times_t = torch.tensor(times_np, dtype=torch.float64, device=_device)

        # Phase factors: shape (N_times, dim)
        phase_factors = torch.exp(-1j * evals_t.unsqueeze(0) * times_t.unsqueeze(1))
        # weighted coefficients: (N_times, dim)
        weighted = phase_factors * coeff_t.unsqueeze(0)
        # psi at all times: (N_times, dim) = weighted @ evecs^T
        psi_all = weighted @ evecs_t.T  # (N_times, dim)

        # Probabilities: (N_times, dim)
        probs_all = torch.abs(psi_all)**2

        # Fidelity: |<psi0|psi(t)>|^2
        fid_all = torch.abs(psi_all @ psi0_t.conj())**2  # (N_times,)

        # Spread: sqrt( sum((x-mu_x)^2 + (y-mu_y)^2) * p )
        mu_x = probs_all @ site_x_t  # (N_times,)
        mu_y = probs_all @ site_y_t
        dx2 = (site_x_t.unsqueeze(0) - mu_x.unsqueeze(1))**2
        dy2 = (site_y_t.unsqueeze(0) - mu_y.unsqueeze(1))**2
        spr_all = torch.sqrt((probs_all * (dx2 + dy2)).sum(dim=1))

        # Localization: sum of probs on masked sites
        loc_all = probs_all[:, local_mask_t].sum(dim=1)

        # Move results to CPU
        fid_h = fid_all.cpu().numpy()
        spr_h = spr_all.cpu().numpy()
        loc_h = loc_all.cpu().numpy()
        probs_cpu = probs_all.cpu().numpy()

        # Snapshots at specific times
        snap_t_vals = [0, T_evolve * 0.25, T_evolve * 0.5, T_evolve]
        snaps = []
        # Build A-site index array for snapshot extraction
        a_indices = np.array([self._site(ixx, iyy, 0)
                              for ixx in range(Lx) for iyy in range(Ly)])
        for st in snap_t_vals:
            tidx = np.argmin(np.abs(times_np - st))
            probs_t = probs_cpu[tidx]
            ad = probs_t[a_indices].reshape(Lx, Ly)
            snaps.append((times_np[tidx], ad.copy()))

        final_fid = fid_h[-1]
        final_loc = loc_h[-1]
        if verbose:
            print(f"    {len(peak_idx)} peaks loaded onto A-sites")
            print(f"    Fidelity:      {final_fid:.4f}  "
                  f"({'PRESERVED' if final_fid > 0.5 else 'decayed'})")
            print(f"    Localization:  {final_loc:.4f}  "
                  f"({'CAGED' if final_loc > 0.8 else 'spread'})")
            print(f"    Spread: {spr_h[0]:.2f} -> {spr_h[-1]:.2f}  "
                  f"(delta={spr_h[-1]-spr_h[0]:+.2f} cells)")

        # Init 2D density (A-sites only)
        init_2d = init_dens[a_indices].reshape(Lx, Ly)

        return {
            'times':         times_np,
            'fidelity':      fid_h,
            'spread':        spr_h,
            'localization':  loc_h,
            'init_2d':       init_2d,
            'snapshots':     snaps,
            'phi':           phi_cage,
            'n_peaks':       len(peak_idx),
            'flat_ok':       flat_ok,
        }
