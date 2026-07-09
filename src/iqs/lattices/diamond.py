"""
2D Diamond-network tight-binding model with flux-dependent Aharonov-Bohm caging.

Unit cell — 5 sites:
  A  : hub site (coordination 8: four diamonds — right/left/up/down —
       meet at each hub, contributing two bridge bonds each)
  B1 : horizontal diamond, upper bridge
  B2 : horizontal diamond, lower bridge
  C1 : vertical diamond, upper bridge
  C2 : vertical diamond, lower bridge

At Φ = π per diamond all five bands are flat (exact AB caging).
At Φ = 0 the bands are dispersive and the state delocalises.

GPU-accelerated via PyTorch for Hamiltonian diagonalisation and
vectorised time evolution.
"""

import numpy as np
import torch
from typing import Tuple, Dict, Any

from iqs.numerics.device import get_device

_device = get_device()


class DiamondNetwork:
    """
    2D diamond-network tight-binding model.

    Parameters
    ----------
    Lx, Ly : int
        Lattice dimensions in unit cells along x and y.
    """

    def __init__(self, Lx: int, Ly: int):
        self.Lx  = Lx
        self.Ly  = Ly
        self.dim = 5 * Lx * Ly

    # ------------------------------------------------------------------
    # Index helpers
    # ------------------------------------------------------------------

    def _site(self, ix: int, iy: int, t: int) -> int:
        """
        Linear index of sublattice *t* (0=A, 1=B1, 2=B2, 3=C1, 4=C2)
        at unit cell (ix, iy).
        """
        return 5 * (ix * self.Ly + iy) + t

    # ------------------------------------------------------------------
    # Hamiltonian
    # ------------------------------------------------------------------

    def build_hamiltonian(self, phi: float = np.pi, J: float = 1.0,
                          disorder_W: float = 0.0,
                          disorder_frac: float = 1.0,
                          disorder_seed: int | None = None) -> np.ndarray:
        """
        Build the tight-binding Hamiltonian with AB phase factors.

        At phi = π: all 5 bands are flat → complete caging.
        Exactly 3 unique eigenvalues: {−2√2, 0, +2√2} J.

        Parameters
        ----------
        phi : float
            AB flux per diamond [rad].  phi = π gives complete caging.
        J : float
            Hopping amplitude.
        disorder_W : float
            On-site disorder strength (uniform in [−W/2, W/2]).
        disorder_frac : float
            Fraction of sites that receive disorder.
        disorder_seed : int or None
            Seed for the disorder realisation.  None (default) draws from
            the global NumPy RNG (legacy behaviour — reproducible only if
            the caller seeds np.random itself).

        Returns
        -------
        H : np.ndarray, shape (dim, dim), dtype complex
        """
        Lx, Ly = self.Lx, self.Ly
        dim = self.dim

        ix_all = np.arange(Lx)
        iy_all = np.arange(Ly)
        IX, IY = np.meshgrid(ix_all, iy_all, indexing='ij')
        IX_f   = IX.ravel()
        IY_f   = IY.ravel()

        # Site indices per cell
        a  = 5 * (IX_f * Ly + IY_f)
        b1 = a + 1
        b2 = a + 2
        c1 = a + 3
        c2 = a + 4

        # Neighbour A sites (periodic boundary conditions)
        a_right = 5 * (((IX_f + 1) % Lx) * Ly + IY_f)
        a_up    = 5 * (IX_f * Ly + ((IY_f + 1) % Ly))

        phase_p = J * np.exp( 1j * phi / 2)
        phase_m = J * np.exp(-1j * phi / 2)

        H = np.zeros((dim, dim), dtype=complex)

        # Horizontal diamond: A ↔ B1 ↔ A_right, A ↔ B2 ↔ A_right
        H[a, b1] = J;  H[b1, a] = J
        H[b1, a_right] = phase_p;  H[a_right, b1] = np.conj(phase_p)
        H[a, b2] = J;  H[b2, a] = J
        H[b2, a_right] = phase_m;  H[a_right, b2] = np.conj(phase_m)

        # Vertical diamond: A ↔ C1 ↔ A_up, A ↔ C2 ↔ A_up
        H[a, c1] = J;  H[c1, a] = J
        H[c1, a_up] = phase_p;  H[a_up, c1] = np.conj(phase_p)
        H[a, c2] = J;  H[c2, a] = J
        H[c2, a_up] = phase_m;  H[a_up, c2] = np.conj(phase_m)

        # On-site disorder (seedable; global RNG when no seed, for
        # backward compatibility with callers that seed np.random)
        if disorder_W > 0:
            rng = (np.random.default_rng(disorder_seed)
                   if disorder_seed is not None else np.random)
            n_dis = int(disorder_frac * dim)
            sites = rng.choice(dim, n_dis, replace=False)
            vals  = rng.uniform(-disorder_W / 2, disorder_W / 2, n_dis)
            H[sites, sites] += vals

        return H

    # ------------------------------------------------------------------
    # Spectrum validation
    # ------------------------------------------------------------------

    def validate_spectrum(self, phi: float = np.pi, J: float = 1.0,
                          abort_on_fail: bool = True
                          ) -> Tuple[bool, np.ndarray, np.ndarray]:
        """
        Verify AB caging by checking that at phi = π the spectrum has
        at most 3 unique eigenvalues.

        Parameters
        ----------
        phi : float
            Flux per diamond [rad].
        J : float
            Hopping amplitude.
        abort_on_fail : bool
            Raise RuntimeError on failure if True.

        Returns
        -------
        ok : bool
        H : np.ndarray
        evals : np.ndarray
        """
        H     = self.build_hamiltonian(phi=phi, J=J, disorder_W=0.0)
        H_t   = torch.tensor(H, dtype=torch.complex128, device=_device)
        evals = torch.linalg.eigvalsh(H_t).cpu().numpy()

        bands  = np.unique(np.round(evals, 2))
        n_uniq = len(bands)
        ok     = n_uniq <= 3

        print(f"\n  DIAMOND SPECTRUM CHECK  ({self.Lx}x{self.Ly}, "
              f"phi={phi/np.pi:.1f}pi)")
        print(f"  Unique eigenvalues: {n_uniq}  "
              f"({'OK (AB Caging)' if ok else 'FAIL'})")

        if not ok and abort_on_fail:
            raise RuntimeError(
                f"Diamond AB caging check failed "
                f"({n_uniq} unique eigenvalues, expected ≤3 at phi=π).")

        return ok, H, evals

    # ------------------------------------------------------------------
    # Time evolution
    # ------------------------------------------------------------------

    def evolve(self, psi0: np.ndarray,
               phi: float = np.pi, J: float = 1.0,
               T: float = 40.0, dt: float = 0.05,
               disorder_W: float = 0.0, disorder_frac: float = 1.0,
               disorder_seed: int | None = None,
               n_peaks: int | None = None,
               verbose: bool = True) -> Dict[str, Any]:
        """
        Exact time evolution of an initial state on the diamond lattice.

        Uses GPU-accelerated diagonalisation and vectorised phase factors.

        Parameters
        ----------
        psi0 : np.ndarray, shape (dim,), dtype complex
            Initial state vector.  Should be normalised; if norm < 1e-12
            a uniform A-site state is substituted automatically.
        phi : float
            AB flux per diamond [rad].
        J : float
            Hopping amplitude.
        T : float
            Total evolution time [units of ℏ/J].
        dt : float
            Time step for the output time axis.
        disorder_W : float
            On-site disorder strength.
        disorder_frac : float
            Fraction of sites that receive disorder.
        disorder_seed : int or None
            Seed for the disorder realisation (see build_hamiltonian).
        n_peaks : int or None
            Number of loaded peaks to record in the result dict.
            Defaults to counting non-zero A-sites in *psi0*.
        verbose : bool

        Returns
        -------
        dict with keys:
            times, fidelity, spread, localization,
            init_2d, snapshots, phi, n_peaks, flat_ok
        """
        Lx, Ly = self.Lx, self.Ly
        dim    = self.dim

        if verbose:
            print(f"\n  STAGE 4: 2D Diamond Caging  "
                  f"[phi={phi/np.pi:.2f}pi, {Lx}x{Ly}"
                  + (f", W={disorder_W:.2f}J" if disorder_W > 0 else "")
                  + "]")

        # ── Spectrum check ────────────────────────────────────────────
        flat_ok = True
        if abs(phi - np.pi) < 0.01 and disorder_W == 0:
            flat_ok, _, _ = self.validate_spectrum(
                phi=phi, J=J, abort_on_fail=False)

        # ── Normalise initial state ───────────────────────────────────
        psi_np = np.array(psi0, dtype=complex)
        norm   = np.linalg.norm(psi_np)
        if norm < 1e-12:
            psi_np[::5] = 1.0 / np.sqrt(Lx * Ly)
        else:
            psi_np = psi_np / norm
        psi0_np = psi_np.copy()

        # ── A-site index array ────────────────────────────────────────
        a_indices = np.array([self._site(ix, iy, 0)
                               for ix in range(Lx) for iy in range(Ly)])

        # ── n_peaks from non-zero A-sites if not provided ─────────────
        if n_peaks is None:
            n_peaks = int(np.sum(np.abs(psi0[a_indices]) > 1e-12))

        if verbose:
            print(f"    {n_peaks} peaks loaded onto A-sites")

        # ── Localisation mask: 3×3 neighbourhood of loaded A-sites ────
        local_mask_np = np.zeros(dim, dtype=bool)
        for ix_ in range(Lx):
            for iy_ in range(Ly):
                if abs(psi0[self._site(ix_, iy_, 0)]) > 1e-12:
                    for ddx in [-1, 0, 1]:
                        for ddy in [-1, 0, 1]:
                            base = 5 * (((ix_ + ddx) % Lx) * Ly
                                        + ((iy_ + ddy) % Ly))
                            local_mask_np[base:base + 5] = True

        # ── Site coordinate tensors ───────────────────────────────────
        ix_all = np.arange(Lx)
        iy_all = np.arange(Ly)
        IX, IY = np.meshgrid(ix_all, iy_all, indexing='ij')
        site_x_np = np.repeat(IX.ravel(), 5).astype(float)
        site_y_np = np.repeat(IY.ravel(), 5).astype(float)

        # ── Move to GPU ───────────────────────────────────────────────
        H_np = self.build_hamiltonian(phi=phi, J=J,
                                       disorder_W=disorder_W,
                                       disorder_frac=disorder_frac,
                                       disorder_seed=disorder_seed)
        H_t          = torch.tensor(H_np,        dtype=torch.complex128, device=_device)
        psi_t        = torch.tensor(psi_np,       dtype=torch.complex128, device=_device)
        psi0_t       = torch.tensor(psi0_np,      dtype=torch.complex128, device=_device)
        site_x_t     = torch.tensor(site_x_np,   dtype=torch.float64,    device=_device)
        site_y_t     = torch.tensor(site_y_np,   dtype=torch.float64,    device=_device)
        local_mask_t = torch.tensor(local_mask_np,                        device=_device)

        # ── GPU diagonalisation ───────────────────────────────────────
        evals_t, evecs_t = torch.linalg.eigh(H_t)

        # ── Vectorised time evolution ─────────────────────────────────
        coeff_t   = evecs_t.conj().T @ psi_t
        times_np  = np.arange(0, T, dt)
        times_t   = torch.tensor(times_np, dtype=torch.float64, device=_device)

        phase_factors = torch.exp(
            -1j * evals_t.unsqueeze(0) * times_t.unsqueeze(1))
        weighted  = phase_factors * coeff_t.unsqueeze(0)
        psi_all   = weighted @ evecs_t.T          # (N_times, dim)

        probs_all = torch.abs(psi_all) ** 2       # (N_times, dim)
        fid_all   = torch.abs(psi_all @ psi0_t.conj()) ** 2  # (N_times,)

        mu_x  = probs_all @ site_x_t
        mu_y  = probs_all @ site_y_t
        dx2   = (site_x_t.unsqueeze(0) - mu_x.unsqueeze(1)) ** 2
        dy2   = (site_y_t.unsqueeze(0) - mu_y.unsqueeze(1)) ** 2
        spr_all = torch.sqrt((probs_all * (dx2 + dy2)).sum(dim=1))

        loc_all = probs_all[:, local_mask_t].sum(dim=1)

        # ── CPU transfer ──────────────────────────────────────────────
        fid_h    = fid_all.cpu().numpy()
        spr_h    = spr_all.cpu().numpy()
        loc_h    = loc_all.cpu().numpy()
        probs_cpu = probs_all.cpu().numpy()

        # ── Snapshots (A-sites only) ──────────────────────────────────
        snap_t_vals = [0, T * 0.25, T * 0.5, T]
        snaps = []
        for st in snap_t_vals:
            tidx   = np.argmin(np.abs(times_np - st))
            probs_t = probs_cpu[tidx]
            ad     = probs_t[a_indices].reshape(Lx, Ly)
            snaps.append((times_np[tidx], ad.copy()))

        final_fid = fid_h[-1]
        final_loc = loc_h[-1]
        if verbose:
            print(f"    Fidelity:      {final_fid:.4f}  "
                  f"({'PRESERVED' if final_fid > 0.5 else 'decayed'})")
            print(f"    Localization:  {final_loc:.4f}  "
                  f"({'CAGED' if final_loc > 0.8 else 'spread'})")
            print(f"    Spread: {spr_h[0]:.2f} -> {spr_h[-1]:.2f}  "
                  f"(delta={spr_h[-1] - spr_h[0]:+.2f} cells)")

        init_2d = (np.abs(psi0_np[a_indices]) ** 2).reshape(Lx, Ly)

        return {
            'times':        times_np,
            'fidelity':     fid_h,
            'spread':       spr_h,
            'localization': loc_h,
            'init_2d':      init_2d,
            'snapshots':    snaps,
            'phi':          phi,
            'n_peaks':      n_peaks,
            'flat_ok':      flat_ok,
        }
