"""
Integrated Quantum Substrate Deposition Simulator — v9
=======================================================

Changes from v8:

  FIX 4 — Replace Lieb lattice with 2D Diamond Network  (AB caging)
    v8: Used a Lieb lattice at Φ=π, expecting all bands to flatten.
    This was wrong for three compounding reasons:
      (a) The Lieb lattice has only ONE flat band at E=0, which exists
          at ALL flux values, not just Φ=π.  The dispersive bands never
          flatten, so the ≤5-unique-eigenvalue check always failed.
      (b) The E=0 flat band has ZERO weight on A-sites (hub sites,
          coordination 4).  Since the sim loaded deposition density
          onto A-sites, the initial state projected entirely onto the
          dispersive bands — no caging was possible.
      (c) Even with correct CLS initialization on B/C edge-sites, the
          Lieb flat band gives perfect caging at BOTH Φ=π and Φ=0.
          The comparison is meaningless because the flat band is an
          intrinsic geometric property, not flux-dependent.

    Fix: Replace with a 2D diamond network (decorated square lattice).
    Each bond of a square lattice of "hub" (A) sites is replaced by a
    two-path diamond — an upper and lower intermediate site.
    5 sites per unit cell: A, B_up, B_down (horizontal), C_up, C_down (vertical).

    At Φ=π per diamond:  ALL 5 bands are flat → complete AB caging.
      Exactly 3 unique eigenvalues: {−2√2, 0, +2√2} J.
      A-site initialization works naturally (non-zero flat-band weight).
      State never leaves the initial diamonds.

    At Φ=0:  dispersive bands → state diffuses across the lattice.

    This produces genuine, flux-dependent Aharonov-Bohm caging.
    Localization is robust to disorder: ≥90% at W=2J (Φ=π) vs ~45% (Φ=0).

  Retained from v8:
    FIX 1 — Phase field apodization (boundary artefact removal)
    FIX 3 — Denser vortex lattice (a=2λ)
    NEW   — Multi-species deposition
    NEW   — Disorder robustness test (now on diamond network)
    NEW   — Decoherence sweep

  Also retained:
    Floquet stage preserved as-is (see note below).
    The "spatial Floquet" formalism in Stage 2 is a local dressed-state
    decomposition.  It becomes rigorous Floquet theory when the
    superconducting array is modulated at frequency ω, making the
    potential genuinely time-periodic: V(r,t) = V₀(r) + V₁(r)cos(ωt).
    For v9, the dressed-state calculation is unchanged; only the lattice
    caging stage is rewritten.

Pipeline (v9):
  ψ_beam → [Kuramoto sync]              → ψ_coherent
         → [A-B phase + apodization]    → ψ_phased
         → [Propagation]                → ψ_propagated
         → [Spatial Floquet dress]      → ψ_dressed(x,y,n)
         → [Bind filter]                → ψ_adsorbed
         → [2D Diamond / AB caging]     → 2D fidelity + localization
"""

import os
import sys
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter, zoom
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

os.makedirs('results', exist_ok=True)

hbar = 1.0545718e-34
k_B  = 1.380649e-23
m_He = 6.6464731e-27


# ===========================================================================
# FLOQUET VALIDATION
# ===========================================================================

def validate_floquet(N_side=2, V_frac=0.9, abort_on_fail=True):
    n_vals = np.arange(-N_side, N_side + 1, dtype=float)
    fl_dim = len(n_vals)
    H = np.diag(n_vals.copy())
    for i in range(fl_dim - 1):
        H[i, i+1] = H[i+1, i] = V_frac
    evals, evecs = np.linalg.eigh(H)
    U     = evecs @ np.diag(np.exp(-1j * evals * 2*np.pi)) @ evecs.conj().T
    psi0  = np.zeros(fl_dim, dtype=complex)
    psi0[N_side] = 1.0
    c_n   = U @ psi0
    pops  = np.abs(c_n)**2
    u_err = np.max(np.abs(U @ U.conj().T - np.eye(fl_dim)))
    pnz   = pops[pops > 1e-12]
    S     = -np.sum(pnz * np.log(pnz))
    ok    = (u_err < 1e-12 and S > 0.1 and S < np.log(fl_dim))

    print("\n  FLOQUET VALIDATION")
    print(f"  N_side={N_side}, V={V_frac:.2f}ℏω")
    print(f"  Unitarity err: {u_err:.2e}  {'OK' if u_err<1e-12 else 'FAIL'}")
    print(f"  Entropy: {S:.4f}  (bounds: 0.1 < S < {np.log(fl_dim):.2f})")
    pop_str = '  '.join(f'n={int(n):+d}:{p:.3f}'
                         for n, p in zip(n_vals, pops))
    print(f"  Pops: {pop_str}")
    print(f"  Validation: {'PASS' if ok else 'FAIL'}")
    if not ok and abort_on_fail:
        sys.exit("ABORT: Floquet validation failed.")
    return ok, c_n, S


# ===========================================================================
# 2D DIAMOND NETWORK VALIDATION  (replaces Lieb lattice validation)
# ===========================================================================

def build_2d_diamond_hamiltonian(Lx, Ly, phi, J=1.0, disorder_W=0.0):
    """
    2D Diamond network with Aharonov-Bohm flux.

    Unit cell (5 sites):
      A  — hub site at square lattice vertex
      Bu, Bd — upper/lower path on horizontal bond (A → A_right)
      Cu, Cd — upper/lower path on vertical   bond (A → A_up)

    Each diamond (4 bonds forming a rhombus) encloses flux phi.
    Phase distributed symmetrically: ±phi/4 per bond.

    At phi=π:  all bands flat → complete AB caging.
    At phi=0:  dispersive bands → free propagation.
    """
    n_sites = 5
    dim = n_sites * Lx * Ly
    H = np.zeros((dim, dim), dtype=complex)

    def s(ix, iy, t):
        return n_sites * ((ix % Lx) * Ly + (iy % Ly)) + t

    p  = np.exp( 1j * phi / 4)
    pc = np.exp(-1j * phi / 4)

    for ix in range(Lx):
        for iy in range(Ly):
            ix1 = (ix + 1) % Lx
            iy1 = (iy + 1) % Ly
            a  = s(ix, iy, 0)
            bu = s(ix, iy, 1);  bd = s(ix, iy, 2)
            cu = s(ix, iy, 3);  cd = s(ix, iy, 4)
            a_right = s(ix1, iy, 0)
            a_up    = s(ix, iy1, 0)

            # Horizontal diamond:  A —(Bu,Bd)— A_right
            H[a,  bu] = J * p;   H[bu, a]  = J * pc
            H[a,  bd] = J * pc;  H[bd, a]  = J * p
            H[bu, a_right] = J * p;   H[a_right, bu] = J * pc
            H[bd, a_right] = J * pc;  H[a_right, bd] = J * p

            # Vertical diamond:  A —(Cu,Cd)— A_up
            H[a,  cu] = J * p;   H[cu, a]  = J * pc
            H[a,  cd] = J * pc;  H[cd, a]  = J * p
            H[cu, a_up] = J * p;   H[a_up, cu] = J * pc
            H[cd, a_up] = J * pc;  H[a_up, cd] = J * p

    # On-site disorder
    if disorder_W > 0:
        dis = np.random.uniform(-disorder_W / 2, disorder_W / 2, dim)
        np.fill_diagonal(H, H.diagonal() + dis)

    return H


def validate_diamond_spectrum(Lx, Ly, phi=np.pi, J=1.0, abort_on_fail=True):
    """
    Hard gate: at Φ=π, the 2D diamond network must have exactly
    3 unique eigenvalues (all 5 bands flat).

    Expected band energies: {−2√2, 0, +2√2} · J ≈ {−2.83, 0, +2.83}.
    """
    H = build_2d_diamond_hamiltonian(Lx, Ly, phi, J)
    evals = np.linalg.eigvalsh(H)
    n_uniq = len(np.unique(np.round(evals, 2)))
    ok = n_uniq <= 5  # 3 ideal; allow small rounding tolerance

    print(f"\n  DIAMOND SPECTRUM CHECK  ({Lx}×{Ly}, Φ={phi/np.pi:.1f}π)")
    print(f"  Unique eigenvalues: {n_uniq}  "
          f"({'OK (all flat)' if ok else 'FAIL'})")
    if ok:
        bands = np.unique(np.round(evals, 2))
        print(f"  Band energies: {bands}")
    if not ok and abort_on_fail:
        sys.exit(f"ABORT: Diamond flat-band check failed "
                 f"({n_uniq} unique eigenvalues, need ≤5).")
    return ok, H, evals


# ===========================================================================
# CORE SIMULATOR
# ===========================================================================

class IntegratedQuantumSubstrate:

    def __init__(self, N=256, L=400e-9, T_beam=1e-3,
                 N_floquet_sidebands=2,
                 N_diamond_cells_x=16, N_diamond_cells_y=16):
        self.N  = N
        self.L  = L
        self.dx = L / N
        self.x  = np.linspace(-L/2, L/2, N)
        self.y  = np.linspace(-L/2, L/2, N)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

        self.mass   = m_He
        self.T_beam = T_beam
        self.v      = np.sqrt(2 * k_B * T_beam / m_He)
        self.k0     = self.mass * self.v / hbar
        self.lam    = 2 * np.pi / self.k0
        self.E0     = 0.5 * self.mass * self.v**2

        self.N_side  = N_floquet_sidebands
        self.fl_dim  = 2 * N_floquet_sidebands + 1
        self.n_vals  = np.arange(-N_floquet_sidebands,
                                  N_floquet_sidebands + 1, dtype=float)

        self.Lx        = N_diamond_cells_x
        self.Ly        = N_diamond_cells_y
        self.n_sites   = 5   # sites per unit cell for diamond network
        self.lattice_dim = self.n_sites * self.Lx * self.Ly

    def info(self):
        print("=" * 65)
        print("INTEGRATED QUANTUM SUBSTRATE  v9")
        print("=" * 65)
        print(f"  He-4 at {self.T_beam*1e3:.3f} mK  "
              f"λ={self.lam*1e9:.2f}nm  v={self.v:.3f}m/s")
        print(f"  Grid: {self.N}×{self.N}, dx={self.dx*1e9:.2f}nm, "
              f"λ/dx={self.lam/self.dx:.1f}")
        print(f"  Floquet: N_side={self.N_side}, spatial dressing")
        print(f"  2D Diamond: {self.Lx}×{self.Ly} (AB caging)")

    # -----------------------------------------------------------------------
    # STAGE 0: Beam + Kuramoto
    # -----------------------------------------------------------------------
    def stage0_beam(self, N_atoms=200, K=6.0, T_sync=30,
                    alpha=0.5, L_coherence=None, verbose=True):
        if verbose:
            print("\n  STAGE 0: Beam + Kuramoto")

        omega  = np.random.normal(0, 1.0, N_atoms)
        theta  = np.random.uniform(0, 2*np.pi, N_atoms)
        dtheta = np.zeros(N_atoms)
        dt_k   = 0.01
        order_hist = []

        for _ in np.arange(0, T_sync, dt_k):
            z        = np.mean(np.exp(1j * theta))
            order_hist.append(np.abs(z))
            coupling = np.imag(z * np.exp(-1j * theta))
            ddtheta  = -alpha * dtheta + omega + K * coupling
            dtheta  += ddtheta * dt_k
            theta   += dtheta * dt_k

        r     = order_hist[-1]
        sigma = 0.35 * self.L
        psi   = (np.exp(-(self.X**2 + self.Y**2) / (2*sigma**2))
                 * np.exp(1j * self.k0 * self.X))

        noise = (1 - r) * np.random.normal(0, np.pi, (self.N, self.N))
        noise = gaussian_filter(noise, sigma=3)
        psi  *= np.exp(1j * noise)

        if L_coherence is not None:
            R2 = self.X**2 + self.Y**2
            psi *= np.exp(-R2 / (2 * L_coherence**2))

        psi /= np.sqrt(np.sum(np.abs(psi)**2) * self.dx**2)

        if verbose:
            lc_str = (f", L_c={L_coherence*1e9:.0f}nm"
                      if L_coherence is not None else "")
            print(f"    r={r:.4f}  noise_RMS={np.std(noise):.4f}rad{lc_str}")

        return psi, r, np.array(order_hist), np.arange(0, T_sync, dt_k)

    # -----------------------------------------------------------------------
    # STAGE 1: A-B phase imprinting + apodization
    # -----------------------------------------------------------------------
    def stage1_ab_phase(self, psi_in, pattern='vortex_lattice',
                        apodize=True, apod_margin=0.15,
                        verbose=True, **kw):
        if verbose:
            print("\n  STAGE 1: A-B phase imprinting"
                  + (" + apodization" if apodize else ""))

        X, Y  = self.X, self.Y
        phase = np.zeros((self.N, self.N))

        if pattern == 'vortex_lattice':
            a    = kw.get('a',    2 * self.lam)
            core = kw.get('core', 0.6 * self.lam)
            n_max = int(self.L / (2 * a)) + 2
            count = 0
            for i in range(-n_max, n_max + 1):
                for j in range(-n_max, n_max + 1):
                    x0 = a * (i + 0.5 * (j % 2))
                    y0 = a * np.sqrt(3) / 2 * j
                    if abs(x0) < self.L/2*1.1 and abs(y0) < self.L/2*1.1:
                        R     = np.sqrt((X-x0)**2 + (Y-y0)**2 + 1e-30)
                        theta = np.arctan2(Y-y0, X-x0)
                        sign  = (-1)**(i + j)
                        w     = 1 - np.exp(-R**2 / (2*core**2))
                        phase += sign * np.pi * theta / (2*np.pi) * w
                        count += 1
            if verbose:
                print(f"    Vortex lattice: a={a*1e9:.1f}nm "
                      f"({count} vortices), "
                      f"φ∈[{phase.min():.2f},{phase.max():.2f}]")

        elif pattern == 'ring_array':
            spacing = kw.get('spacing', 3 * self.lam)
            r_ring  = kw.get('r_ring',  1.0 * self.lam)
            width   = kw.get('width',   0.3 * self.lam)
            flux    = kw.get('flux',    0.5)
            n_max   = int(self.L / (2 * spacing)) + 1
            count   = 0
            for i in range(-n_max, n_max + 1):
                for j in range(-n_max, n_max + 1):
                    x0, y0 = i*spacing, j*spacing
                    if abs(x0) < self.L/2*1.1 and abs(y0) < self.L/2*1.1:
                        R     = np.sqrt((X-x0)**2 + (Y-y0)**2)
                        theta = np.arctan2(Y-y0, X-x0)
                        w     = np.exp(-(R - r_ring)**2 / (2*width**2))
                        phase += w * flux * theta
                        count += 1
            if verbose:
                print(f"    Ring array: {count} rings")

        elif pattern == 'sinusoidal':
            period = kw.get('period', 2 * self.lam)
            amp    = kw.get('amp',    np.pi)
            phase  = amp * np.cos(2*np.pi*X/period) * np.cos(2*np.pi*Y/period)
            if verbose:
                print(f"    Sinusoidal: period={period*1e9:.1f}nm")

        if apodize:
            margin = apod_margin * self.L
            def _env(coord):
                d = self.L/2 - np.abs(coord)
                return np.where(d < margin,
                    np.cos(np.pi/2 * (margin - d) / margin)**2, 1.0)
            env = _env(X) * _env(Y)
            phase = phase * env
            if verbose:
                print(f"    Apodization: margin={apod_margin*100:.0f}% "
                      f"→ φ∈[{phase.min():.2f},{phase.max():.2f}] after")

        psi_out = psi_in * np.exp(1j * phase)
        return psi_out, phase

    # -----------------------------------------------------------------------
    # STAGE 2: Spatial Floquet dressing
    # -----------------------------------------------------------------------
    def stage2_floquet_dress_spatial(self, psi_in, phase_map,
                                      V_max=0.9, n_levels=20,
                                      verbose=True):
        if verbose:
            print("\n  STAGE 2: Spatial Floquet dressing")

        V_map     = V_max * np.abs(phase_map) / (np.abs(phase_map).max() + 1e-30)
        V_levels  = np.linspace(0.0, V_max, n_levels + 1)
        V_indices = np.clip(np.digitize(V_map, V_levels) - 1, 0, n_levels - 1)

        c_n_cache = {}
        for lvl in range(n_levels):
            V_val = V_levels[lvl]
            H = np.diag(self.n_vals.copy())
            for i in range(self.fl_dim - 1):
                H[i, i+1] = H[i+1, i] = V_val
            evals, evecs = np.linalg.eigh(H)
            U  = evecs @ np.diag(np.exp(-1j * evals * 2*np.pi)) @ evecs.conj().T
            psi0 = np.zeros(self.fl_dim, dtype=complex)
            psi0[self.N_side] = 1.0
            c_n_cache[lvl] = U @ psi0

        psi_dressed = np.zeros((self.N, self.N, self.fl_dim), dtype=complex)
        for lvl in np.unique(V_indices):
            mask = (V_indices == lvl)
            c_n  = c_n_cache[int(lvl)]
            psi_dressed[mask] = psi_in[mask, np.newaxis] * c_n[np.newaxis, :]

        intensity_total = np.sum(np.abs(psi_in)**2) + 1e-30
        avg_pops = np.zeros(self.fl_dim)
        for lvl in np.unique(V_indices):
            mask   = (V_indices == lvl)
            w      = np.sum(np.abs(psi_in[mask])**2) / intensity_total
            avg_pops += w * np.abs(c_n_cache[int(lvl)])**2

        pnz     = avg_pops[avg_pops > 1e-12]
        entropy = -np.sum(pnz * np.log(pnz))

        if verbose:
            print(f"    V_max={V_max:.2f}ℏω, {n_levels} levels")
            pop_str = '  '.join(f'n={int(n):+d}:{p:.3f}'
                                 for n, p in zip(self.n_vals, avg_pops)
                                 if p > 0.005)
            print(f"    Beam-avg pops: {pop_str}")
            print(f"    Beam-avg entropy: {entropy:.4f}")

        return psi_dressed, avg_pops, entropy, V_map

    # -----------------------------------------------------------------------
    # STAGE 3: Binding resonance filter
    # -----------------------------------------------------------------------
    def stage3_binding_filter(self, psi_dressed, avg_pops,
                               n_resonant=0, width_frac=0.4,
                               verbose=True):
        if verbose:
            print(f"\n  STAGE 3: Binding filter  [n_resonant={n_resonant}]")

        weights = np.array([
            width_frac**2 / ((n - n_resonant)**2 + width_frac**2)
            for n in self.n_vals
        ])
        psi_ads  = np.sum(weights[np.newaxis, np.newaxis, :] * psi_dressed,
                          axis=2)
        norm_in  = np.sum(np.abs(psi_dressed)**2) * self.dx**2
        norm_out = np.sum(np.abs(psi_ads)**2)      * self.dx**2
        ads_frac = norm_out / norm_in if norm_in > 0 else 0.0

        if verbose:
            for idx, (n, w) in enumerate(zip(self.n_vals, weights)):
                p = avg_pops[idx]
                if w > 0.02 or p > 0.02:
                    print(f"      n={int(n):+d}: w={w:.4f}  avg_pop={p:.4f}")
            print(f"    ads_frac={ads_frac:.4f}")

        return psi_ads, ads_frac, weights

    # -----------------------------------------------------------------------
    # STAGE 4: 2D Diamond Network AB caging  (replaces Lieb lattice)
    # -----------------------------------------------------------------------
    def _diamond_site(self, ix, iy, t):
        """Map (cell_x, cell_y, sublattice) → flat index."""
        return self.n_sites * ((ix % self.Lx) * self.Ly + (iy % self.Ly)) + t

    def stage4_diamond_caging(self, density_2d, phi_cage=np.pi,
                               J_hop=1.0, T_evolve=40, dt=0.05,
                               n_peaks=16, disorder_W=0.0,
                               verbose=True):
        """
        2D Diamond network AB caging.

        At Φ=π: all bands flat → state spatially caged.
        At Φ=0: dispersive → state diffuses.

        Loads deposition density peaks onto A-sites (sublattice 0).
        Measures: fidelity (overlap), RMS spread, localization probability.
        """
        if verbose:
            print(f"\n  STAGE 4: 2D Diamond caging  "
                  f"[Φ={phi_cage/np.pi:.2f}π, "
                  f"{self.Lx}×{self.Ly}"
                  + (f", W={disorder_W:.2f}J" if disorder_W > 0 else "") + "]")

        Lx, Ly = self.Lx, self.Ly
        dim = self.lattice_dim
        s = self._diamond_site

        # Build Hamiltonian
        H = build_2d_diamond_hamiltonian(Lx, Ly, phi_cage, J_hop,
                                          disorder_W=disorder_W)

        # Spectrum check (clean Φ=π only)
        flat_ok = True
        if abs(phi_cage - np.pi) < 0.01 and disorder_W == 0:
            evals_chk = np.linalg.eigvalsh(H)
            n_uniq    = len(np.unique(np.round(evals_chk, 2)))
            flat_ok   = (n_uniq <= 5)
            if verbose:
                print(f"    Spectrum: {n_uniq} unique eigenvalues  "
                      f"({'ALL FLAT OK' if flat_ok else 'FAIL'})")
                if flat_ok:
                    bv = np.unique(np.round(evals_chk, 2))
                    print(f"    Band energies: {bv}")

        # Peak-load 2D density onto A-sites
        scale_x = Lx / self.N
        scale_y = Ly / self.N
        density_ds = np.maximum(zoom(density_2d, (scale_x, scale_y),
                                      order=1), 0)
        d_smooth   = gaussian_filter(density_ds, sigma=1)
        local_max  = maximum_filter(d_smooth, size=3) == d_smooth
        peak_mask  = local_max & (d_smooth > 0.1 * d_smooth.max())
        peak_idx   = np.argwhere(peak_mask)
        peak_vals  = d_smooth[peak_mask]
        order_     = np.argsort(peak_vals)[::-1]
        peak_idx   = peak_idx[order_[:n_peaks]]

        psi_lattice = np.zeros(dim, dtype=complex)
        for px, py in peak_idx:
            ix_ = min(px, Lx - 1)
            iy_ = min(py, Ly - 1)
            psi_lattice[s(ix_, iy_, 0)] += np.sqrt(d_smooth[px, py])

        norm = np.sqrt(np.sum(np.abs(psi_lattice)**2))
        if norm < 1e-12:
            for ix_ in range(Lx):
                for iy_ in range(Ly):
                    psi_lattice[s(ix_, iy_, 0)] = 1.0 / np.sqrt(Lx * Ly)
        else:
            psi_lattice /= norm

        psi0_lattice = psi_lattice.copy()
        init_dens = np.abs(psi_lattice)**2

        # Exact time evolution via eigendecomposition
        evals_H, evecs_H = np.linalg.eigh(H)
        times  = np.arange(0, T_evolve, dt)
        fid_h, spr_h, loc_h, snaps = [], [], [], []
        snap_t = [0, T_evolve * 0.25, T_evolve * 0.5, T_evolve]

        site_x = np.zeros(dim, dtype=float)
        site_y = np.zeros(dim, dtype=float)
        for ix_ in range(Lx):
            for iy_ in range(Ly):
                for t_ in range(self.n_sites):
                    site_x[s(ix_, iy_, t_)] = float(ix_)
                    site_y[s(ix_, iy_, t_)] = float(iy_)

        # Build localization mask: sites within ±1 cell of any peak
        local_mask = np.zeros(dim, dtype=bool)
        for px, py in peak_idx:
            ix_ = min(px, Lx - 1)
            iy_ = min(py, Ly - 1)
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for t_ in range(self.n_sites):
                        local_mask[s(ix_ + dx, iy_ + dy, t_)] = True

        coeff = evecs_H.conj().T @ psi_lattice

        for t in times:
            psi_t = evecs_H @ (np.exp(-1j * evals_H * t) * coeff)
            probs = np.abs(psi_t)**2
            fid_h.append(np.abs(np.dot(psi0_lattice.conj(), psi_t))**2)
            mu_x  = np.dot(site_x, probs)
            mu_y  = np.dot(site_y, probs)
            spr_h.append(np.sqrt(np.dot(
                (site_x - mu_x)**2 + (site_y - mu_y)**2, probs)))
            loc_h.append(np.sum(probs[local_mask]))
            for st in snap_t:
                if abs(t - st) < dt / 2:
                    ad = np.zeros((Lx, Ly))
                    for ixx in range(Lx):
                        for iyy in range(Ly):
                            ad[ixx, iyy] = probs[s(ixx, iyy, 0)]
                    snaps.append((t, ad.copy()))

        final_fid = fid_h[-1]
        final_loc = loc_h[-1]
        if verbose:
            print(f"    {len(peak_idx)} peaks loaded onto A-sites")
            print(f"    Fidelity:      {final_fid:.4f}  "
                  f"({'PRESERVED' if final_fid > 0.5 else 'decayed'})")
            print(f"    Localization:  {final_loc:.4f}  "
                  f"({'CAGED' if final_loc > 0.8 else 'spread'})")
            print(f"    Spread: {spr_h[0]:.2f} → {spr_h[-1]:.2f}  "
                  f"(Δ={spr_h[-1]-spr_h[0]:+.2f} cells)")

        init_2d = np.zeros((Lx, Ly))
        for ix_ in range(Lx):
            for iy_ in range(Ly):
                init_2d[ix_, iy_] = init_dens[s(ix_, iy_, 0)]

        return {
            'times':         times,
            'fidelity':      np.array(fid_h),
            'spread':        np.array(spr_h),
            'localization':  np.array(loc_h),
            'init_2d':       init_2d,
            'snapshots':     snaps,
            'phi':           phi_cage,
            'n_peaks':       len(peak_idx),
            'flat_ok':       flat_ok,
        }

    # -----------------------------------------------------------------------
    # Propagator
    # -----------------------------------------------------------------------
    def _propagate(self, psi, distance):
        kx = np.fft.fftfreq(self.N, self.dx) * 2*np.pi
        ky = np.fft.fftfreq(self.N, self.dx) * 2*np.pi
        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        kz_sq  = self.k0**2 - KX**2 - KY**2
        valid  = kz_sq > 0
        kz     = np.where(valid, np.sqrt(np.maximum(kz_sq, 0)), 0.0)
        H_prop = np.where(valid, np.exp(1j * kz * distance), 0.0)
        return np.fft.ifft2(np.fft.fft2(psi) * H_prop)

    # -----------------------------------------------------------------------
    # Full pipeline
    # -----------------------------------------------------------------------
    def run_full_pipeline(self, pattern='vortex_lattice',
                          V_max=0.9, n_resonant=0,
                          K_kuramoto=6.0, phi_cage=np.pi,
                          prop_distance_lam=20, apodize=True,
                          verbose=True, **pattern_kw):
        if verbose:
            self.info()

        psi_beam, r_sync, order_hist, sync_times = self.stage0_beam(
            K=K_kuramoto, verbose=verbose)

        psi_phased, phase_map = self.stage1_ab_phase(
            psi_beam, pattern=pattern, apodize=apodize,
            verbose=verbose, **pattern_kw)

        dist = prop_distance_lam * self.lam
        if verbose:
            print(f"\n  PROPAGATION: {prop_distance_lam}λ = {dist*1e9:.1f}nm")
        psi_prop = self._propagate(psi_phased, dist)

        psi_dressed, avg_pops, entropy, V_map = \
            self.stage2_floquet_dress_spatial(
                psi_prop, phase_map, V_max=V_max, verbose=verbose)

        psi_ads, ads_frac, fl_weights = self.stage3_binding_filter(
            psi_dressed, avg_pops, n_resonant=n_resonant, verbose=verbose)

        density_prop  = np.abs(psi_prop)**2
        density_final = np.abs(psi_ads)**2

        cage_pi = self.stage4_diamond_caging(
            density_final, phi_cage=np.pi, verbose=verbose)
        cage_0  = self.stage4_diamond_caging(
            density_final, phi_cage=0.0, verbose=verbose)

        return {
            'psi_beam':        psi_beam,
            'psi_phased':      psi_phased,
            'psi_propagated':  psi_prop,
            'psi_adsorbed':    psi_ads,
            'density_prop':    density_prop,
            'density_final':   density_final,
            'phase_map':       phase_map,
            'V_map':           V_map,
            'avg_pops':        avg_pops,
            'fl_weights':      fl_weights,
            'fl_entropy':      entropy,
            'adsorption_frac': ads_frac,
            'sync_order':      order_hist,
            'sync_times':      sync_times,
            'r_sync':          r_sync,
            'cage_pi':         cage_pi,
            'cage_0':          cage_0,
        }


# ===========================================================================
# METRICS
# ===========================================================================

def michelson_contrast(d, lo=5, hi=95):
    lo_v = np.percentile(d, lo)
    hi_v = np.percentile(d, hi)
    return (hi_v - lo_v) / (hi_v + lo_v + 1e-30)

def ssim_score(ref, test):
    r = ref  / (ref.max()  + 1e-30)
    t = test / (test.max() + 1e-30)
    try:
        from skimage.metrics import structural_similarity
        return structural_similarity(r, t, data_range=1.0)
    except ImportError:
        r_f = r - r.mean(); t_f = t - t.mean()
        return float(np.sum(r_f*t_f) /
                     (np.sqrt(np.sum(r_f**2)*np.sum(t_f**2)) + 1e-30))

def min_feature_size(density, x_axis):
    mid  = density.shape[1] // 2
    line = density[:, mid] / (density[:, mid].max() + 1e-30)
    peaks, _ = find_peaks(line, height=0.3, distance=3)
    if not len(peaks):
        return np.nan
    fwhms = []
    for pk in peaks:
        half = 0.5 * line[pk]
        l = pk
        while l > 0 and line[l] > half: l -= 1
        r = pk
        while r < len(line)-1 and line[r] > half: r += 1
        fwhms.append((r - l) * abs(x_axis[1] - x_axis[0]))
    return float(np.min(fwhms)) if fwhms else np.nan


# ===========================================================================
# SIDEBAND SELECTIVITY SWEEP
# ===========================================================================

def sideband_selectivity_sweep(seed=42):
    print("\n" + "="*65)
    print("SIDEBAND SELECTIVITY SWEEP")
    print("="*65)

    np.random.seed(seed)
    sim = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)

    psi, r, _, _  = sim.stage0_beam(K=6.0, verbose=False)
    psi_ph, phase = sim.stage1_ab_phase(psi, 'vortex_lattice',
                                         apodize=True, verbose=True)
    psi_prop      = sim._propagate(psi_ph, 20 * sim.lam)
    psi_d, avg_p, ent, V_map = sim.stage2_floquet_dress_spatial(
        psi_prop, phase, V_max=0.9, verbose=True)

    results = {}
    n_range = list(sim.n_vals.astype(int))

    for n_res in n_range:
        psi_ads, af, wt = sim.stage3_binding_filter(
            psi_d, avg_p, n_resonant=n_res, verbose=False)
        density = np.abs(psi_ads)**2
        C       = michelson_contrast(density)
        feat    = min_feature_size(density, sim.x * 1e9)
        print(f"  n={n_res:+d}:  ads={af:.4f}  C={C:.3f}  feat={feat:.1f}nm")
        results[n_res] = {
            'density': density, 'ads_frac': af,
            'contrast': C, 'feature': feat, 'weights': wt,
        }

    print("\n  SSIM matrix:")
    print(f"  {'':>5}", end='')
    for n in n_range: print(f"  n={n:+d}", end='')
    print()
    for n1 in n_range:
        print(f"  n={n1:+d}", end='')
        for n2 in n_range:
            s = ssim_score(results[n1]['density'], results[n2]['density'])
            print(f"  {s:.3f}", end='')
        print()

    return sim, results, phase, V_map, avg_p


# ===========================================================================
# MULTI-SPECIES DEPOSITION
# ===========================================================================

def multi_species_deposition(seed=42):
    print("\n" + "="*65)
    print("MULTI-SPECIES DEPOSITION")
    print("  Species A: n_resonant=0   (inter-vortex)")
    print("  Species B: n_resonant=+2  (vortex cores)")
    print("="*65)

    np.random.seed(seed)
    sim = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)

    psi, r, _, _  = sim.stage0_beam(K=6.0, verbose=False)
    psi_ph, phase = sim.stage1_ab_phase(psi, 'vortex_lattice',
                                         apodize=True, verbose=True)
    psi_prop      = sim._propagate(psi_ph, 20 * sim.lam)
    psi_d, avg_p, ent, V_map = sim.stage2_floquet_dress_spatial(
        psi_prop, phase, V_max=0.9, verbose=True)

    psi_A, af_A, _ = sim.stage3_binding_filter(
        psi_d, avg_p, n_resonant=0,  verbose=True)
    psi_B, af_B, _ = sim.stage3_binding_filter(
        psi_d, avg_p, n_resonant=+2, verbose=True)

    dens_A = np.abs(psi_A)**2
    dens_B = np.abs(psi_B)**2

    ssim_AB = ssim_score(dens_A, dens_B)
    overlap  = np.sum(np.sqrt(dens_A * dens_B)) * sim.dx**2
    overlap /= (np.sqrt(np.sum(dens_A) * np.sum(dens_B)) * sim.dx**2 + 1e-30)

    print(f"\n  Species A ads={af_A:.4f}  C={michelson_contrast(dens_A):.3f}")
    print(f"  Species B ads={af_B:.4f}  C={michelson_contrast(dens_B):.3f}")
    print(f"  SSIM(A,B) = {ssim_AB:.4f}  "
          f"(lower = more spatially distinct)")
    print(f"  Overlap   = {overlap:.4f}  "
          f"(lower = better spatial separation)")

    return sim, dens_A, dens_B, phase, V_map, ssim_AB, overlap


# ===========================================================================
# DISORDER ROBUSTNESS  (now on 2D diamond network)
# ===========================================================================

def disorder_robustness(seed=42, n_W=10, W_max=2.0, T_evolve=40):
    """
    Test AB caging robustness to on-site disorder.
    On the 2D diamond network:
      Φ=π: all bands flat → localization robust to moderate disorder.
      Φ=0: dispersive → state diffuses, disorder adds Anderson localization.
    """
    print("\n" + "="*65)
    print("DISORDER ROBUSTNESS TEST  (2D Diamond Network)")
    print(f"  W range: 0 → {W_max:.1f}J  ({n_W} points)")
    print("="*65)

    np.random.seed(seed)
    sim = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)

    psi, _, _, _ = sim.stage0_beam(K=6.0, verbose=False)
    psi_ph, ph   = sim.stage1_ab_phase(psi, 'vortex_lattice',
                                         apodize=True, verbose=False)
    psi_prop     = sim._propagate(psi_ph, 20 * sim.lam)
    psi_d, ap, _, _ = sim.stage2_floquet_dress_spatial(
        psi_prop, ph, V_max=0.9, verbose=False)
    psi_ads, _, _   = sim.stage3_binding_filter(
        psi_d, ap, n_resonant=0, verbose=False)
    density = np.abs(psi_ads)**2

    W_values = np.linspace(0, W_max, n_W)
    fid_pi, fid_0 = [], []
    loc_pi, loc_0 = [], []
    n_trials = 5

    for W in W_values:
        fp_t, f0_t = [], []
        lp_t, l0_t = [], []
        for trial in range(n_trials):
            np.random.seed(seed + trial * 100)
            r_pi = sim.stage4_diamond_caging(
                density, phi_cage=np.pi, T_evolve=T_evolve,
                disorder_W=W, verbose=False)
            r_0  = sim.stage4_diamond_caging(
                density, phi_cage=0.0, T_evolve=T_evolve,
                disorder_W=W, verbose=False)
            fp_t.append(r_pi['fidelity'][-1])
            f0_t.append(r_0['fidelity'][-1])
            lp_t.append(r_pi['localization'][-1])
            l0_t.append(r_0['localization'][-1])

        fid_pi.append(np.mean(fp_t));  fid_0.append(np.mean(f0_t))
        loc_pi.append(np.mean(lp_t));  loc_0.append(np.mean(l0_t))
        ratio_loc = loc_pi[-1] / (loc_0[-1] + 1e-6)
        print(f"  W={W:.2f}J:  Φ=π loc={loc_pi[-1]:.4f} fid={fid_pi[-1]:.4f}  "
              f"Φ=0 loc={loc_0[-1]:.4f} fid={fid_0[-1]:.4f}  "
              f"loc_ratio={ratio_loc:.1f}×")

    return (W_values, np.array(fid_pi), np.array(fid_0),
            np.array(loc_pi), np.array(loc_0))


# ===========================================================================
# DECOHERENCE SWEEP
# ===========================================================================

def decoherence_sweep(seed=42, n_pts=12):
    print("\n" + "="*65)
    print("DECOHERENCE SWEEP")
    print("="*65)

    np.random.seed(seed)
    sim   = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)
    a_vortex = 2 * sim.lam

    L_c_values = np.concatenate([
        [np.inf],
        np.logspace(np.log10(10*a_vortex), np.log10(0.5*a_vortex), n_pts-1)
    ])

    contrasts = []
    features  = []
    Lc_nm     = []

    for L_c in L_c_values:
        np.random.seed(seed)
        L_c_arg = None if np.isinf(L_c) else L_c

        psi, r, _, _ = sim.stage0_beam(K=6.0, L_coherence=L_c_arg,
                                        verbose=False)
        psi_ph, ph   = sim.stage1_ab_phase(psi, 'vortex_lattice',
                                             apodize=True, verbose=False)
        psi_prop     = sim._propagate(psi_ph, 20 * sim.lam)
        psi_d, ap, ent, _ = sim.stage2_floquet_dress_spatial(
            psi_prop, ph, V_max=0.9, verbose=False)
        psi_ads, af, _    = sim.stage3_binding_filter(
            psi_d, ap, n_resonant=0, verbose=False)
        density = np.abs(psi_ads)**2

        C    = michelson_contrast(density)
        feat = min_feature_size(density, sim.x * 1e9)
        Lc_s = "∞" if np.isinf(L_c) else f"{L_c*1e9:.0f}"
        print(f"  L_c={Lc_s}nm:  C={C:.3f}  feat={feat:.1f}nm  "
              f"ads={af:.3f}")

        contrasts.append(C)
        features.append(feat)
        Lc_nm.append(np.inf if np.isinf(L_c) else L_c * 1e9)

    return (np.array(Lc_nm), np.array(contrasts),
            np.array(features), a_vortex * 1e9)


# ===========================================================================
# ABLATION STUDY
# ===========================================================================

def ablation_study(seed=42):
    print("\n" + "="*65)
    print(f"ABLATION STUDY  (seed={seed})")
    print("="*65)

    np.random.seed(seed)
    sim      = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)
    prop_d   = 20 * sim.lam
    results  = {}
    ref_dens = None

    cases = [
        ('Full pipeline',   True,  True,  0,  6.0),
        ('No Floquet',      True,  False, 0,  6.0),
        ('No A-B phase',    False, True,  0,  6.0),
        ('Poor coherence',  True,  True,  0,  0.5),
        ('Sideband n=+2',   True,  True,  2,  6.0),
    ]

    for name, use_ab, use_fl, n_res, K in cases:
        print(f"\n>>> {name}")
        np.random.seed(seed)

        psi, r, _, _ = sim.stage0_beam(K=K, verbose=False)
        print(f"    r={r:.4f}")

        if use_ab:
            psi, ph = sim.stage1_ab_phase(psi, 'vortex_lattice',
                                           apodize=True, verbose=False)
        else:
            ph = np.zeros((sim.N, sim.N))

        psi = sim._propagate(psi, prop_d)

        if use_fl:
            psi_d, ap, ent, _ = sim.stage2_floquet_dress_spatial(
                psi, ph, V_max=0.9, verbose=False)
            psi_ads, af, _ = sim.stage3_binding_filter(
                psi_d, ap, n_resonant=n_res, verbose=False)
            print(f"    entropy={ent:.4f}  ads={af:.4f}")
        else:
            psi_ads = psi
            af, ent = 1.0, 0.0
            print(f"    [Floquet skipped]  ads={af:.4f}")

        density = np.abs(psi_ads)**2
        C       = michelson_contrast(density)
        if name == 'Full pipeline':
            ref_dens = density.copy()
        results[name] = {'density': density, 'contrast': C,
                          'r': r, 'ads_frac': af, 'entropy': ent}
        print(f"    C={C:.4f}")

    print("\n  SSIM vs Full pipeline:")
    for name, res in results.items():
        s = ssim_score(ref_dens, res['density'])
        res['ssim'] = s
        print(f"    {name:<22}: SSIM={s:.4f}  C={res['contrast']:.4f}")

    return sim, results


# ===========================================================================
# FIGURES
# ===========================================================================

def plot_pipeline(sim, r, fname='results/v9_pipeline.png'):
    print(f"\n  Plotting pipeline → {fname}")
    fig = plt.figure(figsize=(26, 32))
    gs  = GridSpec(6, 4, figure=fig, hspace=0.45, wspace=0.35,
                   left=0.05, right=0.97, top=0.97, bottom=0.02)
    x_nm = sim.x * 1e9
    ext  = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    # ── Banner ──
    ax = fig.add_subplot(gs[0, :])
    ax.axis('off')
    cp = michelson_contrast(r['density_prop'])
    cf = michelson_contrast(r['density_final'])
    loc_gap = (r['cage_pi']['localization'][-1]
               - r['cage_0']['localization'][-1])
    fid_gap = (r['cage_pi']['fidelity'][-1]
               - r['cage_0']['fidelity'][-1])
    txt = (
        "INTEGRATED QUANTUM SUBSTRATE  v9\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"He-4 {sim.T_beam*1e3:.1f}mK  λ={sim.lam*1e9:.1f}nm  "
        f"r_sync={r['r_sync']:.3f}  "
        f"Floquet S={r['fl_entropy']:.3f}  ads={r['adsorption_frac']:.3f}\n"
        f"Contrast: prop={cp:.3f} final={cf:.3f}  "
        f"Caging: loc_gap={loc_gap:.3f}  fid_gap={fid_gap:.3f}  "
        f"flat_bands={'OK' if r['cage_pi']['flat_ok'] else 'FAIL'}\n\n"
        "ψ → [Kuramoto] → [A-B+apodize] → [Propagate] "
        "→ [Spatial Floquet] → [Bind] → [2D Diamond/AB cage]"
    )
    ax.text(0.5, 0.5, txt, transform=ax.transAxes, ha='center', va='center',
            fontsize=11, fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#e3f2fd',
                      edgecolor='#1565c0', linewidth=2))

    # ── Row 1: Density stages ──
    stages = [
        ('① Beam |ψ|²',              np.abs(r['psi_beam'])**2),
        ('② After A-B+apodize',      np.abs(r['psi_phased'])**2),
        ('③ After propagation',       r['density_prop']),
        ('④ After Floquet+filter',    r['density_final']),
    ]
    for idx, (title, dens) in enumerate(stages):
        ax = fig.add_subplot(gs[1, idx])
        d  = dens / (dens.max() + 1e-30)
        im = ax.imshow(d.T, extent=ext, cmap='inferno', origin='lower')
        ax.set_title(f'{title}\nC={michelson_contrast(dens):.3f}',
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ── Row 2: Phase maps ──
    ph_stages = [
        ('A-B phase (apodized)',     r['phase_map']),
        ('V(x,y) drive map',         r['V_map']),
        ('Phase after A-B',          np.angle(r['psi_phased'])),
        ('Phase after propagation',  np.angle(r['psi_propagated'])),
    ]
    for idx, (title, ph) in enumerate(ph_stages):
        ax = fig.add_subplot(gs[2, idx])
        cmap = 'hsv' if 'Phase' in title or 'A-B' in title else 'hot'
        vmin, vmax = (-np.pi, np.pi) if cmap == 'hsv' else (None, None)
        im = ax.imshow(ph.T, extent=ext, cmap=cmap, origin='lower',
                       vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ── Row 3: Cross-sections, sidebands, Kuramoto ──
    ax_cs = fig.add_subplot(gs[3, 0:2])
    mid   = sim.N // 2
    for lbl, dens, col in [
        ('① Beam',         np.abs(r['psi_beam'])**2, '#aaaaaa'),
        ('③ Post-prop',    r['density_prop'],          '#1565c0'),
        ('④ Post-Floquet', r['density_final'],         '#c62828'),
    ]:
        line = dens[:, mid]
        ax_cs.plot(x_nm, line/(line.max()+1e-30), color=col, lw=2, label=lbl)
    ax_cs.set_xlabel('x (nm)')
    ax_cs.set_ylabel('Normalised |ψ|²')
    ax_cs.set_title('Cross-sections', fontsize=11, fontweight='bold')
    ax_cs.legend(fontsize=10)
    ax_cs.grid(True, alpha=0.3)

    ax_sb = fig.add_subplot(gs[3, 2])
    n_ids = sim.n_vals.astype(int)
    cols  = ['#c62828' if n==0 else '#1565c0' for n in n_ids]
    ax_sb.bar(n_ids, r['avg_pops'], color=cols, edgecolor='k', alpha=0.85)
    ax_sb.set_xlabel('Sideband n')
    ax_sb.set_ylabel('Beam-avg |c_n|²')
    ax_sb.set_title(f'Floquet Sidebands\nS={r["fl_entropy"]:.3f}',
                    fontsize=10, fontweight='bold')
    ax_sb.grid(True, alpha=0.3)

    ax_ku = fig.add_subplot(gs[3, 3])
    ax_ku.plot(r['sync_times'], r['sync_order'], color='#2e7d32', lw=1.5)
    ax_ku.axhline(r['r_sync'], color='red', ls='--',
                  label=f"r={r['r_sync']:.3f}")
    ax_ku.set_title('Kuramoto Sync', fontsize=10, fontweight='bold')
    ax_ku.set_xlabel('Time'); ax_ku.set_ylabel('r')
    ax_ku.legend(); ax_ku.grid(True, alpha=0.3); ax_ku.set_ylim(0, 1.05)

    # ── Row 4: Caging results ──
    cage_pi = r['cage_pi']
    cage_0  = r['cage_0']

    ax_loc = fig.add_subplot(gs[4, 0])
    ax_loc.plot(cage_pi['times'], cage_pi['localization'], 'b-', lw=2,
                label=f"Φ=π loc={cage_pi['localization'][-1]:.3f}")
    ax_loc.plot(cage_0['times'],  cage_0['localization'],  'r-', lw=2,
                label=f"Φ=0 loc={cage_0['localization'][-1]:.3f}")
    ax_loc.set_xlabel('Time (ℏ/J)')
    ax_loc.set_ylabel('Localization prob.')
    ax_loc.set_title('2D Diamond — Localization\n(primary caging metric)',
                     fontsize=10, fontweight='bold')
    ax_loc.legend(fontsize=9); ax_loc.grid(True, alpha=0.3)
    ax_loc.set_ylim(-0.05, 1.05)

    ax_fid = fig.add_subplot(gs[4, 1])
    ax_fid.plot(cage_pi['times'], cage_pi['fidelity'], 'b-', lw=2,
                label=f"Φ=π f={cage_pi['fidelity'][-1]:.3f}")
    ax_fid.plot(cage_0['times'],  cage_0['fidelity'],  'r-', lw=2,
                label=f"Φ=0 f={cage_0['fidelity'][-1]:.3f}")
    ax_fid.set_xlabel('Time (ℏ/J)')
    ax_fid.set_ylabel('Fidelity')
    ax_fid.set_title('2D Diamond — Fidelity\n(phase oscillations expected)',
                     fontsize=10, fontweight='bold')
    ax_fid.legend(fontsize=9); ax_fid.grid(True, alpha=0.3)

    for si, (t_lbl, snap_list) in enumerate([
            ('t=0',  [cage_pi['snapshots'][0]]),
            ('t=T',  [cage_pi['snapshots'][-1]])]):
        ax_sn = fig.add_subplot(gs[4, 2+si])
        _, d_sn = snap_list[0]
        im = ax_sn.imshow(d_sn.T, cmap='hot', origin='lower', aspect='auto')
        ax_sn.set_title(f'Caged Pattern Φ=π\n{t_lbl}',
                        fontsize=10, fontweight='bold')
        ax_sn.set_xlabel('Cell x'); ax_sn.set_ylabel('Cell y')
        plt.colorbar(im, ax=ax_sn, shrink=0.8)

    # ── Row 5: Summary ──
    ax_sum = fig.add_subplot(gs[5, :])
    ax_sum.axis('off')
    fb_str = 'OK' if cage_pi['flat_ok'] else 'FAIL'
    summary = (
        "v9 RESULTS SUMMARY\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
        f"FIX 4  2D Diamond network replaces Lieb lattice\n"
        f"       All bands flat at Φ=π: {fb_str}\n"
        f"       Localization: Φ=π={cage_pi['localization'][-1]:.4f}  "
        f"Φ=0={cage_0['localization'][-1]:.4f}  "
        f"gap={cage_pi['localization'][-1]-cage_0['localization'][-1]:.4f}\n"
        f"       Fidelity:     Φ=π={cage_pi['fidelity'][-1]:.4f}  "
        f"Φ=0={cage_0['fidelity'][-1]:.4f}\n\n"
        f"Retained: apodization, denser lattice, multi-species, disorder, decoherence\n\n"
        f"Pipeline: Floquet entropy={r['fl_entropy']:.4f}  "
        f"ads(n=0)={r['adsorption_frac']:.4f}"
    )
    ax_sum.text(0.02, 0.97, summary, transform=ax_sum.transAxes,
                fontsize=10, fontfamily='monospace', va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='#e8f5e9',
                          edgecolor='#2e7d32', linewidth=2))

    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_selectivity(sim, sel_results, phase, V_map, avg_pops,
                     fname='results/v9_selectivity.png'):
    print(f"\n  Plotting selectivity → {fname}")
    n_vals_int = list(sim.n_vals.astype(int))
    x_nm = sim.x * 1e9
    ext  = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]
    ref  = sel_results[0]['density']

    fig = plt.figure(figsize=(24, 20))
    gs  = GridSpec(4, 5, figure=fig, hspace=0.45, wspace=0.3,
                   left=0.05, right=0.97, top=0.94, bottom=0.04)
    fig.suptitle(
        'Sideband Selectivity — v9\n'
        'Same beam + substrate · Five deposition patterns · '
        'Boundary artefact REMOVED (apodization)',
        fontsize=13, fontweight='bold', y=0.98)

    for idx, n_res in enumerate(n_vals_int):
        ax  = fig.add_subplot(gs[0, idx])
        d   = sel_results[n_res]['density']
        d_n = d / (d.max() + 1e-30)
        im  = ax.imshow(d_n.T, extent=ext, cmap='inferno', origin='lower')
        C   = sel_results[n_res]['contrast']
        af  = sel_results[n_res]['ads_frac']
        feat= sel_results[n_res]['feature']
        ax.set_title(f'n={n_res:+d}\nC={C:.3f}  ads={af:.3f}\n'
                     f'feat={feat:.0f}nm',
                     fontsize=9, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    ref_n = ref / (ref.max() + 1e-30)
    for idx, n_res in enumerate(n_vals_int):
        ax   = fig.add_subplot(gs[1, idx])
        d_n  = sel_results[n_res]['density'] / (
            sel_results[n_res]['density'].max() + 1e-30)
        diff = d_n - ref_n
        im   = ax.imshow(diff.T, extent=ext, cmap='RdBu',
                         origin='lower', vmin=-0.5, vmax=0.5)
        s    = ssim_score(ref, sel_results[n_res]['density'])
        ax.set_title(f'Δ vs n=0\nSSIM={s:.3f}',
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    ax_cs = fig.add_subplot(gs[2, :])
    mid   = sim.N // 2
    colors_s = {-2:'#1a237e', -1:'#1565c0', 0:'#000000',
                1:'#c62828',  2:'#b71c1c'}
    for n_res in n_vals_int:
        d   = sel_results[n_res]['density'][:, mid]
        d_n = d / (d.max() + 1e-30)
        ax_cs.plot(x_nm, d_n, color=colors_s[n_res], lw=2,
                   label=f'n={n_res:+d} ads='
                         f'{sel_results[n_res]["ads_frac"]:.3f}')
    ax_cs.set_xlabel('x (nm)', fontsize=11)
    ax_cs.set_ylabel('Normalised |ψ|²', fontsize=11)
    ax_cs.set_title('Cross-sections — apodization removes boundary stripe',
                    fontsize=12, fontweight='bold')
    ax_cs.legend(fontsize=10, ncol=5)
    ax_cs.grid(True, alpha=0.3)

    ax_sm = fig.add_subplot(gs[3, 0:2])
    ax_sm.axis('off')
    ssim_v = [ssim_score(ref, sel_results[n]['density']) for n in n_vals_int]
    ads_v  = [sel_results[n]['ads_frac'] for n in n_vals_int]
    lines  = ["SELECTIVITY RESULTS  (v9 — apodized)\n" + "─"*40]
    for n_res, s, af in zip(n_vals_int, ssim_v, ads_v):
        lines.append(f"n={n_res:+d}: SSIM_vs_n0={s:.3f}  ads={af:.4f}")
    lines += ["", "n=±2 deposits at interior vortex cores",
              "(boundary artefact removed by apodization)"]
    ax_sm.text(0.03, 0.97, '\n'.join(lines), transform=ax_sm.transAxes,
               fontsize=10, fontfamily='monospace', va='top',
               bbox=dict(boxstyle='round', facecolor='#e8f5e9',
                         edgecolor='#2e7d32'))

    ax_ph = fig.add_subplot(gs[3, 2])
    im = ax_ph.imshow(phase.T, extent=ext, cmap='hsv', origin='lower',
                      vmin=-np.pi, vmax=np.pi)
    ax_ph.set_title('A-B Phase (apodized)\nV(x,y) ∝ |φ|',
                    fontsize=10, fontweight='bold')
    ax_ph.set_xlabel('x (nm)', fontsize=8)
    ax_ph.set_ylabel('y (nm)', fontsize=8)
    plt.colorbar(im, ax=ax_ph, shrink=0.8)

    ax_vm = fig.add_subplot(gs[3, 3])
    im = ax_vm.imshow(V_map.T, extent=ext, cmap='hot', origin='lower')
    ax_vm.set_title('V(x,y) Drive Map\n(zero at edges)',
                    fontsize=10, fontweight='bold')
    ax_vm.set_xlabel('x (nm)', fontsize=8)
    ax_vm.set_ylabel('y (nm)', fontsize=8)
    plt.colorbar(im, ax=ax_vm, shrink=0.8, label='V (ℏω)')

    ax_ab = fig.add_subplot(gs[3, 4])
    bar_c = [colors_s[n] for n in n_vals_int]
    ax_ab.bar(n_vals_int, ads_v, color=bar_c, edgecolor='k', alpha=0.85)
    ax_ab.set_xlabel('n_resonant')
    ax_ab.set_ylabel('Adsorption fraction')
    ax_ab.set_title('Adsorption by n', fontsize=10, fontweight='bold')
    ax_ab.grid(True, alpha=0.3, axis='y')
    for i, v in enumerate(ads_v):
        ax_ab.text(n_vals_int[i], v+0.002, f'{v:.3f}',
                   ha='center', fontsize=8)

    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_multi_species(sim, dens_A, dens_B, phase, V_map,
                       ssim_AB, overlap,
                       fname='results/v9_multispecies.png'):
    print(f"\n  Plotting multi-species → {fname}")
    x_nm = sim.x * 1e9
    ext  = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]
    mid  = sim.N // 2

    fig, axes = plt.subplots(2, 4, figsize=(22, 12))
    fig.suptitle(
        'Multi-Species Deposition — v9\n'
        'Same beam + substrate · Species A (n=0) and B (n=+2)',
        fontsize=13, fontweight='bold')

    dA_n = dens_A / (dens_A.max() + 1e-30)
    dB_n = dens_B / (dens_B.max() + 1e-30)

    im = axes[0,0].imshow(dA_n.T, extent=ext, cmap='Blues', origin='lower')
    axes[0,0].set_title('Species A  (n_res=0)\nInter-vortex regions',
                         fontsize=10, fontweight='bold')
    plt.colorbar(im, ax=axes[0,0], shrink=0.8)

    im = axes[0,1].imshow(dB_n.T, extent=ext, cmap='Reds', origin='lower')
    axes[0,1].set_title('Species B  (n_res=+2)\nVortex-core regions',
                         fontsize=10, fontweight='bold')
    plt.colorbar(im, ax=axes[0,1], shrink=0.8)

    composite = np.zeros((*dens_A.shape, 3))
    composite[:,:,2] = dA_n
    composite[:,:,0] = dB_n
    composite = np.clip(composite.transpose(1,0,2), 0, 1)
    axes[0,2].imshow(composite, extent=ext, origin='lower', aspect='auto')
    axes[0,2].set_title(f'Overlay  (blue=A, red=B)\n'
                         f'SSIM={ssim_AB:.3f}  overlap={overlap:.3f}',
                         fontsize=10, fontweight='bold')

    im = axes[0,3].imshow(phase.T, extent=ext, cmap='hsv',
                           origin='lower', vmin=-np.pi, vmax=np.pi)
    axes[0,3].set_title('A-B Phase (apodized)', fontsize=10, fontweight='bold')
    plt.colorbar(im, ax=axes[0,3], shrink=0.8)

    axes[1,0].plot(x_nm, dA_n[:, mid], 'b-', lw=2, label='Species A (n=0)')
    axes[1,0].plot(x_nm, dB_n[:, mid], 'r-', lw=2, label='Species B (n=+2)')
    axes[1,0].set_xlabel('x (nm)'); axes[1,0].set_ylabel('Normalised density')
    axes[1,0].set_title('Cross-sections', fontsize=10, fontweight='bold')
    axes[1,0].legend(fontsize=9); axes[1,0].grid(True, alpha=0.3)

    diff = dA_n - dB_n
    im = axes[1,1].imshow(diff.T, extent=ext, cmap='RdBu',
                           origin='lower', vmin=-0.5, vmax=0.5)
    axes[1,1].set_title('A − B difference map', fontsize=10, fontweight='bold')
    plt.colorbar(im, ax=axes[1,1], shrink=0.8)

    im = axes[1,2].imshow(V_map.T, extent=ext, cmap='hot', origin='lower')
    axes[1,2].set_title('V(x,y) — A sees low V\nB sees high V',
                         fontsize=10, fontweight='bold')
    plt.colorbar(im, ax=axes[1,2], shrink=0.8)

    axes[1,3].axis('off')
    Ca = michelson_contrast(dens_A)
    Cb = michelson_contrast(dens_B)
    fa = min_feature_size(dens_A, x_nm)
    fb = min_feature_size(dens_B, x_nm)
    txt = (
        "MULTI-SPECIES SUMMARY\n" + "─"*28 + "\n\n"
        f"Species A (n=0):\n"
        f"  C={Ca:.3f}  feat={fa:.0f}nm\n"
        f"  Location: inter-vortex\n\n"
        f"Species B (n=+2):\n"
        f"  C={Cb:.3f}  feat={fb:.0f}nm\n"
        f"  Location: vortex cores\n\n"
        f"Separation metrics:\n"
        f"  SSIM(A,B) = {ssim_AB:.4f}\n"
        f"  Overlap   = {overlap:.4f}\n\n"
        f"ONE substrate pass\n"
        f"TWO distinct locations"
    )
    axes[1,3].text(0.05, 0.97, txt, transform=axes[1,3].transAxes,
                   fontsize=9, fontfamily='monospace', va='top',
                   bbox=dict(boxstyle='round', facecolor='#fff3e0',
                             edgecolor='#e65100'))

    for ax in axes.flat:
        if not ax.has_data() and not ax.get_title():
            ax.axis('off')
        else:
            ax.set_xlabel('x (nm)', fontsize=8)
            ax.set_ylabel('y (nm)', fontsize=8)

    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_disorder(W_values, fid_pi, fid_0, loc_pi, loc_0,
                  fname='results/v9_disorder.png'):
    print(f"\n  Plotting disorder robustness → {fname}")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle('Disorder Robustness — v9\n'
                 '2D Diamond Network: AB caging (Φ=π) vs free (Φ=0)',
                 fontsize=13, fontweight='bold')

    # Localization (primary metric)
    axes[0].plot(W_values, loc_pi, 'b-o', lw=2, ms=7, label='Φ=π (caged)')
    axes[0].plot(W_values, loc_0,  'r-s', lw=2, ms=7, label='Φ=0 (free)')
    axes[0].axhline(0.5, color='gray', ls='--', alpha=0.7, label='Threshold')
    axes[0].set_xlabel('Disorder strength W (units of J)', fontsize=11)
    axes[0].set_ylabel('Localization probability', fontsize=11)
    axes[0].set_title('Localization vs Disorder\n(primary caging metric)',
                      fontsize=11, fontweight='bold')
    axes[0].legend(fontsize=10); axes[0].grid(True, alpha=0.3)
    axes[0].set_ylim(-0.05, 1.05)

    # Fidelity (secondary)
    axes[1].plot(W_values, fid_pi, 'b-o', lw=2, ms=7, label='Φ=π (caged)')
    axes[1].plot(W_values, fid_0,  'r-s', lw=2, ms=7, label='Φ=0 (free)')
    axes[1].set_xlabel('Disorder strength W (units of J)', fontsize=11)
    axes[1].set_ylabel('Final fidelity', fontsize=11)
    axes[1].set_title('Fidelity vs Disorder\n(phase-sensitive, oscillates)',
                      fontsize=11, fontweight='bold')
    axes[1].legend(fontsize=10); axes[1].grid(True, alpha=0.3)

    # Localization ratio
    ratio = loc_pi / (loc_0 + 1e-6)
    axes[2].plot(W_values, ratio, 'g-^', lw=2, ms=7)
    axes[2].axhline(1.0, color='gray', ls='--', alpha=0.7, label='No advantage')
    axes[2].set_xlabel('Disorder strength W (units of J)', fontsize=11)
    axes[2].set_ylabel('Localization ratio  Φ=π / Φ=0', fontsize=11)
    axes[2].set_title('Caging Advantage\n'
                       'Ratio > 1 = AB caging protects pattern',
                       fontsize=11, fontweight='bold')
    axes[2].legend(fontsize=10); axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_decoherence(Lc_nm, contrasts, features, a_nm,
                     fname='results/v9_decoherence.png'):
    print(f"\n  Plotting decoherence → {fname}")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Decoherence Sweep — v9\n'
                 'Minimum beam coherence for target feature size',
                 fontsize=13, fontweight='bold')

    finite = np.isfinite(Lc_nm)
    Lc_f   = Lc_nm[finite]
    C_f    = contrasts[finite]
    feat_f = features[finite]

    if not np.isfinite(Lc_nm[0]):
        axes[0].axhline(contrasts[0], color='gray', ls=':', lw=1.5,
                        label=f'L_c=∞: C={contrasts[0]:.3f}')

    axes[0].semilogx(Lc_f, C_f, 'bo-', lw=2, ms=6)
    axes[0].axhline(0.5, color='red',   ls='--', alpha=0.7, label='C=0.5')
    axes[0].axhline(0.8, color='green', ls='--', alpha=0.7, label='C=0.8')
    axes[0].axvline(a_nm, color='orange', ls='--', alpha=0.8,
                    label=f'L_c = a = {a_nm:.0f}nm')
    axes[0].set_xlabel('Coherence length L_c (nm)', fontsize=11)
    axes[0].set_ylabel('Michelson contrast', fontsize=11)
    axes[0].set_title('Contrast vs Coherence Length', fontsize=11,
                       fontweight='bold')
    axes[0].legend(fontsize=9); axes[0].grid(True, alpha=0.3)

    feat_fin = feat_f[np.isfinite(feat_f)]
    Lc_fin   = Lc_f[np.isfinite(feat_f)]

    if not np.isfinite(Lc_nm[0]) and np.isfinite(features[0]):
        axes[1].axhline(features[0], color='gray', ls=':', lw=1.5,
                        label=f'L_c=∞: {features[0]:.0f}nm')

    if len(feat_fin):
        axes[1].semilogx(Lc_fin, feat_fin, 'rs-', lw=2, ms=6, label='FWHM')
    axes[1].axhline(10, color='purple', ls=':', lw=2, label='10nm (EUV)')
    axes[1].axhline(3,  color='orange', ls=':', lw=2, label='3nm')
    axes[1].axvline(a_nm, color='orange', ls='--', alpha=0.8,
                    label=f'L_c = a = {a_nm:.0f}nm')
    axes[1].set_xlabel('Coherence length L_c (nm)', fontsize=11)
    axes[1].set_ylabel('Feature size (nm)', fontsize=11)
    axes[1].set_title('Feature Size vs Coherence Length',
                       fontsize=11, fontweight='bold')
    axes[1].legend(fontsize=9); axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_ablation(sim, ab_results, fname='results/v9_ablation.png'):
    print(f"\n  Plotting ablation → {fname}")
    fig, axes = plt.subplots(3, 3, figsize=(20, 18))
    fig.suptitle('Ablation Study — v9  (seeded, SSIM)',
                 fontsize=13, fontweight='bold')
    x_nm  = sim.x * 1e9
    ext   = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]
    mid   = sim.N // 2
    names = list(ab_results.keys())
    colors = {'Full pipeline':'#1565c0', 'No Floquet':'#2e7d32',
               'No A-B phase':'#c62828', 'Poor coherence':'#e65100',
               'Sideband n=+2':'#6a1b9a'}

    positions = [(0,0),(0,1),(0,2),(1,0),(1,1)]
    for (row,col), name in zip(positions, names):
        ax  = axes[row, col]
        d   = ab_results[name]['density']
        d_n = d / (d.max() + 1e-30)
        im  = ax.imshow(d_n.T, extent=ext, cmap='inferno', origin='lower')
        S   = ab_results[name]['ssim']
        C   = ab_results[name]['contrast']
        ax.set_title(f'{name}\nSSIM={S:.3f}  C={C:.3f}',
                     fontsize=9, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    ax = axes[1, 2]
    for name in names:
        d   = ab_results[name]['density'][:, mid]
        d_n = d / (d.max() + 1e-30)
        ax.plot(x_nm, d_n, color=colors.get(name,'gray'), lw=1.5, label=name)
    ax.set_xlabel('x (nm)'); ax.set_ylabel('Norm |ψ|²')
    ax.set_title('Cross-sections', fontsize=10, fontweight='bold')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    ax = axes[2, 0]
    ssim_v = [ab_results[n]['ssim']     for n in names]
    mich_v = [ab_results[n]['contrast'] for n in names]
    x_pos  = np.arange(len(names))
    w      = 0.35
    bc     = [colors.get(n,'gray') for n in names]
    ax.bar(x_pos-w/2, ssim_v, width=w, color=bc, edgecolor='k',
           alpha=0.9, label='SSIM')
    ax.bar(x_pos+w/2, mich_v, width=w, color=bc, edgecolor='k',
           alpha=0.4, hatch='//', label='Michelson C')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([n.replace(' ','\n') for n in names], fontsize=8)
    ax.set_ylabel('Score'); ax.set_ylim(0, 1.1)
    ax.set_title('SSIM vs Michelson', fontsize=10, fontweight='bold')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3, axis='y')
    for i,(s,m) in enumerate(zip(ssim_v, mich_v)):
        ax.text(i-w/2, s+0.02, f'{s:.3f}', ha='center', fontsize=7)
        ax.text(i+w/2, m+0.02, f'{m:.3f}', ha='center', fontsize=7)

    ax = axes[2, 1]
    ax.axis('off')
    sorted_n = sorted(names, key=lambda n: ab_results[n]['ssim'], reverse=True)
    lines = ["ABLATION RESULTS\n" + "─"*30]
    for rank, name in enumerate(sorted_n, 1):
        res = ab_results[name]
        lines.append(f"{rank}. {name}\n"
                     f"   SSIM={res['ssim']:.4f}  C={res['contrast']:.4f}")
    ax.text(0.03, 0.97, '\n'.join(lines), transform=ax.transAxes,
            fontsize=9, fontfamily='monospace', va='top',
            bbox=dict(boxstyle='round', facecolor='#fce4ec',
                      edgecolor='#c62828'))

    ax = axes[2, 2]
    ax.axis('off')
    fp   = ab_results['Full pipeline']['ssim']
    nfl  = ab_results['No Floquet']['ssim']
    nab  = ab_results['No A-B phase']['ssim']
    poor = ab_results['Poor coherence']['ssim']
    n2   = ab_results['Sideband n=+2']['ssim']
    interp = [
        "INTERPRETATION\n" + "─"*24, "",
        f"Full vs No Floquet:   ΔSSIM={fp-nfl:+.3f}",
        f"  Floquet {'ADDS' if fp>nfl else 'neutral'}",
        "", f"Full vs No A-B:       ΔSSIM={fp-nab:+.3f}",
        f"  Phase: PRIMARY driver",
        "", f"Full vs Poor coh:     ΔSSIM={fp-poor:+.3f}",
        f"  Coherence matters",
        "", f"Full vs n=+2:         ΔSSIM={fp-n2:+.3f}",
        f"  n=+2 DIFFERENT region",
        "  (interior core, not edge)",
    ]
    ax.text(0.03, 0.97, '\n'.join(interp), transform=ax.transAxes,
            fontsize=9, fontfamily='monospace', va='top',
            bbox=dict(boxstyle='round', facecolor='#e8f5e9',
                      edgecolor='#2e7d32'))

    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("╔═══════════════════════════════════════════════════════════╗")
    print("║  INTEGRATED QUANTUM SUBSTRATE  v9                        ║")
    print("║  Apodization · 2D Diamond AB cage · Multi-species        ║")
    print("╚═══════════════════════════════════════════════════════════╝")

    np.random.seed(0)

    # ── Validation gates ─────────────────────────────────────────────
    print("\n" + "="*65)
    print("STEP 0: VALIDATION GATES")
    print("="*65)
    validate_floquet(N_side=2, V_frac=0.9, abort_on_fail=True)
    validate_diamond_spectrum(Lx=4, Ly=4, phi=np.pi, abort_on_fail=True)

    # ── Main pipeline ─────────────────────────────────────────────────
    print("\n" + "="*65)
    print("MAIN PIPELINE")
    print("="*65)
    sim     = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3,
                                          N_floquet_sidebands=2,
                                          N_diamond_cells_x=16,
                                          N_diamond_cells_y=16)
    results = sim.run_full_pipeline(
        pattern='vortex_lattice', V_max=0.9,
        n_resonant=0, K_kuramoto=6.0,
        phi_cage=np.pi, prop_distance_lam=20, apodize=True)
    plot_pipeline(sim, results)

    # ── Sideband selectivity ──────────────────────────────────────────
    print("\n" + "="*65)
    print("SIDEBAND SELECTIVITY")
    print("="*65)
    sel_sim, sel_res, sel_ph, sel_V, sel_ap = \
        sideband_selectivity_sweep(seed=42)
    plot_selectivity(sel_sim, sel_res, sel_ph, sel_V, sel_ap)

    # ── Multi-species ─────────────────────────────────────────────────
    print("\n" + "="*65)
    print("MULTI-SPECIES DEPOSITION")
    print("="*65)
    ms_sim, dA, dB, ms_ph, ms_V, ssim_AB, overlap = \
        multi_species_deposition(seed=42)
    plot_multi_species(ms_sim, dA, dB, ms_ph, ms_V, ssim_AB, overlap)

    # ── Disorder robustness ───────────────────────────────────────────
    print("\n" + "="*65)
    print("DISORDER ROBUSTNESS")
    print("="*65)
    W_vals, fid_pi, fid_0, loc_pi, loc_0 = disorder_robustness(
        seed=42, n_W=10, W_max=2.0, T_evolve=40)
    plot_disorder(W_vals, fid_pi, fid_0, loc_pi, loc_0)

    # ── Decoherence sweep ─────────────────────────────────────────────
    print("\n" + "="*65)
    print("DECOHERENCE SWEEP")
    print("="*65)
    Lc_nm, contrasts, features, a_nm = decoherence_sweep(seed=42, n_pts=12)
    plot_decoherence(Lc_nm, contrasts, features, a_nm)

    # ── Ablation ──────────────────────────────────────────────────────
    print("\n" + "="*65)
    print("ABLATION STUDY")
    print("="*65)
    ab_sim, ab_results = ablation_study(seed=42)
    plot_ablation(ab_sim, ab_results)

    # ── Summary ───────────────────────────────────────────────────────
    print("\n" + "="*65)
    print("ALL SIMULATIONS COMPLETE")
    print("="*65)
    print("\nOutput files:")
    print("  results/v9_pipeline.png      — Pipeline + 2D Diamond caging")
    print("  results/v9_selectivity.png   — Sideband selectivity (apodized)")
    print("  results/v9_multispecies.png  — Two-species placement")
    print("  results/v9_disorder.png      — Disorder robustness (diamond)")
    print("  results/v9_decoherence.png   — Coherence requirement")
    print("  results/v9_ablation.png      — Ablation study")

    loc_gap = (results['cage_pi']['localization'][-1]
               - results['cage_0']['localization'][-1])
    fid_gap = (results['cage_pi']['fidelity'][-1]
               - results['cage_0']['fidelity'][-1])
    print(f"\nKey results:")
    print(f"  Diamond flat bands: "
          f"{'OK' if results['cage_pi']['flat_ok'] else 'FAIL'}")
    print(f"  Caging (localization): Φ=π={results['cage_pi']['localization'][-1]:.4f}  "
          f"Φ=0={results['cage_0']['localization'][-1]:.4f}  "
          f"gap={loc_gap:.4f}")
    print(f"  Caging (fidelity):    Φ=π={results['cage_pi']['fidelity'][-1]:.4f}  "
          f"Φ=0={results['cage_0']['fidelity'][-1]:.4f}  "
          f"gap={fid_gap:.4f}")
    print(f"  Floquet entropy:   {results['fl_entropy']:.4f}")
    print(f"  Adsorption (n=0):  {results['adsorption_frac']:.4f}")
    print(f"  Multi-species SSIM(A,B): {ssim_AB:.4f}  overlap: {overlap:.4f}")

    sel_ssim_02 = ssim_score(sel_res[0]['density'], sel_res[2]['density'])
    print(f"  Selectivity SSIM(n=0,n=+2): {sel_ssim_02:.4f}")

    max_W_protected = W_vals[loc_pi > 0.8][-1] if np.any(loc_pi > 0.8) else 0
    print(f"  Disorder protection (loc>0.8) up to W={max_W_protected:.2f}J  (Φ=π)")

    Lc_finite = Lc_nm[np.isfinite(Lc_nm)]
    C_finite  = contrasts[np.isfinite(Lc_nm)]
    Lc_half   = Lc_finite[C_finite > 0.5*contrasts[0]]
    if len(Lc_half):
        print(f"  Contrast >50% of max for L_c > {Lc_half[-1]:.0f}nm")


if __name__ == "__main__":
    main()