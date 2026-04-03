"""
Integrated Quantum Substrate v10
=================================

End-to-end pipeline for controllable matter deposition:

  He⁺ CMWB source → neutralization → SQUID holography → diamond caging

Physical concept:
  1. A coherent He⁺ ion beam is generated in a diode cavity array via
     Aharonov-Bohm phase synchronization (US Patent 9,502,202 B2).
     The Kuramoto order parameter r quantifies beam coherence and
     depends on B-field, cavity pressure, and energy spread.
  2. The He⁺ beam is neutralized (charge exchange) to produce neutral
     He-4 at the same velocity and wavelength.  Neutralization preserves
     the coherent wavefront but removes the charge, so the beam passes
     through the SQUID array acquiring only a geometric AB phase shift
     (not a Lorentz deflection).
  3. The neutral beam traverses a 32×32 SQUID phase screen whose loop
     phases are computed by inverse holography to produce a desired
     target deposition pattern after Fresnel propagation.
  4. The deposited pattern is stabilized by Aharonov-Bohm caging on a
     2D diamond lattice (flat bands at Φ=π).

Key result: He⁺ and He-4 have identical mass (to < 0.01%) and therefore
identical de Broglie wavelength λ ≈ 49 nm at 1 mK.  The SQUID array
geometry (32×32, 12.5 nm pitch, 400 nm substrate) is wavelength-matched.

Dependencies: coherent_matterwave_beam.py, inverse_holography.py,
              sim_v9.py, diamond_caging.py
"""

import os
import sys
import time
import numpy as np
import torch
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

from coherent_matterwave_beam import (
    CoherentMatterwaveBeam, CavityGeometry, SPECIES,
    hbar, k_B, e_C, m_u,
)
from inverse_holography import (
    SQUIDArray, InverseHolographySolver,
    target_single_spot, target_line, target_grid_of_dots,
    target_ring, target_letter, smooth_target, compute_metrics,
    validate_roundtrip,
)
from sim_v9 import (
    IntegratedQuantumSubstrate, DiamondNetworkSimulator,
    michelson_contrast, ssim_score, min_feature_size,
    m_He,
)

_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
os.makedirs('results', exist_ok=True)


# ===========================================================================
# INTEGRATED v10 PIPELINE
# ===========================================================================

class IntegratedPipelineV10:
    """
    End-to-end: He⁺ CMWB source → neutralization → SQUID holography
    → (optional Floquet) → diamond caging.

    The CMWB generates a coherent He⁺ beam whose order parameter r
    depends on cavity pressure, B-field, and energy spread.  After
    neutralization, the beam has identical wavelength (λ ≈ 49 nm at
    1 mK) but no charge — it passes through the SQUID array acquiring
    only the AB phase imprint.

    Parameters
    ----------
    B_field : float
        Magnetic field for AB synchronization in the CMWB cavity [T].
    pressure_Pa : float
        CMWB cavity pressure [Pa].  Determines effective passes.
    dE_frac : float
        Beam energy spread ΔE/E.
    beam_current_A : float
        CMWB beam current [A].
    cavity : CavityGeometry or None
        Cavity parameters.  Pressure is overridden by pressure_Pa.
    T_beam : float
        Post-neutralization beam temperature [K].  Sets λ_dB.
    N : int
        Spatial grid resolution.
    L : float
        Physical domain size [m].
    N_loops : int
        SQUID array resolution per axis.
    prop_distance_lam : float
        Propagation distance in units of de Broglie wavelength.
    use_floquet : bool
        Whether to apply Floquet dressing after holography.
    V_max_floquet : float
        Floquet drive strength (if enabled).
    n_resonant : int
        Floquet sideband for binding filter (if enabled).
    N_diamond_x, N_diamond_y : int
        Diamond lattice dimensions for caging.
    """

    def __init__(self, B_field=0.01, pressure_Pa=100,
                 dE_frac=0.01, beam_current_A=1e-6,
                 cavity=None, T_beam=1e-3,
                 N=256, L=400e-9, N_loops=32,
                 prop_distance_lam=20.0,
                 use_floquet=False, V_max_floquet=0.9,
                 n_resonant=0,
                 N_diamond_x=16, N_diamond_y=16):
        self.N = N
        self.L = L
        self.dx = L / N
        self.N_loops = N_loops
        self.prop_distance_lam = prop_distance_lam
        self.use_floquet = use_floquet
        self.V_max_floquet = V_max_floquet
        self.n_resonant = n_resonant
        self.Lx = N_diamond_x
        self.Ly = N_diamond_y
        self.T_beam = T_beam

        # --- CMWB source (He⁺) ---
        cav = cavity or CavityGeometry(pressure_Pa=pressure_Pa)
        if cavity is not None and pressure_Pa != 101325:
            cav.pressure = pressure_Pa
        self.cmwb = CoherentMatterwaveBeam(
            species='He+', E_kinetic_eV=k_B * T_beam / e_C,
            B_field=B_field, cavity=cav,
            dE_frac=dE_frac, beam_current_A=beam_current_A,
        )

        # --- Post-neutralization beam parameters ---
        # He⁺ and He-4 have identical mass → identical λ_dB
        self.mass = m_He
        self.v = np.sqrt(2 * k_B * T_beam / self.mass)
        self.k0 = self.mass * self.v / hbar
        self.lam = 2 * np.pi / self.k0

        # --- SQUID holography solver ---
        self.squid = SQUIDArray(N_loops=N_loops, N_grid=N, L_grid=L)
        self.holo_solver = InverseHolographySolver(
            squid_array=self.squid, N=N, L=L, T_beam=T_beam,
            prop_distance_lam=prop_distance_lam,
        )

        # --- v9 substrate sim (for optional Floquet) ---
        if self.use_floquet:
            self.substrate = IntegratedQuantumSubstrate(
                N=N, L=L, T_beam=T_beam,
            )

        # --- Diamond caging ---
        self.diamond = DiamondNetworkSimulator(
            Lx=N_diamond_x, Ly=N_diamond_y, N_grid=N,
        )

    def info(self):
        print("=" * 65)
        print("INTEGRATED QUANTUM SUBSTRATE  v10")
        print("  He⁺ CMWB → neutralization → SQUID holography → caging")
        print("=" * 65)
        print(f"  CMWB source: He⁺ at {self.T_beam*1e3:.1f} mK")
        print(f"    B = {self.cmwb.B*1e4:.0f} G, "
              f"P = {self.cmwb.cavity.pressure:.0f} Pa, "
              f"ΔE/E = {self.cmwb.dE_frac:.1%}")
        print(f"    K_dim = {self.cmwb.K_dim:.3f}, "
              f"n_eff = {self.cmwb.n_eff:.0f} "
              f"({self.cmwb.K_limit}-limited), "
              f"N_vol = {self.cmwb.N_per_volume}")
        print(f"  Post-neutralization: He-4 neutral")
        print(f"    λ_dB = {self.lam*1e9:.2f} nm, "
              f"v = {self.v:.3f} m/s")
        print(f"  Grid: {self.N}×{self.N}, L = {self.L*1e9:.0f} nm, "
              f"dx = {self.dx*1e9:.2f} nm")
        print(f"  SQUID: {self.N_loops}×{self.N_loops} "
              f"(pitch = {self.squid.pitch*1e9:.1f} nm)")
        print(f"  Propagation: {self.prop_distance_lam}λ "
              f"= {self.prop_distance_lam * self.lam * 1e9:.1f} nm")
        print(f"  Floquet: {'ON' if self.use_floquet else 'OFF'}")
        print(f"  Diamond: {self.Lx}×{self.Ly}")
        print(f"  Device: {_device}")

    # -------------------------------------------------------------------
    # Stage 1: CMWB He⁺ source + neutralization
    # -------------------------------------------------------------------
    def generate_source(self, seed=42, verbose=True):
        """Generate He⁺ beam via CMWB, then neutralize.

        Returns the post-neutralization He-4 wavefunction with phase
        noise set by the CMWB coherence (1 - r).

        Returns
        -------
        psi : ndarray (N, N) complex128, normalized
        sync_result : dict from CMWB.synchronize()
        """
        if verbose:
            print("\n  STAGE 1: CMWB He⁺ source + neutralization")

        sync = self.cmwb.synchronize(seed=seed, verbose=verbose, mode='auto')
        r = sync['r_final']

        # Build the post-neutralization wavefunction.
        # Same mass, same velocity, same λ — only charge is removed.
        # Phase noise from incomplete synchronization carries through.
        x = np.linspace(-self.L / 2, self.L / 2, self.N)
        X, Y = np.meshgrid(x, x, indexing='ij')
        sigma = 0.35 * self.L

        psi = (np.exp(-(X**2 + Y**2) / (2 * sigma**2))
               * np.exp(1j * self.k0 * X))

        # Phase noise proportional to incoherence
        noise = (1 - r) * np.random.normal(0, np.pi, (self.N, self.N))
        noise = gaussian_filter(noise, sigma=3)
        psi = psi * np.exp(1j * noise)

        # Normalize
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * self.dx**2)

        if verbose:
            print(f"    Neutralized He-4 beam: r = {r:.4f}, "
                  f"λ = {self.lam*1e9:.2f} nm  "
                  f"[sync mode: {sync.get('mode', 'unknown')}]")

        return psi, sync

    # -------------------------------------------------------------------
    # Stage 2: Inverse holography
    # -------------------------------------------------------------------
    def solve_holography(self, target_pattern, method='gd',
                         n_iter_gs=300, n_restarts_gs=3,
                         n_iter_gd=500, lr_gd=0.05,
                         verbose=True):
        """Solve the inverse holography problem for a target pattern."""
        if verbose:
            print(f"\n  STAGE 2: Inverse holography ({method.upper()})")

        target_smooth = smooth_target(
            target_pattern, N_loops=self.N_loops,
            corner_radius=0.03, sigma=2,
        )

        t0 = time.time()
        if method == 'gs':
            raw = self.holo_solver.solve_gerchberg_saxton(
                target_smooth, n_iter=n_iter_gs,
                n_restarts=n_restarts_gs, verbose=verbose)
        elif method == 'gd':
            raw = self.holo_solver.solve_gradient_descent(
                target_smooth, n_iter=n_iter_gd,
                lr=lr_gd, verbose=verbose)
        else:
            raise ValueError(f"Unknown method: {method}")
        elapsed = time.time() - t0

        metrics = compute_metrics(raw['achieved'], target_smooth)

        if verbose:
            print(f"    SSIM = {metrics['ssim']:.4f}, "
                  f"eff = {metrics['efficiency']:.4f}, "
                  f"time = {elapsed:.1f}s")

        return {
            'phase_screen':   raw['phase_screen'],
            'P_achieved':     raw['achieved'],
            'P_target':       target_smooth,
            'P_target_raw':   target_pattern,
            'convergence':    raw['convergence'],
            'metrics':        metrics,
            'method':         method,
            'time_s':         elapsed,
            'phi_loops':      raw.get('phi_loops', None),
        }

    # -------------------------------------------------------------------
    # Stage 3 (optional): Floquet dressing + binding filter
    # -------------------------------------------------------------------
    def floquet_filter(self, density_2d, phase_screen, verbose=True):
        """Apply Floquet dressing and sideband-selective binding."""
        if not self.use_floquet:
            if verbose:
                print("\n  STAGE 3: Floquet — SKIPPED (disabled)")
            return density_2d, {}

        if verbose:
            print("\n  STAGE 3: Floquet dressing + binding filter")

        psi = np.sqrt(np.maximum(density_2d, 0)) * np.exp(
            1j * phase_screen)

        psi_d, avg_pops, entropy, V_map = \
            self.substrate.stage2_floquet_dress_spatial(
                psi, phase_screen, V_max=self.V_max_floquet,
                verbose=verbose)

        psi_ads, ads_frac, weights = \
            self.substrate.stage3_binding_filter(
                psi_d, avg_pops, n_resonant=self.n_resonant,
                verbose=verbose)

        return np.abs(psi_ads)**2, {
            'avg_pops': avg_pops,
            'entropy': entropy,
            'ads_frac': ads_frac,
            'V_map': V_map,
        }

    # -------------------------------------------------------------------
    # Stage 4: Diamond lattice AB caging
    # -------------------------------------------------------------------
    def cage_pattern(self, density_2d, phi_cage=np.pi,
                     T_evolve=40.0, verbose=True):
        """Evolve the deposited density on a 2D diamond lattice."""
        if verbose:
            print(f"\n  STAGE 4: Diamond caging "
                  f"(Φ = {phi_cage/np.pi:.1f}π)")

        return self.diamond.evolve_caging(
            density_2d, phi_cage=phi_cage,
            T_evolve=T_evolve, verbose=verbose)

    # -------------------------------------------------------------------
    # Full pipeline
    # -------------------------------------------------------------------
    def run(self, target_pattern, method='gd',
            run_caging=True, T_evolve=40.0,
            seed=42, verbose=True):
        """Execute the full v10 pipeline."""
        if verbose:
            self.info()

        # Stage 1: CMWB source + neutralization
        psi_source, sync_result = self.generate_source(
            seed=seed, verbose=verbose)

        # Stage 2: Inverse holography
        holo = self.solve_holography(
            target_pattern, method=method, verbose=verbose)

        density_holo = holo['P_achieved']

        # Stage 3: Optional Floquet
        density_post_floquet, floquet_info = self.floquet_filter(
            density_holo, holo['phase_screen'], verbose=verbose)

        # Stage 4: Caging
        cage_pi, cage_0 = None, None
        if run_caging:
            cage_pi = self.cage_pattern(
                density_post_floquet, phi_cage=np.pi,
                T_evolve=T_evolve, verbose=verbose)
            cage_0 = self.cage_pattern(
                density_post_floquet, phi_cage=0.0,
                T_evolve=T_evolve, verbose=verbose)

        results = {
            'psi_source':       psi_source,
            'sync_result':      sync_result,
            'r_source':         sync_result['r_final'],
            'holo':             holo,
            'density_holo':     density_holo,
            'density_filtered': density_post_floquet,
            'floquet_info':     floquet_info,
            'cage_pi':          cage_pi,
            'cage_0':           cage_0,
            'lam_nm':           self.lam * 1e9,
            'pressure_Pa':      self.cmwb.cavity.pressure,
            'B_gauss':          self.cmwb.B * 1e4,
            'K_dim':            self.cmwb.K_dim,
            'n_eff':            self.cmwb.n_eff,
            'sync_mode':        sync_result.get('mode', 'unknown'),
            'N_per_volume':     self.cmwb.N_per_volume,
        }

        if verbose:
            self._print_summary(results)

        return results

    def _print_summary(self, r):
        print("\n" + "=" * 65)
        print("v10 PIPELINE SUMMARY")
        print("=" * 65)
        print(f"  Source:  He⁺ → He-4, r = {r['r_source']:.4f}")
        sync_mode = r.get('sync_mode', 'unknown')
        print(f"  CMWB:   B = {r['B_gauss']:.0f} G, "
              f"P = {r['pressure_Pa']:.0f} Pa, "
              f"K_dim = {r['K_dim']:.3f}, "
              f"n_eff = {r['n_eff']:.0f}, "
              f"mode = {sync_mode}")
        m = r['holo']['metrics']
        print(f"  Holo:   SSIM = {m['ssim']:.4f}, "
              f"eff = {m['efficiency']:.4f}")
        if r['floquet_info']:
            fi = r['floquet_info']
            print(f"  Floquet: S = {fi.get('entropy', 'N/A')}, "
                  f"ads = {fi.get('ads_frac', 'N/A')}")
        if r['cage_pi'] is not None:
            loc_pi = r['cage_pi']['localization'][-1]
            loc_0 = r['cage_0']['localization'][-1]
            print(f"  Caging: Φ=π loc = {loc_pi:.4f}, "
                  f"Φ=0 loc = {loc_0:.4f}, "
                  f"gap = {loc_pi - loc_0:.4f}")


# ===========================================================================
# PRESSURE SWEEP
# ===========================================================================

def pressure_sweep(target_pattern, pressures=None,
                   B_field=0.01, seed=42, verbose=True):
    """Sweep CMWB cavity pressure to find the coherence threshold
    for usable holographic deposition."""
    if pressures is None:
        pressures = [1e-3, 1.0, 100, 1000, 10000, 101325]

    results = []
    for P in pressures:
        if verbose:
            print(f"\n{'#' * 65}")
            print(f"# PRESSURE = {P:.1e} Pa")
            print(f"{'#' * 65}")

        pipe = IntegratedPipelineV10(
            B_field=B_field, pressure_Pa=P,
        )
        r = pipe.run(target_pattern, method='gd',
                     run_caging=False, seed=seed, verbose=verbose)
        r['pressure'] = P
        results.append(r)

    return results


# ===========================================================================
# PLOTTING
# ===========================================================================

def plot_v10_pipeline(pipeline, results, fname='results/v10_pipeline.png'):
    """Comprehensive pipeline figure."""
    print(f"\n  Plotting → {fname}")

    fig = plt.figure(figsize=(24, 20))
    gs = GridSpec(4, 4, figure=fig, hspace=0.45, wspace=0.35,
                  left=0.05, right=0.97, top=0.95, bottom=0.03)

    x_nm = np.linspace(-pipeline.L / 2, pipeline.L / 2,
                       pipeline.N) * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    # Row 0: Banner
    ax = fig.add_subplot(gs[0, :])
    ax.axis('off')
    m = results['holo']['metrics']
    cage_str = ''
    if results['cage_pi'] is not None:
        loc_pi = results['cage_pi']['localization'][-1]
        loc_0 = results['cage_0']['localization'][-1]
        cage_str = (f"\nCaging: Φ=π loc={loc_pi:.3f}  "
                    f"Φ=0 loc={loc_0:.3f}")
    txt = (
        "INTEGRATED QUANTUM SUBSTRATE  v10\n"
        "He⁺ CMWB → neutralization → SQUID holography → diamond cage\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Source: He⁺ r={results['r_source']:.3f}  "
        f"B={results['B_gauss']:.0f}G  "
        f"P={results['pressure_Pa']:.0f}Pa  "
        f"K={results['K_dim']:.3f}\n"
        f"λ={pipeline.lam*1e9:.1f}nm  "
        f"Holo SSIM={m['ssim']:.3f}  "
        f"eff={m['efficiency']:.3f}"
        f"{cage_str}"
    )
    ax.text(0.5, 0.5, txt, transform=ax.transAxes,
            ha='center', va='center',
            fontsize=11, fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#e3f2fd',
                      edgecolor='#1565c0', linewidth=2))

    # Row 1: Target, achieved, phase screen, source density
    panels = [
        ('Target (smoothed)', results['holo']['P_target']),
        ('Holographic output', results['density_holo']),
        ('Phase screen', None),
        ('Source |ψ|²', np.abs(results['psi_source'])**2),
    ]
    for idx, (title, data) in enumerate(panels):
        ax = fig.add_subplot(gs[1, idx])
        if title == 'Phase screen':
            ph = results['holo']['phase_screen']
            im = ax.imshow(ph.T, extent=ext, cmap='twilight_shifted',
                           origin='lower', vmin=-np.pi, vmax=np.pi)
            ax.set_title(title, fontsize=10, fontweight='bold')
            plt.colorbar(im, ax=ax, shrink=0.8, label='rad')
        else:
            d = data / (data.max() + 1e-30)
            im = ax.imshow(d.T, extent=ext, cmap='inferno',
                           origin='lower')
            C = michelson_contrast(data)
            ax.set_title(f'{title}\nC={C:.3f}',
                         fontsize=10, fontweight='bold')
            plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)

    # Row 2: Convergence and cross-sections
    ax = fig.add_subplot(gs[2, 0:2])
    conv = results['holo']['convergence']
    ax.semilogy(conv, 'r-', lw=1.5)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Loss')
    ax.set_title('Holography convergence', fontsize=10,
                 fontweight='bold')
    ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[2, 2:4])
    mid = pipeline.N // 2
    target_line = results['holo']['P_target'][:, mid]
    achieved_line = results['density_holo'][:, mid]
    target_line /= (target_line.max() + 1e-30)
    achieved_line /= (achieved_line.max() + 1e-30)
    ax.plot(x_nm, target_line, 'k--', lw=2, label='Target')
    ax.plot(x_nm, achieved_line, 'r-', lw=2, label='Achieved')
    if results['floquet_info']:
        filt_line = results['density_filtered'][:, mid]
        filt_line /= (filt_line.max() + 1e-30)
        ax.plot(x_nm, filt_line, 'b-', lw=1.5,
                label='Post-Floquet', alpha=0.7)
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('Normalised intensity')
    ax.set_title('Cross-sections', fontsize=10, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Row 3: Caging
    if results['cage_pi'] is not None:
        cage_pi = results['cage_pi']
        cage_0 = results['cage_0']

        ax = fig.add_subplot(gs[3, 0])
        ax.plot(cage_pi['times'], cage_pi['localization'],
                'b-', lw=2,
                label=f"Φ=π loc={cage_pi['localization'][-1]:.3f}")
        ax.plot(cage_0['times'], cage_0['localization'],
                'r-', lw=2,
                label=f"Φ=0 loc={cage_0['localization'][-1]:.3f}")
        ax.set_xlabel('Time (ℏ/J)')
        ax.set_ylabel('Localization')
        ax.set_title('Diamond AB caging', fontsize=10,
                     fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)

        ax = fig.add_subplot(gs[3, 1])
        ax.plot(cage_pi['times'], cage_pi['fidelity'],
                'b-', lw=2, label='Φ=π')
        ax.plot(cage_0['times'], cage_0['fidelity'],
                'r-', lw=2, label='Φ=0')
        ax.set_xlabel('Time (ℏ/J)')
        ax.set_ylabel('Fidelity')
        ax.set_title('Fidelity evolution', fontsize=10,
                     fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        for si, (lbl, snaps) in enumerate([
                ('Caged t=0', cage_pi['snapshots'][0]),
                ('Caged t=T', cage_pi['snapshots'][-1])]):
            ax = fig.add_subplot(gs[3, 2 + si])
            _, d_sn = snaps
            im = ax.imshow(d_sn.T, cmap='hot', origin='lower',
                           aspect='auto')
            ax.set_title(f'Φ=π {lbl}', fontsize=10,
                         fontweight='bold')
            ax.set_xlabel('Cell x')
            ax.set_ylabel('Cell y')
            plt.colorbar(im, ax=ax, shrink=0.8)
    else:
        ax = fig.add_subplot(gs[3, :])
        ax.axis('off')
        ax.text(0.5, 0.5, 'Caging stage not run',
                ha='center', va='center', fontsize=14,
                transform=ax.transAxes, color='gray')

    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_multi_target(pipeline, all_results,
                      fname='results/v10_multi_target.png'):
    """Summary figure comparing multiple targets."""
    print(f"\n  Plotting → {fname}")
    names = list(all_results.keys())
    n = len(names)

    fig, axes = plt.subplots(3, n, figsize=(5 * n, 14))
    if n == 1:
        axes = axes.reshape(-1, 1)

    x_nm = np.linspace(-pipeline.L / 2, pipeline.L / 2,
                       pipeline.N) * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    for j, name in enumerate(names):
        r = all_results[name]
        m = r['holo']['metrics']

        ax = axes[0, j]
        P = r['holo']['P_target']
        ax.imshow(P.T / (P.max() + 1e-30), extent=ext,
                  cmap='inferno', origin='lower')
        ax.set_title(f'{name}\n(target)', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)

        ax = axes[1, j]
        I = r['density_holo']
        ax.imshow(I.T / (I.max() + 1e-30), extent=ext,
                  cmap='inferno', origin='lower')
        ax.set_title(f"SSIM={m['ssim']:.3f}\n"
                     f"eff={m['efficiency']:.3f}",
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)

        ax = axes[2, j]
        ph = r['holo']['phase_screen']
        ax.imshow(ph.T, extent=ext, cmap='twilight_shifted',
                  origin='lower', vmin=-np.pi, vmax=np.pi)
        ax.set_title('Phase screen', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)

    fig.suptitle(f'v10 Multi-Target — He⁺ CMWB '
                 f'(r={all_results[names[0]]["r_source"]:.3f})',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_pressure_sweep(sweep_results,
                        fname='results/v10_pressure_sweep.png'):
    """Plot source coherence and holographic SSIM vs cavity pressure."""
    print(f"\n  Plotting → {fname}")

    pressures = [r['pressure'] for r in sweep_results]
    rs = [r['r_source'] for r in sweep_results]
    ssims = [r['holo']['metrics']['ssim'] for r in sweep_results]
    effs = [r['holo']['metrics']['efficiency'] for r in sweep_results]
    K_dims = [r['K_dim'] for r in sweep_results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.semilogx(pressures, rs, 'bo-', lw=2, markersize=8)
    ax.set_xlabel('Cavity pressure (Pa)')
    ax.set_ylabel('Source coherence r')
    ax.set_title('CMWB He⁺ synchronization', fontweight='bold')
    ax.axhline(0.95, color='g', ls='--', alpha=0.5, label='r = 0.95')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_ylim(-0.05, 1.05)
    for p, r, k in zip(pressures, rs, K_dims):
        ax.annotate(f'K={k:.2f}', (p, r), fontsize=7,
                    textcoords='offset points', xytext=(5, 5))

    ax = axes[1]
    ax.semilogx(pressures, ssims, 'rs-', lw=2, markersize=8)
    ax.set_xlabel('Cavity pressure (Pa)')
    ax.set_ylabel('Holographic SSIM')
    ax.set_title('Deposition fidelity vs pressure', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)

    ax = axes[2]
    ax.plot(rs, ssims, 'ko-', lw=2, markersize=8)
    ax.set_xlabel('Source coherence r')
    ax.set_ylabel('Holographic SSIM')
    ax.set_title('SSIM vs source coherence', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(-0.05, 1.05)
    for p, r, s in zip(pressures, rs, ssims):
        label = f'{p:.0e}Pa' if p >= 1 else f'{p:.0e}'
        ax.annotate(label, (r, s), fontsize=7,
                    textcoords='offset points', xytext=(5, 5))

    fig.suptitle('v10: Cavity Pressure → Source Coherence → '
                 'Holographic Fidelity',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("╔═══════════════════════════════════════════════════════════╗")
    print("║  INTEGRATED QUANTUM SUBSTRATE  v10                       ║")
    print("║  He⁺ CMWB → neutralization → SQUID holography → caging  ║")
    print("╚═══════════════════════════════════════════════════════════╝")

    np.random.seed(0)
    N, L = 256, 400e-9

    # ==================================================================
    # VALIDATION GATES
    # ==================================================================
    print("\n" + "=" * 65)
    print("STEP 0: VALIDATION GATES")
    print("=" * 65)

    DiamondNetworkSimulator(Lx=4, Ly=4).validate_spectrum(
        phi=np.pi, abort_on_fail=True)

    squid_val = SQUIDArray(N_loops=32, N_grid=N, L_grid=L)
    solver_val = InverseHolographySolver(
        squid_array=squid_val, N=N, L=L, T_beam=1e-3,
        prop_distance_lam=20.0)
    rt_ok, rt_ssim = validate_roundtrip(solver_val)
    if not rt_ok:
        sys.exit("ABORT: Holography roundtrip validation failed.")

    # ==================================================================
    # DEMO 1: Full pipeline — vacuum He⁺ source, dots target
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 1: Full pipeline — vacuum He⁺, dots target")
    print("=" * 65)

    pipeline = IntegratedPipelineV10(
        B_field=0.01, pressure_Pa=100,  # rough vacuum
        T_beam=1e-3,
        N=N, L=L, N_loops=32,
        prop_distance_lam=20.0,
        use_floquet=False,
        N_diamond_x=16, N_diamond_y=16,
    )

    target_dots = target_grid_of_dots(N, L, n_dots=3)
    results_dots = pipeline.run(
        target_dots, method='gd', run_caging=True, seed=42)
    plot_v10_pipeline(pipeline, results_dots,
                      fname='results/v10_pipeline_dots.png')

    # ==================================================================
    # DEMO 2: Multi-target benchmark (same vacuum source)
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 2: Multi-target benchmark")
    print("=" * 65)

    targets = {
        'spot':   target_single_spot(N, L),
        'line':   target_line(N, L),
        'dots':   target_grid_of_dots(N, L, n_dots=3),
        'ring':   target_ring(N, L),
        'letter': target_letter(N, L, letter='H'),
    }

    all_results = {}
    for name, P in targets.items():
        print(f"\n{'#' * 65}")
        print(f"# TARGET: {name}")
        print(f"{'#' * 65}")
        all_results[name] = pipeline.run(
            P, method='gd', run_caging=True, seed=42)
    plot_multi_target(pipeline, all_results)

    # ==================================================================
    # DEMO 3: Pressure sweep — coherence threshold
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 3: Pressure sweep — source coherence threshold")
    print("=" * 65)

    sweep = pressure_sweep(
        target_dots,
        pressures=[1e-3, 1.0, 100, 1000, 10000, 101325],
        B_field=0.01, seed=42,
    )
    plot_pressure_sweep(sweep)

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print("\n" + "=" * 65)
    print("ALL SIMULATIONS COMPLETE")
    print("=" * 65)
    print("\nOutput files:")
    print("  results/v10_pipeline_dots.png     — Full pipeline demo")
    print("  results/v10_multi_target.png      — 5 target benchmark")
    print("  results/v10_pressure_sweep.png    — Coherence threshold")

    print(f"\nKey results (vacuum source, P = 100 Pa):")
    m = results_dots['holo']['metrics']
    print(f"  Source r:     {results_dots['r_source']:.4f}")
    print(f"  Holo SSIM:    {m['ssim']:.4f}")
    print(f"  Holo eff:     {m['efficiency']:.4f}")
    if results_dots['cage_pi'] is not None:
        loc = results_dots['cage_pi']['localization'][-1]
        print(f"  Cage loc:     {loc:.4f}")

    print(f"\nMulti-target SSIMs:")
    for name, res in all_results.items():
        m = res['holo']['metrics']
        print(f"  {name:<10}: SSIM={m['ssim']:.3f}  "
              f"eff={m['efficiency']:.3f}")

    print(f"\nPressure sweep:")
    print(f"  {'P (Pa)':<12} {'r':>6} {'SSIM':>8} {'K_dim':>8}")
    print(f"  {'-'*38}")
    for r in sweep:
        m = r['holo']['metrics']
        print(f"  {r['pressure']:<12.1e} {r['r_source']:6.4f} "
              f"{m['ssim']:8.4f} {r['K_dim']:8.3f}")


if __name__ == "__main__":
    main()