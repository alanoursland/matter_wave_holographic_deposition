"""
v11 Design Studies (fable5 Phase 4)
====================================

Three studies on the corrected pipeline:

  T14 — verification that transport-aware target conditioning stops the
        optimizer from chasing undeliverable content (efficiency check).
  T13 — species/velocity trade study: λ_dB, deliverable bandwidth, SSIM
        at the fixed 32×32 array, AB coupling, image-charge sensitivity,
        and space-charge-limited current, over (species, E_k).
  T15 — caging reframed as per-adatom migration suppression: single
        A-site initial states, escape probability vs flux Φ, disorder W,
        and flux error δΦ (the hardware tolerance spec).

Run:  python v11_design_studies.py
"""

import os
import time
import numpy as np
import torch
import matplotlib.pyplot as plt

from iqs.sources.coherent_matterwave import (
    CoherentMatterwaveBeam, CavityGeometry, SPECIES,
    hbar, k_B, e_C, eps_0,
)
from iqs.holography import (
    SQUIDArray, InverseHolographySolver,
    target_grid_of_dots, target_line, target_letter,
    smooth_target, compute_metrics,
)
from iqs.lattices.diamond import DiamondNetwork
from iqs.numerics.device import get_device

_device = get_device()
os.makedirs('results', exist_ok=True)

# Fixed hardware geometry (the v10 stage, held constant across studies)
N_GRID   = 256
L_AP     = 400e-9      # aperture / substrate window [m]
Z_PROP   = 978.2e-9    # array → substrate distance [m]
N_LOOPS  = 32


def _solve_dots_ssim(lam, n_iter=300, seed=0, transport_aware=True):
    """Clean-solver dots SSIM at wavelength lam on the fixed geometry."""
    np.random.seed(seed)
    torch.manual_seed(seed)
    k0 = 2 * np.pi / lam
    T_equiv = (hbar * k0)**2 / (2 * (SPECIES['He+']['mass']) * k_B)
    # solver derives k0 from T_beam for He+; construct via T_beam that
    # reproduces this k0 exactly for the He+ mass it uses internally.
    squid = SQUIDArray(N_loops=N_LOOPS, N_grid=N_GRID, L_grid=L_AP)
    solver = InverseHolographySolver(
        squid_array=squid, N=N_GRID, L=L_AP, T_beam=T_equiv,
        prop_distance_lam=Z_PROP / lam)
    target = target_grid_of_dots(N_GRID, L_AP, n_dots=3)
    kwargs = dict(N_loops=N_LOOPS, corner_radius=0.03, sigma=2)
    if transport_aware:
        kwargs.update(k0=solver.k0, L=L_AP, z=solver.z)
    tgt = smooth_target(target, **kwargs)
    r = solver.solve_gradient_descent(tgt, n_iter=n_iter, lr=0.05,
                                      verbose=False)
    return compute_metrics(r['achieved'], tgt)


# ===========================================================================
# T14 — conditioning verification
# ===========================================================================

def t14_conditioning_comparison(n_iter=500, seed=0, verbose=True):
    """Old (array-Nyquist-only) vs new (transport-aware) conditioning."""
    if verbose:
        print("\n" + "=" * 65)
        print("T14: TARGET CONDITIONING — array-only vs transport-aware")
        print("=" * 65)

    squid = SQUIDArray(N_loops=N_LOOPS, N_grid=N_GRID, L_grid=L_AP)
    solver = InverseHolographySolver(
        squid_array=squid, N=N_GRID, L=L_AP, T_beam=1e-3,
        prop_distance_lam=20.0)

    targets = {
        'dots':   target_grid_of_dots(N_GRID, L_AP, n_dots=3),
        'line':   target_line(N_GRID, L_AP),
        'letter': target_letter(N_GRID, L_AP, letter='H'),
    }

    rows = []
    for name, P in targets.items():
        for mode in ('array-only', 'transport-aware'):
            np.random.seed(seed)
            torch.manual_seed(seed)
            kwargs = dict(N_loops=N_LOOPS, corner_radius=0.03, sigma=2)
            if mode == 'transport-aware':
                kwargs.update(k0=solver.k0, L=L_AP, z=solver.z)
            tgt = smooth_target(P, **kwargs)
            r = solver.solve_gradient_descent(tgt, n_iter=n_iter, lr=0.05,
                                              verbose=False)
            m = compute_metrics(r['achieved'], tgt)
            rows.append((name, mode, m['ssim'], m['efficiency']))
            if verbose:
                print(f"  {name:<8} {mode:<16} SSIM={m['ssim']:.4f}  "
                      f"eff={m['efficiency']:.4f}")
    return rows


# ===========================================================================
# T13 — species/velocity trade study
# ===========================================================================

TRADE_SPECIES = ['He+', 'Li+', 'Na+', 'Ca+', 'Rb+']
TRADE_EK_EV = [8.617e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0]
IMAGE_D = 100e-9  # reference ion–plane distance for image-charge energy


def t13_trade_study(E_k_list=None, species_list=None, n_iter=300,
                    verbose=True):
    """Sweep (species, E_k); return list of row dicts + recommendation."""
    if E_k_list is None:
        E_k_list = TRADE_EK_EV
    if species_list is None:
        species_list = TRADE_SPECIES

    if verbose:
        print("\n" + "=" * 65)
        print("T13: SPECIES / VELOCITY TRADE STUDY")
        print(f"     fixed stage: L = {L_AP*1e9:.0f} nm, "
              f"z = {Z_PROP*1e9:.0f} nm, {N_LOOPS}×{N_LOOPS} array")
        print("=" * 65)

    sin_geo = np.sin(np.arctan(L_AP / Z_PROP))
    f_nyq = N_LOOPS / 2.0
    E_image = e_C**2 / (16 * np.pi * eps_0 * IMAGE_D)   # J

    ssim_cache = {}
    rows = []
    for sp in species_list:
        for E_eV in E_k_list:
            sim = CoherentMatterwaveBeam(
                species=sp, E_kinetic_eV=E_eV, B_field=0.01,
                cavity=CavityGeometry(pressure_Pa=1e-3), dE_frac=0.01)
            lam = sim.lam_dB
            # Deliverable bandwidth (field-level cycles/aperture)
            f_transport = min(2 * np.pi / lam * sin_geo * L_AP / (2 * np.pi),
                              f_nyq)
            # SSIM depends only on λ at fixed geometry — cache
            lam_key = round(np.log10(lam), 3)
            if lam_key not in ssim_cache:
                # skip absurdly long wavelengths (> aperture: nothing
                # propagates usefully) to save time
                if lam > L_AP:
                    ssim_cache[lam_key] = {'ssim': 0.0, 'efficiency': 0.0}
                else:
                    ssim_cache[lam_key] = _solve_dots_ssim(lam, n_iter=n_iter)
            m = ssim_cache[lam_key]
            row = {
                'species': sp,
                'E_k_eV': E_eV,
                'v': sim.v,
                'lam_nm': lam * 1e9,
                'f_transport_ca': f_transport,
                'ssim_dots': m['ssim'],
                'eff_dots': m['efficiency'],
                'K_dim_1mPa': sim.K_dim,
                'image_ratio': E_image / sim.E_k,
                'I_max_single_A': abs(sim.charge) * sim.v / Z_PROP,
                'mfp_1mPa_m': sim.mfp_beam,
            }
            rows.append(row)
            if verbose:
                print(f"  {sp:<5} E={E_eV:8.2e} eV  v={sim.v:9.3g} m/s  "
                      f"λ={lam*1e9:9.3g} nm  f={f_transport:5.2f} c/a  "
                      f"SSIM={m['ssim']:.3f}  "
                      f"img={E_image/sim.E_k:8.2e}  "
                      f"I₁={row['I_max_single_A']:.2e} A")

    # Recommendation: per species, the smallest E_k whose transport
    # bandwidth reaches 0.8× the array Nyquist (array becomes the limit —
    # more energy buys nothing at this stage geometry).
    corners = {}
    for sp in species_list:
        sp_rows = [r for r in rows if r['species'] == sp]
        ok = [r for r in sp_rows if r['f_transport_ca'] >= 0.8 * f_nyq]
        if ok:
            corners[sp] = min(ok, key=lambda r: r['E_k_eV'])
    # Among corners prefer smallest image-charge ratio, then largest
    # single-ion current (both favour the fastest corner beam).
    recommended = None
    if corners:
        recommended = min(corners.values(),
                          key=lambda r: (r['image_ratio'],
                                         -r['I_max_single_A']))
    return rows, corners, recommended


def plot_trade_study(rows, corners, recommended,
                     fname='results/v11_trade_study.png'):
    print(f"\n  Plotting → {fname}")
    species = sorted({r['species'] for r in rows},
                     key=lambda s: SPECIES[s]['mass'])
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    for sp in species:
        sr = sorted([r for r in rows if r['species'] == sp],
                    key=lambda r: r['E_k_eV'])
        E = [r['E_k_eV'] for r in sr]
        axes[0, 0].loglog(E, [r['lam_nm'] for r in sr], 'o-', label=sp)
        axes[0, 1].semilogx(E, [r['f_transport_ca'] for r in sr], 'o-',
                            label=sp)
        axes[1, 0].semilogx(E, [r['ssim_dots'] for r in sr], 'o-', label=sp)
        axes[1, 1].loglog(E, [r['image_ratio'] for r in sr], 'o-', label=sp)

    axes[0, 0].set_ylabel('λ_dB (nm)')
    axes[0, 0].axhline(L_AP * 1e9, color='k', ls=':', alpha=0.5)
    axes[0, 0].set_title('de Broglie wavelength')
    axes[0, 1].set_ylabel('deliverable bandwidth (cycles/aperture)')
    axes[0, 1].axhline(0.8 * N_LOOPS / 2, color='k', ls='--', alpha=0.5,
                       label='0.8 × array Nyquist')
    axes[0, 1].set_title('Transport bandwidth (field)')
    axes[1, 0].set_ylabel('SSIM (dots, clean solver)')
    axes[1, 0].set_title('Fidelity at fixed 32×32 array')
    axes[1, 1].set_ylabel('E_image(100 nm) / E_k')
    axes[1, 1].axhline(1.0, color='r', ls='--', alpha=0.5)
    axes[1, 1].set_title('Image-charge sensitivity')
    for ax in axes.ravel():
        ax.set_xlabel('E_k (eV)')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    if recommended:
        axes[1, 0].scatter([recommended['E_k_eV']],
                           [recommended['ssim_dots']],
                           s=250, facecolors='none', edgecolors='g', lw=2,
                           zorder=5)
    fig.suptitle('v11 species/velocity trade study (T13) — '
                 f"L={L_AP*1e9:.0f} nm, z={Z_PROP*1e9:.0f} nm",
                 fontweight='bold')
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close(fig)


# ===========================================================================
# T15 — per-adatom migration suppression
# ===========================================================================

def single_adatom_escape(lattice, phi, disorder_W=0.0, seed=0, T=40.0):
    """Escape probability of one adatom initialized on the central A-site.

    P_escape = 1 − P(within the 3×3 cell neighbourhood) at time T.
    """
    psi0 = np.zeros(lattice.dim, dtype=complex)
    psi0[lattice._site(lattice.Lx // 2, lattice.Ly // 2, 0)] = 1.0
    res = lattice.evolve(psi0, phi=phi, T=T,
                         disorder_W=disorder_W, disorder_seed=seed,
                         verbose=False)
    return 1.0 - float(res['localization'][-1])


def t15_adatom_caging(Lx=16, Ly=16, n_seeds=10, verbose=True):
    """P_escape(Φ), P_escape(W at Φ=π), P_escape(δΦ) curves."""
    if verbose:
        print("\n" + "=" * 65)
        print("T15: PER-ADATOM MIGRATION SUPPRESSION")
        print(f"     single A-site initial state, {Lx}×{Ly} diamond lattice")
        print("=" * 65)
    lat = DiamondNetwork(Lx=Lx, Ly=Ly)

    # (a) Escape vs flux
    phis = np.linspace(0, np.pi, 13)
    esc_phi = [single_adatom_escape(lat, p) for p in phis]
    if verbose:
        print("\n  P_escape(Φ):")
        for p, e in zip(phis, esc_phi):
            print(f"    Φ = {p/np.pi:4.2f}π : P_escape = {e:.4f}")

    # (b) Escape vs disorder at Φ = π (ensemble over realizations)
    Ws = [0.0, 0.2, 0.5, 1.0, 2.0]
    esc_W_mean, esc_W_std = [], []
    for W in Ws:
        vals = [single_adatom_escape(lat, np.pi, disorder_W=W, seed=s)
                for s in range(n_seeds if W > 0 else 1)]
        esc_W_mean.append(np.mean(vals))
        esc_W_std.append(np.std(vals))
    if verbose:
        print("\n  P_escape(W) at Φ = π:")
        for W, m, s in zip(Ws, esc_W_mean, esc_W_std):
            print(f"    W = {W:4.1f} J : P_escape = {m:.4f} ± {s:.4f}")

    # (c) Escape vs flux error δΦ = |Φ − π| (hardware tolerance spec)
    dphis = np.array([0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20]) * np.pi
    esc_dphi = [single_adatom_escape(lat, np.pi - d) for d in dphis]
    # tolerance: largest δΦ with P_escape < 5%
    tol = None
    for d, e in zip(dphis, esc_dphi):
        if e < 0.05:
            tol = d
    if verbose:
        print("\n  P_escape(δΦ) at W = 0:")
        for d, e in zip(dphis, esc_dphi):
            print(f"    δΦ = {d/np.pi:5.3f}π : P_escape = {e:.4f}")
        if tol is not None:
            print(f"\n  Flux tolerance (P_escape < 5%): δΦ ≤ {tol/np.pi:.3f}π "
                  f"= {tol:.4f} rad")

    return {
        'phis': phis, 'esc_phi': np.array(esc_phi),
        'Ws': Ws, 'esc_W_mean': np.array(esc_W_mean),
        'esc_W_std': np.array(esc_W_std),
        'dphis': dphis, 'esc_dphi': np.array(esc_dphi),
        'flux_tolerance_rad': tol,
    }


def plot_adatom_caging(res, fname='results/v11_caging_adatom.png'):
    print(f"\n  Plotting → {fname}")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    ax = axes[0]
    ax.plot(res['phis'] / np.pi, res['esc_phi'], 'bo-', lw=2)
    ax.set_xlabel('Φ / π')
    ax.set_ylabel('P_escape (single adatom)')
    ax.set_title('Escape vs flux', fontweight='bold')

    ax = axes[1]
    ax.errorbar(res['Ws'], res['esc_W_mean'], yerr=res['esc_W_std'],
                fmt='rs-', lw=2, capsize=3)
    ax.set_xlabel('Disorder W (J)')
    ax.set_ylabel('P_escape at Φ = π')
    ax.set_title('Escape vs on-site disorder', fontweight='bold')

    ax = axes[2]
    ax.semilogx(np.maximum(res['dphis'] / np.pi, 1e-4), res['esc_dphi'],
                'g^-', lw=2)
    ax.axhline(0.05, color='k', ls='--', alpha=0.5, label='5% spec')
    ax.set_xlabel('flux error δΦ / π')
    ax.set_ylabel('P_escape')
    ax.set_title('Escape vs flux error (tolerance spec)',
                 fontweight='bold')
    ax.legend()
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)
    fig.suptitle('Per-adatom migration suppression (T15) — '
                 'single A-site initial states', fontweight='bold')
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close(fig)


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("╔═══════════════════════════════════════════════════════════╗")
    print("║  v11 DESIGN STUDIES (fable5 Phase 4: T13/T14/T15)        ║")
    print("╚═══════════════════════════════════════════════════════════╝")
    t0 = time.time()

    rows_t14 = t14_conditioning_comparison()

    rows, corners, rec = t13_trade_study()
    plot_trade_study(rows, corners, rec)

    print("\n  Operating corners (smallest E_k where transport bandwidth "
          "reaches 0.8× array Nyquist):")
    for sp, r in corners.items():
        print(f"    {sp:<5} E = {r['E_k_eV']:.2e} eV, v = {r['v']:.3g} m/s, "
              f"λ = {r['lam_nm']:.3g} nm, SSIM = {r['ssim_dots']:.3f}, "
              f"img = {r['image_ratio']:.2e}, "
              f"I₁ = {r['I_max_single_A']:.2e} A")
    if rec:
        print(f"\n  RECOMMENDED GENERATION-1 CORNER: {rec['species']} at "
              f"{rec['E_k_eV']:.2e} eV "
              f"(v = {rec['v']:.3g} m/s, λ = {rec['lam_nm']:.3g} nm)")
        print(f"    SSIM(dots) = {rec['ssim_dots']:.3f}, "
              f"image ratio = {rec['image_ratio']:.2e}, "
              f"single-ion current ≤ {rec['I_max_single_A']:.2e} A")

    res15 = t15_adatom_caging()
    plot_adatom_caging(res15)

    print(f"\nAll studies complete in {time.time()-t0:.0f}s.")
    print("Output: results/v11_trade_study.png, results/v11_caging_adatom.png")


if __name__ == '__main__':
    main()
