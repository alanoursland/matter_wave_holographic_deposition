"""
T22 — Stage Re-Derivation Sweep
================================

Joint sweep of the stage geometry (propagation distance z, wavelength λ
via He⁺ kinetic energy, array size N_loops) at fixed aperture
L = 400 nm, on the corrected pipeline with transport-aware target
conditioning (T14).

Motivation: every ceiling found since Phase 1 — the 3.1 c/a transport
band, the ~10⁻³ arrived power, the 0.12 SSIM coherence gap (dim
patterns are noise-fragile), T21's ~2.3% actuator floor — traces to the
v10 stage geometry (z = 20λ = 978 nm at λ = 48.9 nm). This study maps
the (fidelity × delivered power) landscape and recommends the stage at
which T18/T19 should be evaluated.

Metrics per grid point (dots target):
  ssim_ref   — achieved vs the FIXED array-band reference target
               (12.8 c/a conditioning): how close the print is to the
               pattern we actually want.  The stage-quality number.
  ssim_self  — achieved vs the point's own transport-conditioned
               target: solver quality at that geometry.
  arrived    — fraction of input beam power landing in the L×L frame
               (the E8 accounting; input is unit-normalized).
  delivered  — arrived × in-mask efficiency: fraction of the beam that
               lands *on the pattern*.

Then the Pareto leaders get an ensemble-noise evaluation
(σ_θ = 0.079 rad, the Q-limited coherence ceiling) to confirm that
brighter stages are also noise-robust (the Phase 2 mechanism).

Run:  python t22_stage_sweep.py
"""

import os
import time
import numpy as np
import torch
import matplotlib.pyplot as plt

from inverse_holography import (
    SQUIDArray, InverseHolographySolver,
    target_grid_of_dots, smooth_target, compute_metrics,
)
from coherent_matterwave_beam import sample_phase_noise
from iqs.constants import hbar, k_B, m_He
from iqs.numerics.metrics import ssim_score
from iqs.numerics.device import get_device

_device = get_device()
os.makedirs('results', exist_ok=True)

N_GRID = 256
L_AP   = 400e-9

# Sweep axes.  λ from He⁺ kinetic energy; z in absolute nm (hardware
# gap), not multiples of λ — the stage is a physical object.
LAMBDAS_NM = [48.9, 14.4, 4.54]          # 1 mK, 1e-6 eV, 1e-5 eV He⁺
Z_NM       = [100, 245, 489, 978]        # array→substrate gap
N_LOOPS_LIST = [16, 32, 64]

SIGMA_THETA = 0.079
XI_PERP     = 50e-9
M_ENSEMBLE  = 50


def T_beam_for_lambda(lam):
    """He⁺ beam temperature whose thermal velocity gives λ_dB = lam."""
    k0 = 2 * np.pi / lam
    return (hbar * k0)**2 / (2 * m_He * k_B)


def solve_point(lam, z, n_loops, n_iter=400, seed=0):
    np.random.seed(seed)
    torch.manual_seed(seed)
    squid = SQUIDArray(N_loops=n_loops, N_grid=N_GRID, L_grid=L_AP)
    solver = InverseHolographySolver(
        squid_array=squid, N=N_GRID, L=L_AP,
        T_beam=T_beam_for_lambda(lam), prop_distance_lam=z / lam)

    raw = target_grid_of_dots(N_GRID, L_AP, n_dots=3)
    # fixed reference: array-band conditioning at the *32-loop* Nyquist,
    # independent of the sweep point — one common yardstick
    ref = smooth_target(raw, N_loops=32, corner_radius=0.03, sigma=2)
    # per-point deliverable-band conditioning (T14)
    cond = smooth_target(raw, N_loops=n_loops, corner_radius=0.03,
                         sigma=2, k0=solver.k0, L=L_AP, z=solver.z)

    r = solver.solve_gradient_descent(cond, n_iter=n_iter, lr=0.05,
                                      verbose=False)
    achieved = r['achieved']
    m_self = compute_metrics(achieved, cond)
    arrived = float(achieved.sum() * solver.dx**2)   # input norm = 1
    row = {
        'lam_nm': lam * 1e9 if lam < 1e-3 else lam,  # tolerate nm input
        'z_nm': z * 1e9 if z < 1e-3 else z,
        'n_loops': n_loops,
        'f_transport_ca': float(
            solver.k0 * np.sin(np.arctan(L_AP / solver.z))
            * L_AP / (2 * np.pi)),
        'ssim_self': m_self['ssim'],
        'ssim_ref': float(ssim_score(ref, achieved)),
        'eff': m_self['efficiency'],
        'arrived': arrived,
        'delivered': arrived * m_self['efficiency'],
    }
    return row, solver, r['phase_screen'], cond


def ensemble_actual(solver, phase_screen, cond_target, seed=0):
    """SSIM of the ensemble-noise deposition vs the conditioned target."""
    rng = np.random.default_rng(seed)
    screen_t = torch.tensor(phase_screen, dtype=torch.float64,
                            device=_device)
    T = SQUIDArray.phase_screen_to_transmission(screen_t)
    prof = solver.psi_in_np
    dens = np.zeros((N_GRID, N_GRID))
    for _ in range(M_ENSEMBLE):
        noise = sample_phase_noise(N_GRID, SIGMA_THETA,
                                   XI_PERP / solver.dx, rng)
        psi = prof * np.exp(1j * noise)
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * solver.dx**2)
        psi_t = torch.tensor(psi, dtype=torch.complex128, device=_device)
        out = solver._propagate_torch(psi_t * T)
        dens += torch.abs(out).cpu().numpy()**2
    dens /= M_ENSEMBLE
    return float(ssim_score(cond_target, dens))


def run_sweep(verbose=True):
    rows = []
    keep = {}
    for lam_nm in LAMBDAS_NM:
        for z_nm in Z_NM:
            for n_loops in N_LOOPS_LIST:
                row, solver, screen, cond = solve_point(
                    lam_nm * 1e-9, z_nm * 1e-9, n_loops)
                row['lam_nm'], row['z_nm'] = lam_nm, z_nm
                rows.append(row)
                keep[(lam_nm, z_nm, n_loops)] = (solver, screen, cond)
                if verbose:
                    print(f"  λ={lam_nm:5.1f} nm  z={z_nm:4.0f} nm  "
                          f"loops={n_loops:2d}  "
                          f"f={row['f_transport_ca']:5.2f} c/a  "
                          f"SSIM_ref={row['ssim_ref']:.3f}  "
                          f"SSIM_self={row['ssim_self']:.3f}  "
                          f"arrived={row['arrived']:.3e}  "
                          f"delivered={row['delivered']:.3e}",
                          flush=True)
    return rows, keep


def pareto_front(rows, xkey='delivered', ykey='ssim_ref'):
    """Points not dominated in (higher x, higher y)."""
    front = []
    for r in rows:
        if not any((o[xkey] >= r[xkey] and o[ykey] >= r[ykey]
                    and (o[xkey] > r[xkey] or o[ykey] > r[ykey]))
                   for o in rows):
            front.append(r)
    return sorted(front, key=lambda r: r[xkey])


def plot_sweep(rows, front, v10_row, fname='results/t22_stage_sweep.png'):
    print(f"\n  Plotting → {fname}")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    markers = {16: 'o', 32: 's', 64: '^'}
    colors = {48.9: '#c62828', 14.4: '#1565c0', 4.54: '#2e7d32'}

    ax = axes[0]
    for r in rows:
        ax.scatter(r['delivered'], r['ssim_ref'],
                   marker=markers[r['n_loops']],
                   color=colors[r['lam_nm']], s=60, alpha=0.8)
    fx = [r['delivered'] for r in front]
    fy = [r['ssim_ref'] for r in front]
    ax.plot(fx, fy, 'k--', alpha=0.5, lw=1.5, label='Pareto front')
    ax.scatter([v10_row['delivered']], [v10_row['ssim_ref']],
               s=260, facecolors='none', edgecolors='r', lw=2,
               label='v10 stage')
    ax.set_xscale('log')
    ax.set_xlabel('delivered-to-mask power fraction')
    ax.set_ylabel('SSIM vs fixed reference target')
    ax.set_title('Fidelity × delivered power', fontweight='bold')
    for lam, c in colors.items():
        ax.scatter([], [], color=c, label=f'λ = {lam} nm')
    for nl, mk in markers.items():
        ax.scatter([], [], color='gray', marker=mk, label=f'{nl}² loops')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3, which='both')

    ax = axes[1]
    for lam_nm in LAMBDAS_NM:
        for n_loops in N_LOOPS_LIST:
            pts = sorted([r for r in rows if r['lam_nm'] == lam_nm
                          and r['n_loops'] == n_loops],
                         key=lambda r: r['z_nm'])
            ax.plot([p['z_nm'] for p in pts],
                    [p['ssim_ref'] for p in pts],
                    marker=markers[n_loops], color=colors[lam_nm],
                    lw=1.5, alpha=0.8)
    ax.set_xlabel('array–substrate gap z (nm)')
    ax.set_ylabel('SSIM vs fixed reference target')
    ax.set_title('Fidelity vs stage gap', fontweight='bold')
    ax.grid(True, alpha=0.3)

    fig.suptitle('T22 — stage re-derivation sweep '
                 f'(dots target, L = {L_AP*1e9:.0f} nm)',
                 fontweight='bold')
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close(fig)


def main():
    print("=" * 70)
    print("T22: STAGE RE-DERIVATION SWEEP")
    print(f"     L = {L_AP*1e9:.0f} nm; λ ∈ {LAMBDAS_NM} nm; "
          f"z ∈ {Z_NM} nm; loops ∈ {N_LOOPS_LIST}")
    print("=" * 70)
    t0 = time.time()

    rows, keep = run_sweep()
    front = pareto_front(rows)
    v10_row = next(r for r in rows if r['lam_nm'] == 48.9
                   and r['z_nm'] == 978 and r['n_loops'] == 32)
    plot_sweep(rows, front, v10_row)

    print("\n  Pareto front (delivered ↑, SSIM_ref ↑):")
    for r in front:
        print(f"    λ={r['lam_nm']:5.1f}  z={r['z_nm']:4.0f}  "
              f"loops={r['n_loops']:2d}  SSIM_ref={r['ssim_ref']:.3f}  "
              f"delivered={r['delivered']:.3e}")

    # Recommendation rule: among stages delivering ≥ 1% of the beam to
    # the mask, take the highest fixed-reference SSIM.
    viable = [r for r in rows if r['delivered'] >= 0.01]
    rec = max(viable, key=lambda r: r['ssim_ref']) if viable else None

    # Ensemble-noise check on the top-3 Pareto leaders by SSIM_ref + v10
    print("\n  Ensemble-noise evaluation (σ_θ = 0.079 rad) of leaders:")
    leaders = sorted(front, key=lambda r: -r['ssim_ref'])[:3]
    for r in leaders + [v10_row]:
        key = (r['lam_nm'], r['z_nm'], r['n_loops'])
        solver, screen, cond = keep[key]
        s_act = ensemble_actual(solver, screen, cond)
        r['ssim_actual'] = s_act
        gap = r['ssim_self'] - s_act
        tag = ' (v10 stage)' if r is v10_row else ''
        print(f"    λ={r['lam_nm']:5.1f}  z={r['z_nm']:4.0f}  "
              f"loops={r['n_loops']:2d}: SSIM clean={r['ssim_self']:.3f} "
              f"actual={s_act:.3f}  gap={gap:.3f}{tag}")

    if rec is not None:
        print(f"\n  RECOMMENDED STAGE: λ = {rec['lam_nm']} nm, "
              f"z = {rec['z_nm']:.0f} nm, {rec['n_loops']}² loops")
        print(f"    SSIM_ref = {rec['ssim_ref']:.3f}, "
              f"delivered = {rec['delivered']:.3e}, "
              f"transport band = {rec['f_transport_ca']:.1f} c/a")
        print(f"    vs v10 stage: SSIM_ref = {v10_row['ssim_ref']:.3f}, "
              f"delivered = {v10_row['delivered']:.3e}")
    print(f"\nDone in {time.time()-t0:.0f}s. "
          f"Output: results/t22_stage_sweep.png")


if __name__ == '__main__':
    main()
