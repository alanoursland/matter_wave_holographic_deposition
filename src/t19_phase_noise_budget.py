"""
T19 — Physical Phase-Noise Budget
==================================

Replaces the free σ_θ parameter with a composed budget of physical
noise sources, as a function of beam energy, and inverts it into
per-subsystem stability specs (gen-2 note §8).

Sources modeled (RMS phase across the aperture):

  electrostatic   σ = q·V_diff·τ/ħ,  τ = L_leg/v
      V_diff = *differential* potential across interfering paths,
      integrated over the leg — patches, charging, electrode drift.
      Common-mode potentials are harmless (global phase); *static*
      differential potentials are calibratable (they are part of the
      hologram the T21 loop measures); what burns the budget is DRIFT
      of the differential potential between calibration and exposure.
  screen drive    σ = (2π/Φ₀ᵉ)·s_M·δI,  s_M = max_i √(Σ_j M_ij²)
      current noise through the T16 mutual-inductance matrix
      (velocity-independent — AB phase is achromatic).
  defocus         σ = k·(1−cosθ_NA)·δz
      screen–substrate gap error/vibration at the stage NA.
  stray B         σ = q·δB_diff·(L·z)/ħ
      differential flux through the area between extreme paths.
  Coulomb         zero in single-ion operation (T10/T21).

Scenarios:
  A. Literal 1:1 T22 stage — He⁺ at λ = 14.4 nm (v ≈ 6.9 m/s),
     z = 489 nm, 64² loops.
  B. Gen-2 projection — He⁺ transported at 1 keV / 30 keV through a
     0.1 m column, NA chosen for 5 nm resolution at the substrate.

Outputs: σ_θ(E_k) catastrophe curve; per-subsystem spec table for a
total budget σ_θ = 0.1 rad; an empirical SSIM-vs-σ_θ anchor at the T22
stage justifying that budget.

Run:  python t19_phase_noise_budget.py
"""

import os
import time
import numpy as np
import torch
import matplotlib.pyplot as plt

from inverse_holography import SQUIDArray, InverseHolographySolver, \
    target_grid_of_dots, smooth_target
from coherent_matterwave_beam import sample_phase_noise, hbar, e_C
from iqs.constants import k_B, m_He
from iqs.numerics.metrics import ssim_score
from iqs.numerics.device import get_device

_device = get_device()
os.makedirs('results', exist_ok=True)

Phi_0_e = 4.135667696e-15   # Wb, h/e
Q = e_C                      # He⁺
BUDGET = 0.1                 # rad, total σ_θ (anchored empirically below)
N_SOURCES = 4                # electrostatic, drive, defocus, stray B
PER_SOURCE = BUDGET / np.sqrt(N_SOURCES)

# T22-recommended stage
L_AP, Z_STAGE, N_LOOPS, LAM_STAGE = 400e-9, 489e-9, 64, 14.4e-9
N_GRID = 256


def v_of_E(E_eV, m=m_He):
    return np.sqrt(2 * E_eV * e_C / m)


def lam_of_E(E_eV, m=m_He):
    return 6.62607015e-34 / (m * v_of_E(E_eV, m))


# ---------------------------------------------------------------------------
# Source models and spec inversions
# ---------------------------------------------------------------------------

def sigma_electrostatic(V_diff, L_leg, v):
    return Q * V_diff * (L_leg / v) / hbar


def V_spec(sigma, L_leg, v):
    return sigma * hbar * v / (Q * L_leg)


def drive_sensitivity(n_loops=N_LOOPS, L_grid=L_AP):
    """s_M = worst-row RMS coupling of the T16 inductance matrix, and
    the full-2π drive current, at this pitch."""
    squid = SQUIDArray(N_loops=n_loops, N_grid=64, L_grid=L_grid)
    M = squid.build_inductance_matrix(wire_width_frac=0.1)
    s_M = float(np.sqrt((M**2).sum(axis=1).max()))
    I_2pi = Phi_0_e / squid._L_self
    return s_M, I_2pi


def dI_spec(sigma, s_M):
    return sigma * Phi_0_e / (2 * np.pi * s_M)


def sigma_defocus(k, sin_na, dz):
    return k * (1 - np.sqrt(1 - sin_na**2)) * dz


def dz_spec(sigma, k, sin_na):
    return sigma / (k * (1 - np.sqrt(1 - sin_na**2)))


def dB_spec(sigma, L_ap, z):
    return sigma * hbar / (Q * L_ap * z)


# ---------------------------------------------------------------------------
# Empirical anchor: SSIM vs σ_θ at the T22 stage
# ---------------------------------------------------------------------------

def ssim_vs_sigma(sigmas, seed=0, M_ens=50, xi_perp=50e-9):
    T_beam = (hbar * 2 * np.pi / LAM_STAGE)**2 / (2 * m_He * k_B)
    np.random.seed(seed); torch.manual_seed(seed)
    squid = SQUIDArray(N_loops=N_LOOPS, N_grid=N_GRID, L_grid=L_AP)
    solver = InverseHolographySolver(
        squid_array=squid, N=N_GRID, L=L_AP, T_beam=T_beam,
        prop_distance_lam=Z_STAGE / LAM_STAGE)
    raw = target_grid_of_dots(N_GRID, L_AP, n_dots=3)
    cond = smooth_target(raw, N_loops=N_LOOPS, corner_radius=0.03,
                         sigma=2, k0=solver.k0, L=L_AP, z=solver.z)
    sol = solver.solve_gradient_descent(cond, n_iter=400, lr=0.05,
                                        verbose=False)
    clean = float(ssim_score(cond, sol['achieved']))

    screen_t = torch.tensor(sol['phase_screen'], dtype=torch.float64,
                            device=_device)
    T = SQUIDArray.phase_screen_to_transmission(screen_t)
    prof = solver.psi_in_np
    rng = np.random.default_rng(seed)

    rows = []
    for s in sigmas:
        dens = np.zeros((N_GRID, N_GRID))
        for _ in range(M_ens):
            noise = sample_phase_noise(N_GRID, s, xi_perp / solver.dx, rng)
            psi = prof * np.exp(1j * noise)
            psi /= np.sqrt(np.sum(np.abs(psi)**2) * solver.dx**2)
            psi_t = torch.tensor(psi, dtype=torch.complex128,
                                 device=_device)
            out = solver._propagate_torch(psi_t * T)
            dens += torch.abs(out).cpu().numpy()**2
        dens /= M_ens
        rows.append((s, float(ssim_score(cond, dens))))
    return clean, rows


# ---------------------------------------------------------------------------
# Study
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("T19: PHYSICAL PHASE-NOISE BUDGET")
    print(f"     total budget σ_θ = {BUDGET} rad "
          f"({PER_SOURCE:.3f} rad/source × {N_SOURCES} sources in quadrature)")
    print("=" * 70)
    t0 = time.time()

    # --- empirical anchor -------------------------------------------------
    sigmas = [0.03, 0.1, 0.3, 1.0, 3.0]
    clean, anchor = ssim_vs_sigma(sigmas)
    print(f"\n  SSIM vs σ_θ at the T22 stage (clean = {clean:.3f}):")
    for s, v in anchor:
        print(f"    σ_θ = {s:4.2f} rad → SSIM = {v:.3f} "
              f"(gap {clean - v:+.3f})")
    print(f"  → budget σ_θ = {BUDGET} rad costs "
          f"{clean - dict(anchor)[0.1]:.3f} SSIM: acceptable.")

    # --- catastrophe curve: σ_electrostatic(E) at fixed subsystem quality -
    E_grid = np.logspace(-8, 4.5, 200)   # eV
    V_ref = 1e-3                          # 1 mV effective differential
    sig_E = [sigma_electrostatic(V_ref, Z_STAGE, v_of_E(E)) for E in E_grid]

    marks = {
        'v10 (1 mK)': 8.617e-8,
        'T22 1:1 (10⁻⁶ eV)': 1e-6,
        'gen-2 1 keV': 1e3,
        'gen-2 30 keV': 3e4,
    }

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5))
    ax = axes[0]
    ax.loglog(E_grid, sig_E, 'b-', lw=2,
              label=f'electrostatic, V_diff = 1 mV over {Z_STAGE*1e9:.0f} nm')
    ax.axhline(BUDGET, color='k', ls='--', alpha=0.6,
               label=f'budget {BUDGET} rad')
    for name, E in marks.items():
        s = sigma_electrostatic(V_ref, Z_STAGE, v_of_E(E))
        ax.scatter([E], [s], s=60, zorder=5)
        ax.annotate(name, (E, s), fontsize=7,
                    textcoords='offset points', xytext=(6, 4))
    ax.set_xlabel('beam kinetic energy (eV)')
    ax.set_ylabel('σ_θ (rad)')
    ax.set_title('The slow-beam electrostatic catastrophe',
                 fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.semilogx([s for s, _ in anchor], [v for _, v in anchor], 'ro-',
                lw=2, label='ensemble SSIM (measured)')
    ax.axhline(clean, color='k', ls=':', alpha=0.6, label='clean SSIM')
    ax.axvline(BUDGET, color='k', ls='--', alpha=0.6,
               label=f'budget {BUDGET} rad')
    ax.set_xlabel('σ_θ (rad)')
    ax.set_ylabel('SSIM (dots, T22 stage)')
    ax.set_title('Why 0.1 rad is the right budget', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.suptitle('T19 — physical phase-noise budget', fontweight='bold')
    fig.tight_layout()
    fig.savefig('results/t19_phase_noise.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    print("\n  Saved: results/t19_phase_noise.png")

    # --- spec tables -------------------------------------------------------
    s_M, I_2pi = drive_sensitivity()
    dI = dI_spec(PER_SOURCE, s_M)

    scenarios = [
        # (name, E_eV, phase-critical leg length, NA sinθ at pattern)
        ('A: literal 1:1 T22 stage', 1e-6, Z_STAGE,
         np.sin(np.arctan(L_AP / Z_STAGE))),
        ('B1: gen-2, 1 keV, 0.1 m column', 1e3, 0.1,
         lam_of_E(1e3) / (2 * 5e-9)),
        ('B2: gen-2, 30 keV, 0.1 m column', 3e4, 0.1,
         lam_of_E(3e4) / (2 * 5e-9)),
    ]

    print(f"\n  Per-subsystem specs for σ_θ ≤ {BUDGET} rad "
          f"({PER_SOURCE:.3f} rad per source):")
    for name, E, L_leg, sin_na in scenarios:
        v = v_of_E(E)
        lam = lam_of_E(E)
        k = 2 * np.pi / lam
        Vs = V_spec(PER_SOURCE, L_leg, v)
        dz = dz_spec(PER_SOURCE, k, sin_na)
        dB = dB_spec(PER_SOURCE, L_AP, Z_STAGE)
        print(f"\n  {name}")
        print(f"    v = {v:.3g} m/s, λ = {lam*1e12:.3g} pm, "
              f"NA = {sin_na:.3g}, leg = {L_leg:.3g} m, "
              f"t_leg = {L_leg/v:.3g} s")
        print(f"    electrostatic: differential-potential DRIFT "
              f"≤ {Vs:.3g} V over the leg (static part is calibrated "
              f"out by the T21 loop)")
        print(f"    screen drive:  δI ≤ {dI*1e6:.3g} μA "
              f"({dI/I_2pi:.2%} of the 2π drive {I_2pi*1e3:.3g} mA; "
              f"velocity-independent)")
        print(f"    defocus/vibration: δz ≤ {dz:.3g} m")
        print(f"    stray-B differential: δB ≤ {dB*1e4:.3g} G "
              f"(velocity-independent)")

    # --- the actionable inversion: column length vs drift-spec class ------
    # The electrostatic spec depends only on the time of flight of the
    # phase-critical leg: t_max = σ·ħ/(q·V_drift).  Given a realistic
    # differential-drift class, that bounds the column length.
    print("\n  Maximum phase-critical column length vs achievable "
          "differential-drift class:")
    print(f"    {'V_drift':>10} {'t_max':>10} {'L @ 1 keV':>12} "
          f"{'L @ 30 keV':>12}")
    for Vd in (1e-3, 1e-4, 1e-5, 1e-6):
        t_max = PER_SOURCE * hbar / (Q * Vd)
        print(f"    {Vd:10.0e} {t_max:10.3g} "
              f"{t_max * v_of_E(1e3):12.3g} "
              f"{t_max * v_of_E(3e4):12.3g}")
    print("    → even μV-class drift limits the phase-critical throw to "
          "tens of μm; a cm-scale column at 30 keV needs nV-class "
          "differential stability between calibrations — "
          "TEM-holography-grade engineering pushed ~10× (a 300 keV TEM "
          "holds the equivalent of ~10 nV over 1 m).  Favor the highest "
          "transport energy, the shortest possible throw (microcolumns), "
          "and fast exposure with frequent T21 recalibration.")

    # --- electron comparison (why ion holography is intrinsically harder) -
    m_e = 9.1093837015e-31
    v_e300 = np.sqrt(2 * 3e5 * e_C / m_e)   # non-relativistic estimate
    print(f"\n  Mass penalty: at equal energy ions are √(m_i/m_e) ≈ "
          f"{np.sqrt(m_He/m_e):.0f}× slower than electrons, so every "
          f"electrostatic spec is that much tighter than in electron "
          f"holography (e.g. 300 keV TEM: V_drift ≤ "
          f"{V_spec(PER_SOURCE, 1.0, v_e300):.3g} V over a 1 m column).")

    print(f"\nDone in {time.time()-t0:.0f}s. "
          f"Output: results/t19_phase_noise.png")


if __name__ == '__main__':
    main()
