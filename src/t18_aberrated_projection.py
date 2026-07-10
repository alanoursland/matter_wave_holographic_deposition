"""
T18 — Aberrated Projection Stage
=================================

Replaces the ideal free-space stage with a wave-optical model of the
gen-2 projection column's aberrations, evaluated on the T22-recommended
pattern (17.6 c/a over a 400 nm field ≈ 11 nm features).

Model
-----
The column delivers the demagnified image of the phase plate to the
substrate; the image-side field just before landing carries the T22
pattern's transverse spectrum on a wave at the LANDING energy E_land
(the final deceleration to soft-landing energy is part of the imaging
optics, LEEM-style).  Aberrations multiply the angular spectrum by
exp(iW):

    defocus     W = Δf · k⊥²/(2k)
    spherical   W = C_s · k⊥⁴/(4k³)
    chromatic   Δf(δE) = C_c · δE/E_land, intensities averaged over the
                source energy distribution (different energies add
                incoherently in the time-averaged arrival pattern)

Angles at the substrate are tiny (θ_max = k⊥_max/k_land ~ 10⁻³–10⁻⁴ rad
because λ_land is picometers), so third-order axial theory is adequate;
off-axis (field) aberrations are NOT modeled — see caveats.

Static aberrations (defocus, spherical, even a calibrated chromatic
mean) are hologram-precompensatable in the T21 loop; the irreducible
term is chromatic SPREAD.  The sweep therefore maps SSIM vs
(E_land, C_c, ΔE), with a C_s sensitivity check.

Validation gate: the wave model must reproduce the known linear scaling
of the chromatic blur disc d_c = C_c·θ·(ΔE/E) of charged-particle
columns.

Run:  python t18_aberrated_projection.py
"""

import os
import time
import numpy as np
import torch
import matplotlib.pyplot as plt

from inverse_holography import SQUIDArray, InverseHolographySolver, \
    target_grid_of_dots, smooth_target
from coherent_matterwave_beam import hbar, e_C
from iqs.constants import k_B, m_He
from iqs.numerics.metrics import ssim_score
from iqs.numerics.device import get_device

_device = get_device()
os.makedirs('results', exist_ok=True)

# T22-recommended stage (pattern definition)
N_GRID, L_AP, N_LOOPS = 256, 400e-9, 64
LAM_STAGE, Z_STAGE = 14.4e-9, 489e-9

H_PLANCK = 6.62607015e-34


def lam_of_E(E_eV, m=m_He):
    return H_PLANCK / np.sqrt(2 * m * E_eV * e_C)


# ---------------------------------------------------------------------------
# Ideal image field (the T22 pattern, solved once)
# ---------------------------------------------------------------------------

def ideal_image_field(seed=0):
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
    screen_t = torch.tensor(sol['phase_screen'], dtype=torch.float64,
                            device=_device)
    T = SQUIDArray.phase_screen_to_transmission(screen_t)
    psi_out = solver._propagate_torch(solver.psi_in_t * T)
    return psi_out.cpu().numpy(), cond


# ---------------------------------------------------------------------------
# Aberration transfer
# ---------------------------------------------------------------------------

def _k_grids(N, L):
    k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    KX, KY = np.meshgrid(k, k, indexing='ij')
    return KX**2 + KY**2


def aberrated_intensity(psi, E_land_eV, C_c, dE_fwhm_eV, C_s=0.0,
                        defocus=0.0, n_energy=21, L=L_AP):
    """Time-averaged arrival intensity through the aberrated column.

    Chromatic treatment: incoherent average of |image|² over a Gaussian
    energy distribution (FWHM dE_fwhm_eV), each energy seeing a defocus
    C_c·δE/E_land.
    """
    N = psi.shape[0]
    K2 = _k_grids(N, L)
    k_land = 2 * np.pi / lam_of_E(E_land_eV)
    PSI = np.fft.fft2(psi)

    # static terms
    W_static = C_s * K2**2 / (4 * k_land**3) + defocus * K2 / (2 * k_land)

    if dE_fwhm_eV > 0:
        sig = dE_fwhm_eV / 2.355
        dEs = np.linspace(-3 * sig, 3 * sig, n_energy)
        w = np.exp(-dEs**2 / (2 * sig**2))
        w /= w.sum()
    else:
        dEs, w = np.array([0.0]), np.array([1.0])

    I = np.zeros((N, N))
    for dE, wj in zip(dEs, w):
        df = C_c * dE / E_land_eV
        W = W_static + df * K2 / (2 * k_land)
        img = np.fft.ifft2(PSI * np.exp(1j * W))
        I += wj * np.abs(img)**2
    return I


# ---------------------------------------------------------------------------
# Validation gate: chromatic blur scaling of a focused spot
# ---------------------------------------------------------------------------

def _spot_d50(I, L=L_AP):
    """FW50 diameter (disc containing 50% of the intensity) of a spot.

    The standard probe-size measure in charged-particle optics: with a
    Gaussian energy spread the chromatic PSF keeps a sharp
    diffraction-limited core (the near-zero-δE components) and puts the
    blur in the tails, so FWHM under-reports it badly.
    """
    N = I.shape[0]
    Is = np.fft.fftshift(I)
    i0, j0 = np.unravel_index(np.argmax(Is), Is.shape)
    y, x = np.ogrid[:N, :N]
    r = np.sqrt((y - i0)**2 + (x - j0)**2).ravel()
    vals = Is.ravel()
    order = np.argsort(r)
    csum = np.cumsum(vals[order])
    r50 = r[order][np.searchsorted(csum, 0.5 * csum[-1])]
    return 2 * r50 * (L / N)


def validate_chromatic_scaling(E_land_eV=1.0, verbose=True):
    """d_c must scale linearly in C_c (and ΔE): double C_c → ~double
    blur of a focused spot, in the chromatic-dominated regime.
    Parameters chosen so the blur is ≫ the ideal spot (~11 nm) but well
    inside the periodic frame (400 nm)."""
    K2 = _k_grids(N_GRID, L_AP)
    k_land = 2 * np.pi / lam_of_E(E_land_eV)
    k_na = 2 * np.pi * 17.6 / L_AP          # stage bandwidth
    disc = (K2 <= k_na**2).astype(complex)  # ideal converging aperture
    psi_spot = np.fft.ifft2(disc)

    dE = 0.1                                 # eV FWHM
    # diffraction core (no chromatic), then two chromatic points; the
    # chromatic blur is the quadrature excess over the core
    d50_core = _spot_d50(
        aberrated_intensity(psi_spot, E_land_eV, 0.0, dE_fwhm_eV=0.0))
    blurs = []
    for C_c in (0.5e-3, 1.0e-3):
        d50 = _spot_d50(
            aberrated_intensity(psi_spot, E_land_eV, C_c, dE_fwhm_eV=dE))
        blurs.append(np.sqrt(max(d50**2 - d50_core**2, 0.0)))
    ratio = blurs[1] / blurs[0]
    ok_scaling = 1.6 < ratio < 2.4
    # Barth–Kruit FW50 convention for a Gaussian energy distribution:
    # d_c(FW50) ≈ 0.34 · C_c · α · ΔE_FWHM/E
    theta = k_na / k_land
    d_bk = 0.34 * 0.5e-3 * theta * dE / E_land_eV
    ok_coef = abs(blurs[0] / d_bk - 1) < 0.4
    ok = ok_scaling and ok_coef
    if verbose:
        print(f"  CHROMATIC-SCALING GATE: blur(1 mm)/blur(0.5 mm) "
              f"= {ratio:.2f} (expect ≈2); "
              f"blur = {blurs[0]*1e9:.1f} nm vs Barth–Kruit FW50 "
              f"{d_bk*1e9:.1f} nm  [{'OK' if ok else 'FAIL'}]")
    return ok, blurs[0], d_bk


# ---------------------------------------------------------------------------
# Study
# ---------------------------------------------------------------------------

E_LANDS = [1.0, 10.0]                       # eV
C_CS    = [1e-4, 1e-3, 1e-2]                # m
DES     = [0.5, 0.1, 0.01, 0.001]           # eV FWHM


def main():
    print("=" * 70)
    print("T18: ABERRATED PROJECTION STAGE (wave-optical, image side)")
    print(f"     pattern: T22 stage (17.6 c/a over {L_AP*1e9:.0f} nm "
          f"≈ 11 nm features)")
    print("=" * 70)
    t0 = time.time()

    ok, _, _ = validate_chromatic_scaling()
    if not ok:
        raise RuntimeError("Chromatic scaling gate failed")

    psi, cond = ideal_image_field()
    ssim_ideal = float(ssim_score(cond, np.abs(psi)**2))
    print(f"\n  Ideal image SSIM (no aberrations): {ssim_ideal:.3f}")

    # C_s sensitivity check at the gen-2 NA
    for E_land in E_LANDS:
        theta = (2 * np.pi * 17.6 / L_AP) / (2 * np.pi / lam_of_E(E_land))
        I_cs = aberrated_intensity(psi, E_land, C_c=0.0, dE_fwhm_eV=0.0,
                                   C_s=1.0)
        s = float(ssim_score(cond, I_cs))
        print(f"  C_s = 1 m at E_land = {E_land:.0f} eV "
              f"(θ_max = {theta*1e3:.2f} mrad): SSIM = {s:.3f} "
              f"(ideal {ssim_ideal:.3f}) — spherical "
              f"{'negligible' if ssim_ideal - s < 0.005 else 'MATTERS'}")

    # chromatic sweep
    print(f"\n  Chromatic sweep (SSIM; analytic blur d_c in nm):")
    print(f"  {'E_land':>7} {'C_c':>8}" +
          "".join(f"  ΔE={d:>5g} eV" for d in DES))
    results = {}
    for E_land in E_LANDS:
        k_land = 2 * np.pi / lam_of_E(E_land)
        theta = (2 * np.pi * 17.6 / L_AP) / k_land
        for C_c in C_CS:
            cells = []
            for dE in DES:
                I = aberrated_intensity(psi, E_land, C_c, dE)
                s = float(ssim_score(cond, I))
                d_c = C_c * theta * dE / E_land
                results[(E_land, C_c, dE)] = (s, d_c)
                cells.append(f"  {s:.3f}/{d_c*1e9:6.1f}")
            print(f"  {E_land:7.0f} {C_c*1e3:6.1f}mm" + "".join(cells))

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, E_land in zip(axes, E_LANDS):
        for C_c in C_CS:
            ss = [results[(E_land, C_c, d)][0] for d in DES]
            ax.semilogx(DES, ss, 'o-', lw=2, label=f'C_c = {C_c*1e3:g} mm')
        ax.axhline(ssim_ideal, color='k', ls=':', alpha=0.6,
                   label='no aberration')
        ax.axhline(0.8, color='g', ls='--', alpha=0.5, label='SSIM 0.8')
        ax.set_xlabel('source energy spread ΔE FWHM (eV)')
        ax.set_ylabel('SSIM at substrate')
        ax.set_title(f'E_land = {E_land:.0f} eV', fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1)
    fig.suptitle('T18 — chromatic aberration of the decel/projection '
                 'stage (T22 pattern)', fontweight='bold')
    fig.tight_layout()
    fig.savefig('results/t18_aberrated_projection.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    print(f"\n  Saved: results/t18_aberrated_projection.png")

    # feasible-set report
    print("\n  Feasible set (SSIM ≥ 0.8):")
    for (E_land, C_c, dE), (s, d_c) in sorted(results.items()):
        if s >= 0.8:
            print(f"    E_land = {E_land:4.0f} eV, C_c = {C_c*1e3:5.1f} mm, "
                  f"ΔE ≤ {dE:g} eV   (d_c = {d_c*1e9:.2f} nm, "
                  f"SSIM = {s:.3f})")

    print(f"\nDone in {time.time()-t0:.0f}s. "
          f"Output: results/t18_aberrated_projection.png")


if __name__ == '__main__':
    main()
