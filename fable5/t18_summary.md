# T18 — Aberrated Projection Stage: Results

Executed 2026-07-09. Study: [src/t18_aberrated_projection.py](../src/t18_aberrated_projection.py).
Log: [t18_aberration.txt](t18_aberration.txt). Figure:
`results/t18_aberrated_projection.png`. 106 tests pass (5 new in
`test_t18.py`).

## Model and validation

Wave-optical, image-side: the T22 pattern's field (17.6 c/a over
400 nm ≈ 11 nm features) is carried by a wave at the *landing* energy;
column aberrations multiply its angular spectrum by exp(iW) with
defocus, spherical (C_s·k⊥⁴/4k³), and chromatic terms — chromatic as an
incoherent average over a Gaussian source-energy distribution, each
energy seeing defocus C_c·δE/E_land.

**Validation gate passed:** the model reproduces the standard
charged-particle chromatic blur — linear scaling in C_c (ratio 2.24 for
2× C_c) and the Barth–Kruit FW50 coefficient d₅₀ ≈ 0.34·C_c·α·ΔE/E
within 12%. (Instructive detour: FWHM is the wrong metric — a Gaussian
energy spread keeps a sharp diffraction core and hides the blur in the
tails, which is exactly why probe optics uses FW50.)

## Results

- **Spherical aberration is irrelevant.** Because λ_land is picometers,
  the image-side NA for 11 nm features is sub-milliradian; even
  C_s = 1 m changes nothing (θ⁴ ~ 10⁻¹³). Only chromatic and defocus
  matter — and defocus is static (T21-calibratable; tolerance already
  in T19).
- **Chromatic sweep (SSIM at substrate; ideal 0.897):**

  | E_land | C_c | ΔE = 0.5 eV | 0.1 eV | 0.01 eV |
  |---|---|---|---|---|
  | 1 eV | 0.1 mm | 0.886 | 0.897 | 0.897 |
  | 1 eV | 1 mm | 0.236 | 0.834 | 0.897 |
  | 1 eV | 10 mm | 0.112 | 0.135 | 0.834 |
  | 10 eV | 1 mm | **0.896** | 0.897 | 0.897 |
  | 10 eV | 10 mm | 0.731 | 0.894 | 0.897 |

- **SSIM is unharmed until d_c ≈ the 11 nm feature scale**, exactly as
  the analytic blur column predicts — the wave model and the standard
  formulas agree throughout.

## Verdict: aberrations do not kill gen-2

The feasible set is broad. Two comfortable operating classes:

1. **E_land = 10 eV, C_c ≤ 1 mm: works with the RAW source spread**
   (0.5 eV GFIS-class, no monochromator; SSIM 0.896). Physics:
   d_c = C_c·θ·ΔE/E_land ∝ E_land^{−3/2} — landing 10× hotter buys
   ~32× less blur. The cost is landing energy (sputtering/displacement
   risk — T20's question).
2. **E_land = 1 eV needs either a short-gap objective (C_c ≈ 0.1 mm)
   or monochromation to ΔE ≲ 10 meV** (affordable in single-ion mode
   where flux is nearly free).

Since cathode-lens C_c scales with the deceleration gap, C_c ≤ 1 mm ⇔
mm-scale final gap — **the third independent requirement to converge on
the microcolumn** (after throughput parallelism §5 and the T19
phase-stability throw limit). Short columns are simultaneously the
throughput, stability, and aberration answer.

Combined with T19: the binding constraint of the gen-2 machine remains
**electrostatic phase drift**, not aberrations.

## Caveats

- Axial model: off-axis (field) aberrations — coma, field curvature,
  distortion across a die-scale field — are not modeled; at sub-mrad NA
  over a 400 nm field they are negligible, but field stitching across a
  full die is real engineering left to hardware design.
- C_c is a parameter, not derived from a lens design; the mm-scale
  values assumed for short decel gaps are typical of immersion
  objectives but need a real column model to confirm.
- No Coulomb blur in-column — exact, not approximate, in single-ion
  operation (T10).
