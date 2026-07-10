# T20 — Landing Stage: Results

Executed 2026-07-09. Study: [src/t20_landing_stage.py](../src/t20_landing_stage.py).
Log: [t20_landing.txt](t20_landing.txt). Figure:
`results/t20_landing_stage.png`. 113 tests pass (7 new in `test_t20.py`).

## Model

Semi-empirical landing chain replacing diamond caging in the
chip-printing configuration (T15 methodology retained — P_escape vs
control parameter with a tolerance spec):

- **Impact:** Bohdansky sputter threshold E_th(m₁, m₂, U_s); surface
  displacement energy 15 eV for Si(100).
- **Neutralization:** locally deposited E = f_dep·(IP − φ), f_dep swept
  0.5–1 (literature range is wide; f_dep = 1 is worst case).
- **Thermal site-locking:** Arrhenius escape over a τ = 10³ s hold,
  spec P_escape < 10⁻⁴ → T_max = E_a/(k_B·ln(ν₀τ/P)).
- **Templates** set E_a and the placement quantum (0.38 nm on Si(100)):
  beam chooses *which sites*, chemistry chooses *where exactly*.

## Results

**Impact/neutralization budget on Si(100) (E_local vs thresholds, eV):**

| species | E_th sputter | E_land = 10 eV, f_dep = 1 | verdict |
|---|---|---|---|
| He⁺ | 19.1 | 29.8 | **damage at every corner ≥ f=1 or 10 eV** |
| P⁺ | 39.1 | 15.7 | marginal at worst case (vs 15 eV displacement) |
| Si⁺ | 37.6 | 13.4 | safe |
| Al⁺ | 37.0 | 11.2 | safe |
| Ga⁺ | 54.1 | 11.2 | safe |
| In⁺ | 66.0 | 11.0 | safe |

**Site-locking temperatures (P_escape < 10⁻⁴ over 10³ s):**

| binding | E_a | T_max |
|---|---|---|
| physisorbed / metal terrace | 0.3 eV | 76 K |
| Si adatom on Si(100) | 0.7 eV | 176 K |
| strong chemisorption | 1.2 eV | 302 K — RT threshold |
| covalent / templated (H:Si) | 2.5 eV | 630 K |

## The corner (gen-2 acceptance check met)

**Si⁺/Al⁺/Ga⁺ (deposits) or P⁺ (dopant, keep f_dep margin) at
E_land = 10 eV onto templated/chemisorbing Si at ≤ 300 K.**

- E_land = 10 eV is impact-safe for every deposit species even at
  worst-case neutralization — so **T18's monochromator-free operating
  point is confirmed landing-safe**. The E_land trade closes: 10 eV is
  right on both sides.
- RT operation requires E_a ≥ 1.2 eV (strong chemisorption or covalent
  template); bare-adatom deposition needs a cooled stage (≤ 176 K).
  The quantitative version of the gen-2 §3.4 claim.
- **He⁺ is formally excluded as a deposit** (24.6 eV neutralization
  release exceeds the sputter threshold of its own landing site) — the
  v10 species was always a λ-placeholder; now it's a closed question.
  He⁺ remains fine as a probe/imaging species.

## Caveats

- f_dep and E_a values are representative literature-class numbers, not
  measurements for a specific surface prep; the framework inverts
  cleanly when real numbers exist.
- No modeling of transient local heating kinetics (the eV deposit
  thermalizes over a few atoms in ps — treated as a lump energy
  comparison), reflection/sticking coefficients (≈1 for chemisorbing
  species at these energies), or chemistry between deposit and
  template beyond the barrier parameter.
- Placement precision inherits the template pitch (0.38 nm), consistent
  with the stochastic-plus-quantization division of labor; edge
  statistics remain T21's department.
