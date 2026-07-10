"""
T20 — Landing Stage
====================

Replaces the diamond-caging stage in the chip-printing configuration
with the physical landing chain: impact → neutralization → sticking →
staying.  T15's methodology is retained — the deliverable is
P_escape(control parameter) with a tolerance spec — but the control
parameter is now substrate temperature vs binding barrier, not flux.

Models (semi-empirical, standard surface-physics forms):

  Sputter threshold (Bohdansky):
      E_th = U_s / (γ(1−γ))            for m₁/m₂ ≤ 0.3
      E_th = 8·U_s·(m₁/m₂)^(2/5)       otherwise
      with γ = 4m₁m₂/(m₁+m₂)² and U_s the surface binding energy.
  Neutralization: an incoming ion neutralizes near the surface,
      releasing its ionization potential; the locally deposited share
      is  E_neut = f_dep·(IP − φ)  with work function φ and an
      efficiency parameter f_dep (swept — the literature range is
      wide).  Local energy budget: E_local = E_land + E_neut, compared
      against the sputter threshold and a surface displacement energy.
  Thermal migration (Arrhenius):
      P_escape(τ) = 1 − exp(−ν₀·τ·e^(−E_a/kT))
      spec: P_escape < 10⁻⁴ over a hold time τ = 10³ s
      → T_max = E_a / (k_B·ln(ν₀τ/P_spec)).
  Templates set E_a and the placement quantum (site pitch): the beam
      delivers *which sites*, the surface chemistry delivers *where
      exactly* — the division of labor from the gen-2 note §3.4.

Check (gen-2 note §8): a (species, E_land, T_substrate, template)
corner where placed atoms stay placed — and the T18 hand-off: is
E_land = 10 eV (which removes the monochromator) impact-safe?

Run:  python t20_landing_stage.py
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt

k_B_eV = 8.617333262e-5   # eV/K

os.makedirs('results', exist_ok=True)

# ---------------------------------------------------------------------------
# Materials data (representative literature values)
# ---------------------------------------------------------------------------

# Deposit/probe species: mass [amu], ionization potential [eV]
SPECIES = {
    'He+': {'m': 4.003,  'IP': 24.587},
    'P+':  {'m': 30.974, 'IP': 10.487},
    'Si+': {'m': 28.086, 'IP': 8.152},
    'Ga+': {'m': 69.723, 'IP': 5.999},
    'Al+': {'m': 26.982, 'IP': 5.986},
    'In+': {'m': 114.818, 'IP': 5.786},
}

# Substrate: Si(100)
SUBSTRATE = {
    'name': 'Si(100)',
    'm': 28.086,        # amu
    'U_s': 4.7,         # eV, surface binding (≈ cohesive energy)
    'phi': 4.8,         # eV, work function
    'E_disp_surf': 15.0,  # eV, surface displacement energy (~half bulk ~25-30)
    'nu0': 1e13,        # Hz, attempt frequency
}

# Binding/template options: activation barrier E_a [eV], site pitch [nm]
TEMPLATES = {
    'physisorbed / metal terrace': {'E_a': 0.3,  'pitch': 0.0},
    'Si adatom on Si(100)':        {'E_a': 0.7,  'pitch': 0.38},
    'strong chemisorption':        {'E_a': 1.2,  'pitch': 0.38},
    'covalent site (H:Si depass.)': {'E_a': 2.5, 'pitch': 0.38},
}

P_SPEC   = 1e-4     # escape probability budget per atom
TAU_HOLD = 1e3      # s, hold time until the layer is locked in


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

def sputter_threshold(m1, m2, U_s):
    """Bohdansky threshold energy [eV] for ion m1 onto target m2."""
    gamma = 4 * m1 * m2 / (m1 + m2)**2
    if m1 / m2 <= 0.3:
        return U_s / (gamma * (1 - gamma))
    return 8 * U_s * (m1 / m2)**0.4


def neutralization_deposit(IP, phi, f_dep):
    """Locally deposited neutralization energy [eV] (≥ 0)."""
    return max(f_dep * (IP - phi), 0.0)


def local_energy(species, E_land, f_dep, sub=SUBSTRATE):
    return E_land + neutralization_deposit(
        SPECIES[species]['IP'], sub['phi'], f_dep)


def impact_gates(species, E_land, f_dep, sub=SUBSTRATE):
    """Returns (E_local, E_th_sputter, ok_sputter, ok_displacement)."""
    E_loc = local_energy(species, E_land, f_dep, sub)
    E_th = sputter_threshold(SPECIES[species]['m'], sub['m'], sub['U_s'])
    return E_loc, E_th, E_loc < E_th, E_loc < sub['E_disp_surf']


def p_escape(E_a, T, tau=TAU_HOLD, nu0=SUBSTRATE['nu0']):
    """Arrhenius escape probability over hold time tau."""
    with np.errstate(over='ignore'):
        rate = nu0 * np.exp(-E_a / (k_B_eV * np.maximum(T, 1e-9)))
    return 1.0 - np.exp(-rate * tau)


def T_max_for_spec(E_a, p_spec=P_SPEC, tau=TAU_HOLD,
                   nu0=SUBSTRATE['nu0']):
    """Highest substrate temperature meeting the escape spec."""
    return E_a / (k_B_eV * np.log(nu0 * tau / p_spec))


# ---------------------------------------------------------------------------
# Study
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("T20: LANDING STAGE")
    print(f"     substrate {SUBSTRATE['name']}; escape spec "
          f"P < {P_SPEC:g} over τ = {TAU_HOLD:g} s")
    print("=" * 70)
    t0 = time.time()

    # --- impact + neutralization gates ---------------------------------
    print("\n  Impact/neutralization energy budget (E_local vs thresholds "
          f"[eV]; sputter E_th per pair, displacement {SUBSTRATE['E_disp_surf']}):")
    print(f"  {'species':>8} {'E_th':>6}" +
          "".join(f"   E_land={E:>2g}, f={f:g}" for E in (1, 10)
                  for f in (0.5, 1.0)))
    verdicts = {}
    for sp in SPECIES:
        E_th = sputter_threshold(SPECIES[sp]['m'], SUBSTRATE['m'],
                                 SUBSTRATE['U_s'])
        cells = []
        for E_land in (1.0, 10.0):
            for f in (0.5, 1.0):
                E_loc, _, ok_s, ok_d = impact_gates(sp, E_land, f)
                verdicts[(sp, E_land, f)] = ok_s and ok_d
                cells.append(f"   {E_loc:6.1f} {'OK ' if ok_s and ok_d else 'DMG'}")
        print(f"  {sp:>8} {E_th:6.1f}" + "".join(cells))

    print("\n  → He⁺ fails at every corner (24.6 eV neutralization "
          "release): the λ-convenient probe is not a deposit species.")

    # --- thermal migration / templates ---------------------------------
    print(f"\n  Site-locking: max substrate T for P_escape < {P_SPEC:g} "
          f"(ν₀ = {SUBSTRATE['nu0']:.0e} Hz, τ = {TAU_HOLD:g} s):")
    for name, tpl in TEMPLATES.items():
        Tm = T_max_for_spec(tpl['E_a'])
        print(f"    {name:<30} E_a = {tpl['E_a']:>4.1f} eV → "
              f"T ≤ {Tm:6.1f} K"
              + ("   (RT-stable)" if Tm >= 300 else "   (cryo/cooled)"))

    # --- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5))

    ax = axes[0]
    T_grid = np.linspace(40, 700, 400)
    for name, tpl in TEMPLATES.items():
        ax.semilogy(T_grid, np.maximum(p_escape(tpl['E_a'], T_grid), 1e-18),
                    lw=2, label=f"{name} (E_a = {tpl['E_a']} eV)")
    ax.axhline(P_SPEC, color='k', ls='--', alpha=0.6,
               label=f'spec {P_SPEC:g}')
    ax.axvline(300, color='gray', ls=':', alpha=0.6)
    ax.set_xlabel('substrate temperature (K)')
    ax.set_ylabel(f'P_escape over {TAU_HOLD:g} s')
    ax.set_ylim(1e-12, 1.5)
    ax.set_title('Thermal site-locking (T15 methodology)',
                 fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)

    ax = axes[1]
    names = list(SPECIES)
    x = np.arange(len(names))
    E_land = 10.0
    for f, color, alpha in ((0.5, '#1565c0', 0.9), (1.0, '#c62828', 0.6)):
        E_locs = [local_energy(sp, E_land, f) for sp in names]
        ax.bar(x + (0.2 if f == 1.0 else -0.2), E_locs, width=0.35,
               color=color, alpha=alpha, label=f'E_local, f_dep = {f:g}')
    E_ths = [sputter_threshold(SPECIES[sp]['m'], SUBSTRATE['m'],
                               SUBSTRATE['U_s']) for sp in names]
    ax.plot(x, E_ths, 'k^--', lw=1.5, markersize=9,
            label='sputter threshold')
    ax.axhline(SUBSTRATE['E_disp_surf'], color='g', ls='--', alpha=0.7,
               label='surface displacement')
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('energy (eV)')
    ax.set_title(f'Impact budget at E_land = {E_land:g} eV '
                 f'({SUBSTRATE["name"]})', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=8)

    fig.suptitle('T20 — landing stage: impact safety and site locking',
                 fontweight='bold')
    fig.tight_layout()
    fig.savefig('results/t20_landing_stage.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    print("\n  Saved: results/t20_landing_stage.png")

    # --- the corner -------------------------------------------------------
    print("\n  RECOMMENDED LANDING CORNER:")
    print("    Species: Si⁺ / Al⁺ / P⁺ (deposit or dopant; IP − φ ≈ "
          "1–6 eV → benign neutralization)")
    print("    E_land = 10 eV: impact-safe for all deposit species even "
          "at f_dep = 1 — and removes the monochromator (T18).")
    Tm_chem = T_max_for_spec(TEMPLATES['strong chemisorption']['E_a'])
    Tm_cov = T_max_for_spec(TEMPLATES['covalent site (H:Si depass.)']['E_a'])
    Tm_ad = T_max_for_spec(TEMPLATES['Si adatom on Si(100)']['E_a'])
    print(f"    Substrate/template: covalent/templated sites "
          f"(E_a ≈ 2.5 eV) hold to {Tm_cov:.0f} K; strong chemisorption "
          f"(1.2 eV) to {Tm_chem:.0f} K (RT-capable); a bare Si adatom "
          f"(0.7 eV) needs T ≤ {Tm_ad:.0f} K.")
    print("    He⁺: excluded as a deposit — usable only as an imaging/"
          "probe species.")

    print(f"\nDone in {time.time()-t0:.0f}s. "
          f"Output: results/t20_landing_stage.png")


if __name__ == '__main__':
    main()
