"""
Physical constants used across the iqs library.

All values are in SI units unless noted otherwise.
"""

# ── Quantum mechanics ─────────────────────────────────────────────────────────
hbar = 1.0545718e-34    # reduced Planck constant  [J·s]
h    = 6.62607015e-34   # Planck constant          [J·s]

# ── Thermodynamics ────────────────────────────────────────────────────────────
k_B  = 1.380649e-23     # Boltzmann constant       [J/K]

# ── Electromagnetism ─────────────────────────────────────────────────────────
e_C  = 1.602176634e-19  # elementary charge        [C]
mu_0 = 4e-7 * 3.14159265358979323846  # vacuum permeability  [H/m]
Phi_0     = 2.067833848e-15  # magnetic flux quantum  h/(2e)  [Wb]
Phi_0_e   = 4.135667696e-15  # single-charge flux quantum  h/e  [Wb]
c    = 2.99792458e8     # speed of light           [m/s]

# ── Particle masses ───────────────────────────────────────────────────────────
m_e  = 9.1093837015e-31  # electron mass            [kg]
m_u  = 1.66053906660e-27 # atomic mass unit         [kg]
m_He = 6.6464731e-27     # helium-4 atomic mass     [kg]  (= 4.0026 * m_u)
