# Numerical Verification

The checks behind the quantitative claims in this review. All runnable with
numpy + scipy only (no torch), using the same constants and geometry as the
codebase (N = 256, L = 400 nm, T = 1 mK He⁺, z = 20λ, σ_beam = 0.35 L).

## §1 Tilt factor vs the evanescent cutoff (E1)

```python
import numpy as np
N, L = 256, 400e-9
hbar=1.0545718e-34; k_B=1.380649e-23; m_He=6.6464731e-27
v  = np.sqrt(2*k_B*1e-3/m_He);  k0 = m_He*v/hbar;  lam = 2*np.pi/k0
z  = 20*lam;  dx = L/N
x  = np.linspace(-L/2, L/2, N);  X, Y = np.meshgrid(x, x, indexing='ij')
kx = np.fft.fftfreq(N, dx)*2*np.pi
KX, KY = np.meshgrid(kx, kx, indexing='ij')
valid = (k0**2 - KX**2 - KY**2) > 0          # propagator's propagating set

def propagating_fraction(psi):
    F = np.fft.fft2(psi)
    return (np.abs(F)**2 * valid).sum() / (np.abs(F)**2).sum()

g = np.exp(-(X**2+Y**2)/(2*(0.35*L)**2))
print(propagating_fraction(g.astype(complex)))      # 1.0000
print(propagating_fraction(g*np.exp(1j*k0*X)))      # 0.8334
```

**Result:** plain Gaussian keeps 100.0% of its power in the propagating
band; the tilted beam (as built in sim_v9/sim_v10/build_beam) keeps 83.3% —
16.7% is zeroed before any phase screen. λ = 48.91 nm, k0 = 1.285×10⁸ rad/m,
z = 978.2 nm, matching the report's appendix values.

Also: the surviving spectral half spans kz ∈ (0, ~√(2k0Δk)) ≈ 4.2×10⁷ rad/m
for the beam's Δk ≈ 1/σ ≈ 7×10⁶ rad/m, so kz·z varies by ~40 rad across the
bump — the tilted beam is strongly distorted by propagation alone.

## §2 Applied vs nominal phase noise (E2)

```python
from scipy.ndimage import gaussian_filter
rng = np.random.default_rng(0)
for r in [0.9969, 0.95, 0.87]:
    noise = (1-r)*rng.normal(0, np.pi, (256, 256))
    print(r, (1-r)*np.pi, gaussian_filter(noise, sigma=3).std())
```

| r | nominal RMS | applied RMS | physically consistent σ_θ = √(−2 ln r) |
|---|---|---|---|
| 0.9969 | 0.0097 | 0.0009 | 0.0788 |
| 0.95   | 0.1571 | 0.0150 | 0.3203 |
| 0.87   | 0.4084 | 0.0387 | 0.5278 |

Smoothing factor ≈ 1/(2σ√π) ≈ 0.095 at σ = 3 px, consistent with the
measured ratios. The v10 report's §6.2 "0.42 rad RMS" case ran at 0.039 rad.

## §3 Space charge, cyclotron radius, Langevin mfp (M1, M3, E6)

```python
I, e, v, m, B = 1e-6, 1.602e-19, 2.038, 6.646e-27, 0.01
lam_line = I/v                       # 4.91e-7 C/m  →  3.06e12 ions/m
print(lam_line/e * 978e-9)           # 3.0e6 ions in the 978 nm flight path
r_b = 140e-9
n = I/(e*v*np.pi*r_b**2)             # 4.97e25 m^-3 beam density
d = n**(-1/3)                        # 2.72 nm mean spacing
print(2.307e-28/d/e)                 # 0.53 eV Coulomb energy per ion
                                     # vs k_B·1mK = 8.6e-8 eV kinetic
print(m*v/(e*B))                     # 8.5 um cyclotron radius vs 100 um AK gap
k_L = 2.3e-15                        # typical ion–neutral Langevin rate m^3/s
print(k_L/v)                         # 1.1e-15 m^2 cross-section
n_gas = 101325/(1.381e-23*300)
print(1/(n_gas*k_L/v))               # 3.6e-11 m = 0.036 nm mfp at 1 atm
```

**Results:** ~3×10⁶ ions in flight simultaneously; Coulomb/kinetic ratio
≈ 6×10⁶; cyclotron radius 8.5 μm ≪ 100 μm gap; Langevin cross-section
1.1×10⁻¹⁵ m² vs the code's 3×10⁻²⁴ m² default (~8.6 orders); mfp at 1 atm
0.036 nm vs the code's mm-scale claim.

Single-ion-at-a-time current bound: transit ≈ 0.5 μs over the 978 nm leg
(longer including approach), so ≤ ~2×10⁶ ions/s ≈ 0.3 pA keeps occupancy ≤ 1.

## §4 Bandwidth bookkeeping (M5)

```python
pitch = L/32
print(2*np.pi/pitch / k0)     # 3.91  — pitch spatial frequency vs k0
print(k0*L/(2*np.pi))         # 8.2   — propagating cycles/aperture (field)
# array Nyquist = 32/2 = 16 cycles/aperture
# evanescent decay length at pitch frequency:
kxp = 2*np.pi/pitch
print(1/np.sqrt(kxp**2-k0**2))  # ≈ 2.05e-9 m
```

**Results:** field-level free-space cutoff 8.2 cycles/aperture < array
Nyquist 16; pitch-scale phase features decay in ~2 nm, i.e. are invisible at
z = 978 nm. Intensity (an autocorrelation of the field spectrum) can in
principle contain up to 2×8.2 ≈ 16.4 c/a, but only via interference of
near-±90° components — at severe efficiency cost, so the practical intensity
band is well below that.

## §5 Kuramoto critical coupling (E5)

For Gaussian g(ω) with std σ: K_c = 2/(π g(0)) = σ·√(8/π) ≈ **1.5958 σ**.
The code uses K_c = 2σ (the Lorentzian value, where γ is the HWHM), together
with r = √(1 − K_c/K), which is exact only for the Lorentzian. The two modes
(ODE: Gaussian + inertia; analytic: Lorentzian-form, first-order) therefore
disagree by construction; neither has been validated against the other.

## §6 Diamond Hamiltonian spot-check (E7 — confirming correctness)

Hopping matrix elements in `build_hamiltonian` (phase convention
`phase_p = J·e^{+iφ/2}` on B1→A_right, `phase_m = J·e^{−iφ/2}` on B2→A_right):
traversing A→B1→A′→B2→A multiplies ⟨B1|H|A⟩·⟨A′|H|B1⟩·⟨B2|H|A′⟩·⟨A|H|B2⟩
= J · Je^{−iφ/2} · Je^{−iφ/2} · J → net phase −φ: flux magnitude φ per
plaquette as intended (sign is gauge/orientation). Each A hub carries 8
bonds (B1,B2,C1,C2 of its own cell + B pair of the left cell + C pair of the
lower cell), so the caged eigenvalues are ±√8·J = ±2√2·J plus the flat E = 0
manifold — exactly the {−2√2, 0, +2√2}J spectrum the validation gate checks.
Only the docstring ("coordination 4") is wrong.
