"""
Coherent Matterwave Beam Model
===================================

Isolated simulation of the coherent matterwave beam generation system
described in US Patent 9,502,202 B2 (Lockheed Martin, 2016).

Models a microscopic diode cavity array that produces a coherent beam
of charged particles via Aharonov-Bohm phase synchronization.  The AB
effect shifts each particle's quantum phase through exposure to a
magnetic vector potential — without exchanging energy, preserving the
monochromaticity required for coherence.

Phase synchronization follows the Kuramoto coupled-oscillator model.
For a monochromatic beam the critical coupling threshold is zero: the
system always synchronizes regardless of coupling strength, though
stronger coupling (higher B-field) synchronizes faster.

Model choice (fable5 E5/T9): both synchronization paths implement the
*first-order* Kuramoto model with a *Lorentzian* natural-frequency
distribution of half-width γ = omega_spread_dim.  For that combination
the steady state is exact and closed-form: K_c = 2γ and
r = √(1 − K_c/K) above threshold, so the analytic and ODE modes agree
by construction (validated by `validate_kuramoto_modes`).

The simulation supports arbitrary ionized atomic species.  Mass and
charge enter through the de Broglie wavelength, velocity at fixed
kinetic energy, AB phase prefactor, and synchronization timescale.

Reference: notes/sources/coherent_matterwave_beam_patent.md
"""

import numpy as np

from iqs.constants import c, e_C, eps_0, hbar, k_B, m_e, m_u, mu_0, Phi_0_e

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
# Langevin (polarization-capture) rate constant for slow-ion–neutral
# scattering, k_L = 2πq√(α/μ)/(4πε₀)^(1/2).  ~2×10⁻¹⁵ m³/s is typical
# for ion–atom systems (fable5 E6/T8).  The velocity-independent rate
# ν = n·k_L means the effective cross-section σ_L = k_L/v grows as the
# beam slows — slow ions are enormous scattering targets.
K_LANGEVIN = 2e-15  # m³/s


# ---------------------------------------------------------------------------
# Species library
# ---------------------------------------------------------------------------
SPECIES = {
    'electron': {'mass': m_e,          'charge': 1, 'symbol': 'e⁻',   'label': 'Electron'},
    'He+':      {'mass': 4.0026*m_u,   'charge': 1, 'symbol': 'He⁺',  'label': 'Helium ion'},
    'Li+':      {'mass': 6.941*m_u,    'charge': 1, 'symbol': 'Li⁺',  'label': 'Lithium ion'},
    'Na+':      {'mass': 22.990*m_u,   'charge': 1, 'symbol': 'Na⁺',  'label': 'Sodium ion'},
    'Rb+':      {'mass': 86.909*m_u,   'charge': 1, 'symbol': 'Rb⁺',  'label': 'Rubidium ion'},
    'Cs+':      {'mass': 132.905*m_u,  'charge': 1, 'symbol': 'Cs⁺',  'label': 'Cesium ion'},
    'Ca+':      {'mass': 40.078*m_u,   'charge': 1, 'symbol': 'Ca⁺',  'label': 'Calcium ion'},
    'Sr+':      {'mass': 87.62*m_u,    'charge': 1, 'symbol': 'Sr⁺',  'label': 'Strontium ion'},
    'Yb+':      {'mass': 173.045*m_u,  'charge': 1, 'symbol': 'Yb⁺',  'label': 'Ytterbium ion'},
}


def add_species(name, mass_kg, charge_e, symbol=None, label=None):
    """Register a custom species.  charge_e in units of elementary charge."""
    SPECIES[name] = {
        'mass':   mass_kg,
        'charge': charge_e,
        'symbol': symbol or name,
        'label':  label or name,
    }


# ---------------------------------------------------------------------------
# Partial-coherence noise model (fable5 E2/T6)
# ---------------------------------------------------------------------------
def coherence_sigma_theta(r):
    """Phase-noise RMS from the order parameter.

    For Gaussian phase disorder, r = ⟨e^{iθ}⟩ = exp(−σ_θ²/2), so
    σ_θ = √(−2 ln r).  Replaces the ad-hoc (1−r)·π mapping.
    """
    r = float(np.clip(r, 1e-12, 1.0))
    return float(np.sqrt(max(0.0, -2.0 * np.log(r))))


def sample_phase_noise(N_grid, sigma_theta, xi_pixels, rng):
    """One correlated Gaussian phase-screen realization.

    White Gaussian noise smoothed to transverse correlation length
    xi_pixels, then rescaled so the *applied* RMS equals sigma_theta
    exactly (the old code's smoothing silently reduced the RMS ~10×).
    """
    from scipy.ndimage import gaussian_filter
    noise = rng.standard_normal((N_grid, N_grid))
    if xi_pixels > 0:
        noise = gaussian_filter(noise, sigma=xi_pixels)
    rms = noise.std()
    if rms > 0 and sigma_theta > 0:
        noise *= sigma_theta / rms
    else:
        noise = np.zeros_like(noise)
    return noise


# ---------------------------------------------------------------------------
# Cavity geometry
# ---------------------------------------------------------------------------
class CavityGeometry:
    """Microscopic diode cavity from the patent.

    The cathode-to-anode gap must be shorter than the particle mean free
    path for ballistic transport.  The patent specifies ~0.1 mm at
    atmospheric pressure (electron mfp ~9.3 mm).
    """

    def __init__(self, AK_gap=0.1e-3, width=1e-6, depth=0.1e-6,
                 n_cavities=10, n_channels=4, pressure_Pa=101325, T_gas=300,
                 cavity_Q=1000):
        self.AK_gap     = AK_gap       # cathode-to-anode distance [m]
        self.width      = width         # cavity width (perpendicular to E) [m]
        self.depth      = depth         # cavity depth (perpendicular to E and B) [m]
        self.n_cavities = n_cavities
        self.n_channels = n_channels
        self.pressure   = pressure_Pa
        self.T_gas      = T_gas
        self.cavity_Q   = cavity_Q     # cavity quality factor (max bounces)

    def number_density(self):
        """Ideal-gas number density n = p / (k_B T)  [m⁻³]."""
        return self.pressure / (k_B * self.T_gas)

    def mean_free_path(self, v_beam=None, k_L=K_LANGEVIN, cross_section=None):
        """Mean free path of a directed beam particle through stationary gas.

        λ = 1/(n σ) with n = p/(k_B T).  No √2: that factor applies to
        collisions among identical thermal particles, not to a beam
        traversing a stationary gas (fable5 E6/T8).

        Parameters
        ----------
        v_beam : float or None
            Beam speed [m/s].  If given, uses the Langevin capture
            cross-section σ_L = k_L / v_beam — the correct physics for a
            slow ion, where the induced-dipole attraction dominates.
        k_L : float
            Langevin rate constant [m³/s].
        cross_section : float or None
            Explicit cross-section [m²]; overrides the Langevin form.

        One of v_beam or cross_section must be supplied.
        """
        if cross_section is None:
            if v_beam is None:
                raise ValueError(
                    "mean_free_path needs v_beam (Langevin) or an explicit "
                    "cross_section")
            cross_section = k_L / v_beam
        return 1.0 / (self.number_density() * cross_section)

    def is_ballistic(self, v_beam=None, k_L=K_LANGEVIN, cross_section=None):
        mfp = self.mean_free_path(v_beam=v_beam, k_L=k_L,
                                  cross_section=cross_section)
        return self.AK_gap < mfp, mfp

    def effective_passes(self, v_beam=None, k_L=K_LANGEVIN,
                         cross_section=None):
        """Number of effective cavity passes before decoherence.

        The particle bounces in the resonant cavity, accumulating AB phase
        each pass.  The effective number of passes is limited by either:
          - cavity_Q (geometric/reflective quality factor), or
          - mfp / AK_gap (mean free path limits — gas collisions randomize phase)

        Returns (n_eff, limit_type) where limit_type is 'Q' or 'collision'.
        """
        mfp = self.mean_free_path(v_beam=v_beam, k_L=k_L,
                                  cross_section=cross_section)
        n_collision = mfp / self.AK_gap
        if self.cavity_Q <= n_collision:
            return self.cavity_Q, 'Q'
        else:
            return n_collision, 'collision'

    def info(self, v_beam=None):
        if v_beam is not None:
            ballistic, mfp = self.is_ballistic(v_beam=v_beam)
            n_eff, limit = self.effective_passes(v_beam=v_beam)
            print(f"  Cavity: AK gap = {self.AK_gap*1e3:.2f} mm, "
                  f"mfp = {mfp*1e3:.3g} mm  "
                  f"{'[ballistic]' if ballistic else '[WARNING: not ballistic]'}")
            print(f"  n_eff = {n_eff:.3g} passes ({limit}-limited)")
        else:
            print(f"  Cavity: AK gap = {self.AK_gap*1e3:.2f} mm "
                  f"(mfp requires beam velocity)")
        print(f"  {self.n_cavities} cavities × {self.n_channels} channels, "
              f"Q = {self.cavity_Q}")


# ---------------------------------------------------------------------------
# Main simulator
# ---------------------------------------------------------------------------
class CoherentMatterwaveBeam:
    """Simulates coherent matterwave beam generation via AB-effect
    Kuramoto synchronization in a diode cavity array.

    Parameters
    ----------
    species : str
        Key into SPECIES dict (e.g. 'electron', 'He+', 'Rb+').
    E_kinetic_eV : float
        Kinetic energy per particle [eV].
    B_field : float
        External magnetic field strength [T].
    cavity : CavityGeometry or None
        Cavity parameters.  Uses defaults if None.
    """

    def __init__(self, species='electron', E_kinetic_eV=1.0,
                 B_field=0.01, cavity=None,
                 dE_frac=0.01, beam_current_A=1e-6):

        sp = SPECIES[species]
        self.species_name = species
        self.symbol   = sp['symbol']
        self.label    = sp['label']
        self.mass     = sp['mass']
        self.charge_e = sp['charge']           # in units of e
        self.charge   = sp['charge'] * e_C     # in Coulombs

        self.E_k   = E_kinetic_eV * e_C        # kinetic energy [J]
        self.v      = np.sqrt(2 * self.E_k / self.mass)
        self.lam_dB = hbar * 2 * np.pi / (self.mass * self.v)
        self.k0     = 2 * np.pi / self.lam_dB

        self.B      = B_field
        self.cavity = cavity or CavityGeometry()

        # AB coupling constant β.
        # Patent §5.2 states β ≈ 1.02 × 10⁹ s⁻¹ for electrons at 1 T
        # with room-temperature speed v ~ 1e5 m/s.  The coupling scales as:
        #   β ∝ q × B × v
        # (charge couples to vector potential; field sets potential magnitude;
        # velocity sets traversal rate).  We calibrate from the patent's
        # stated value and scale for arbitrary species/conditions.
        beta_ref = 1.02e9   # s⁻¹, electrons at 1 T, v ~ 1e5 m/s
        B_ref    = 1.0      # T
        v_ref    = 1e5      # m/s (room-temp electron speed)
        q_ref    = e_C      # C

        self.beta = beta_ref * (self.charge / q_ref) * (self.B / B_ref) * (self.v / v_ref)
        self.tau_sync = 2.0 / self.beta if self.beta > 0 else np.inf

        # --- Species-dependent Kuramoto parameters ---
        self.dE_frac = dE_frac
        self.beam_current_A = beam_current_A

        # (a) Physical natural frequency spread [s⁻¹]
        # ΔE/E = dE_frac → Δv/v = dE_frac/2 (since E ∝ v²)
        # Spread in AB phase rate: Δω_AB = β × (Δv/v)
        self.omega_spread_phys = self.beta * self.dE_frac / 2

        # (b) Dimensionless frequency spread
        # σ_ω_dim = omega_spread_phys / β = dE_frac / 2
        self.omega_spread_dim_default = self.dE_frac / 2

        # (c) Effective particle count per synchronization volume
        # The cavity array has n_cavities × n_channels independent volumes.
        # Each synchronizes its own sub-ensemble; inter-cavity coherence comes
        # from the array geometry (patent §4), not from Kuramoto coupling.
        # N_per_volume = (I_beam / q) × transit_time / n_volumes
        self.n_volumes = self.cavity.n_cavities * self.cavity.n_channels
        self.transit_time = self.cavity.AK_gap / self.v
        N_total = (beam_current_A / self.charge) * self.transit_time
        self.N_per_volume = max(10, int(N_total / self.n_volumes))

        # (e) Dimensionless coupling from AB phase per transit
        # φ_single = β × τ_transit = AB phase accumulated per single pass [rad]
        # This is species-independent for same charge, B, cavity geometry
        # (velocity cancels: β ∝ v, τ_transit ∝ 1/v)
        self.phi_single = self.beta * self.transit_time

        # Multi-pass cavity resonance: particle bounces n_eff times,
        # accumulating AB phase each pass.  n_eff is limited by either
        # cavity Q or collision mean free path.  The mfp uses the species
        # Langevin cross-section σ_L = k_L/v (fable5 E6/T8).
        n_eff, K_limit = self.cavity.effective_passes(v_beam=self.v)
        self.n_eff = n_eff
        self.K_limit = K_limit
        self.sigma_langevin = K_LANGEVIN / self.v

        # K_phys = φ_single × n_eff = total accumulated AB phase [rad]
        # K_dim = K_phys / (2π) puts it in natural Kuramoto units
        self.K_phys = self.phi_single * self.n_eff
        self.K_dim = self.K_phys / (2 * np.pi)

        # (d) Single-transit scattering probability
        # In a ballistic single-pass device, the relevant dephasing is the
        # probability of one or more scattering events during a single transit.
        # p_scatter = 1 - exp(-AK_gap / mfp) ≈ AK_gap / mfp for AK_gap << mfp
        mfp = self.cavity.mean_free_path(v_beam=self.v)
        self.mfp_beam = mfp
        self.p_scatter = 1.0 - np.exp(-self.cavity.AK_gap / mfp)

    def space_charge_check(self, path_length, beam_radius,
                           current_A=None, verbose=True,
                           abort_on_fail=False):
        """Space-charge validation gate (fable5 T10/M1).

        From (I, v, beam radius, path length) compute the number of ions
        simultaneously in flight and the Coulomb/kinetic energy ratio.
        The single-particle propagation picture used by the whole
        pipeline requires occupancy ≤ 1 (one ion at a time, pattern
        accumulating statistically — T11); above that the gate fails
        loudly, in the same spirit as the spectrum/roundtrip gates.

        Returns a dict: ok, N_in_flight, coulomb_kinetic_ratio,
        I_max_single_A (largest current with occupancy ≤ 1), t_flight_s.
        """
        q = abs(self.charge)
        I = self.beam_current_A if current_A is None else current_A
        t_flight = path_length / self.v
        rate = I / q                       # ions per second
        N_flight = rate * t_flight         # simultaneous occupancy
        I_max_single = q / t_flight        # occupancy-1 current

        if N_flight > 1:
            # Mean inter-ion spacing in the beam cylinder and the mutual
            # Coulomb energy per ion relative to the kinetic energy.
            n = N_flight / (np.pi * beam_radius**2 * path_length)
            d = n ** (-1.0 / 3.0)
            E_coul = q**2 / (4 * np.pi * eps_0 * d)
            ratio = E_coul / self.E_k
        else:
            ratio = 0.0

        ok = N_flight <= 1
        result = {
            'ok': ok,
            'N_in_flight': float(N_flight),
            'coulomb_kinetic_ratio': float(ratio),
            'I_max_single_A': float(I_max_single),
            't_flight_s': float(t_flight),
            'current_A': float(I),
        }

        if verbose:
            status = 'PASS' if ok else 'FAIL'
            print(f"  SPACE-CHARGE GATE [{status}]: "
                  f"I = {I:.3e} A, t_flight = {t_flight*1e6:.3g} μs → "
                  f"{N_flight:.3g} ions in flight "
                  f"(single-ion limit I ≤ {I_max_single:.3e} A)")
            if not ok:
                print(f"    E_Coulomb/E_kinetic ≈ {ratio:.3g} — the "
                      f"single-particle picture is INVALID at this current; "
                      f"use one-ion-at-a-time statistical accumulation (T11).")
        if abort_on_fail and not ok:
            raise RuntimeError(
                f"Space-charge gate failed: {N_flight:.3g} ions in flight "
                f"(Coulomb/kinetic ≈ {ratio:.3g})")
        return result

    def transport_survival(self, L_path, pressure_Pa=None):
        """Fraction of the beam that traverses L_path without a collision.

        exp(−L_path/mfp) with the species Langevin mfp, evaluated at the
        transport-region pressure (defaults to the cavity pressure).
        Reported separately from cavity dephasing (fable5 T8): even a
        cavity that synchronizes may feed a transport leg that absorbs
        the beam.
        """
        if pressure_Pa is None:
            pressure_Pa = self.cavity.pressure
        n = pressure_Pa / (k_B * self.cavity.T_gas)
        if n <= 0:
            return 1.0
        mfp = 1.0 / (n * self.sigma_langevin)
        return float(np.exp(-L_path / mfp))

    def info(self):
        print("=" * 65)
        print("COHERENT MATTERWAVE BEAM SIMULATOR")
        print("  US Patent 9,502,202 B2 — Lockheed Martin (2016)")
        print("=" * 65)
        print(f"  Species: {self.label} ({self.symbol}), "
              f"q = {self.charge_e}e, m = {self.mass:.4e} kg")
        print(f"  E_k = {self.E_k/e_C:.2f} eV, "
              f"v = {self.v:.3e} m/s")
        print(f"  λ_dB = {self.lam_dB*1e12:.3f} pm")
        print(f"  B = {self.B:.4f} T ({self.B*1e4:.1f} gauss)")
        print(f"  β = {self.beta:.3e} s⁻¹  "
              f"(τ_sync ≈ {self.tau_sync*1e9:.2f} ns)")
        print(f"  ΔE/E = {self.dE_frac:.2%}, "
              f"σ_ω = {self.omega_spread_phys:.3e} s⁻¹ "
              f"(dim: {self.omega_spread_dim_default:.4f})")
        print(f"  I_beam = {self.beam_current_A*1e6:.1f} μA, "
              f"N_per_volume = {self.N_per_volume} "
              f"({self.n_volumes} volumes)")
        print(f"  Transit time = {self.transit_time*1e9:.3f} ns")
        print(f"  σ_Langevin = {self.sigma_langevin:.3e} m², "
              f"mfp = {self.mfp_beam:.3e} m")
        print(f"  p_scatter = {self.p_scatter:.4f} "
              f"(AK_gap/mfp = {self.cavity.AK_gap / self.mfp_beam:.4g})")
        print(f"  phi_single = {self.phi_single:.4f} rad/transit, "
              f"n_eff = {self.n_eff:.3g} ({self.K_limit}-limited)")
        print(f"  K_phys = {self.K_phys:.2f} rad (total), "
              f"K_dim = {self.K_dim:.3f}")
        self.cavity.info(v_beam=self.v)

    # -------------------------------------------------------------------
    # Kuramoto synchronization
    # -------------------------------------------------------------------
    def synchronize(self, N_particles=None, K=None, T_sync_dim=50,
                    alpha=0.5, K_scale=1.0, omega_spread_dim=None,
                    seed=None, verbose=True, mode='auto'):
        """Run Kuramoto synchronization of N_particles in the cavity.

        Both modes implement the same theory (fable5 E5/T9): first-order
        Kuramoto with a Lorentzian natural-frequency distribution of
        half-width γ = omega_spread_dim, for which K_c = 2γ and
        r = √(1 − K_c/K) are exact.  'ode' integrates the dynamics
        (dθ/dt = ω + K·coupling) in dimensionless time (dt=0.01), then
        maps to physical time using τ_sync; 'analytic' evaluates the
        closed form.

        For a perfectly monochromatic beam (patent assumption), set
        omega_spread_dim=0.  The patent claims zero critical coupling
        threshold for monochromatic beams.

        Parameters
        ----------
        N_particles : int or None
            Number of charged particles.  None uses species-dependent
            N_per_volume (from beam current and transit time).
        K : float or None
            Dimensionless Kuramoto coupling strength.  None uses
            physically derived K_dim (AB phase per transit / 2π).
        T_sync_dim : float
            Dimensionless integration time.
        alpha : float
            Unused; retained for API compatibility.  The inertial
            (second-order) model was removed in fable5 T9 — it is a
            different theory (bistability, discontinuous transition)
            than the analytic branch evaluates.
        K_scale : float
            Multiplier on coupling (K_eff = K × K_scale).
        omega_spread_dim : float or None
            Std dev of natural frequency distribution (dimensionless).
            None uses physically derived value from dE_frac.
            0 = perfectly monochromatic (patent assumption).
        seed : int or None
            RNG seed for reproducibility.
        verbose : bool
            Print progress.

        Returns
        -------
        result : dict
            theta, order_hist, times (physical [s]), r_final, beta_eff
        """
        if seed is not None:
            np.random.seed(seed)

        if N_particles is None:
            N_particles = self.N_per_volume
        if omega_spread_dim is None:
            omega_spread_dim = self.omega_spread_dim_default
        if K is None:
            K = self.K_dim

        K_eff = K * K_scale

        # Resolve mode
        if mode == 'auto':
            resolved_mode = 'analytic' if N_particles > 500 else 'ode'
        else:
            resolved_mode = mode

        if resolved_mode == 'analytic':
            return self._synchronize_analytic(
                N=N_particles, K_eff=K_eff, omega_spread_dim=omega_spread_dim,
                T_sync_dim=T_sync_dim, K_scale=K_scale,
                seed=seed, verbose=verbose)

        dt_k  = 0.01   # dimensionless timestep

        N = N_particles
        theta  = np.random.uniform(0, 2*np.pi, N)
        # Lorentzian (Cauchy) natural frequencies with half-width γ =
        # omega_spread_dim — the distribution for which the analytic
        # branch (K_c = 2γ, r = √(1−K_c/K)) is exact (fable5 T9).
        omega  = omega_spread_dim * np.random.standard_cauchy(N)

        # Single-pass dephasing: particles that scatter during transit
        # start with fully randomized phase and are effectively incoherent.
        # They act as a noise floor that the coupling must overcome.
        n_scattered = 0
        if self.p_scatter > 0:
            scattered = np.random.random(N) < self.p_scatter
            n_scattered = scattered.sum()
            if n_scattered > 0:
                theta[scattered] = np.random.uniform(0, 2*np.pi, n_scattered)
                omega[scattered] = (max(omega_spread_dim, 0.5)
                                    * np.random.standard_cauchy(n_scattered))

        n_steps = int(T_sync_dim / dt_k)
        order_hist = np.empty(n_steps)

        if verbose:
            print(f"\n  Kuramoto synchronization (1st order, Lorentzian): "
                  f"N={N}, K_eff={K_eff:.2f}, "
                  f"n_eff={self.n_eff:.3g}, γ_ω={omega_spread_dim:.4f}, "
                  f"p_scat={self.p_scatter:.4f}, T_dim={T_sync_dim}")

        for i in range(n_steps):
            z = np.mean(np.exp(1j * theta))
            order_hist[i] = np.abs(z)

            coupling = np.imag(z * np.exp(-1j * theta))
            theta   += (omega + K_eff * coupling) * dt_k

        # Steady-state estimate: average over the final 10% of the run
        # (instantaneous r fluctuates from the drifting Lorentzian tail).
        r_final = float(np.mean(order_hist[-max(1, n_steps // 10):]))

        # Map dimensionless time to physical time.
        # The Kuramoto converges in ~5-10 dimensionless time units.
        # Physically this corresponds to a few × τ_sync.
        time_scale = self.tau_sync / 5.0   # ~5 dim. units ≈ τ_sync
        times = np.arange(n_steps) * dt_k * time_scale

        if verbose:
            conv_dim = self._find_convergence_dim(order_hist, dt_k)
            conv_phys = conv_dim * time_scale
            print(f"    r_final = {r_final:.6f}  "
                  f"(converged at τ_dim ≈ {conv_dim:.1f}, "
                  f"t_phys ≈ {conv_phys*1e9:.2f} ns)")

        return {
            'theta':       theta % (2*np.pi),
            'order_hist':  order_hist,
            'times':       times,
            'r_final':     r_final,
            'beta_eff':    self.beta * K_scale,
            'N':           N,
            'T_max':       times[-1],
            'p_scatter':   self.p_scatter,
            'n_scattered': int(n_scattered),
            'mode':        'ode',
        }

    def _find_convergence_dim(self, order_hist, dt_k, threshold=0.95):
        """Dimensionless time at which order parameter first exceeds threshold."""
        idx = np.where(order_hist >= threshold)[0]
        if len(idx) == 0:
            return len(order_hist) * dt_k
        return idx[0] * dt_k

    def _synchronize_analytic(self, N, K_eff, omega_spread_dim,
                              T_sync_dim, K_scale, seed, verbose):
        """Compute steady-state Kuramoto order parameter analytically.

        Exact closed form for first-order Kuramoto with a Lorentzian
        frequency distribution of half-width γ = omega_spread_dim
        (fable5 T9): K_c = 2γ, r = √(1 − K_c/K) above threshold.  Adds
        finite-size and scattering corrections and generates synthetic
        time series / phase distributions for downstream compatibility.
        """
        # Critical coupling for Lorentzian(γ) frequency distribution
        K_c = 2 * omega_spread_dim

        # Steady-state order parameter
        if K_eff > K_c:
            r_inf = np.sqrt(1 - K_c / K_eff)
        else:
            r_inf = 0.0

        # Finite-size correction
        noise_floor = 1.0 / np.sqrt(N)
        r_inf = max(r_inf, noise_floor)
        r_inf = min(r_inf, 1.0 - noise_floor)

        # Scattering correction
        if self.p_scatter > 0:
            r_inf = r_inf * (1 - self.p_scatter) + (self.p_scatter / np.sqrt(N))

        # Synthetic time series: exponential approach
        if K_eff > K_c:
            tau_converge = 5.0 / (K_eff - K_c + 0.1)
        else:
            tau_converge = 50.0

        dt_k = 0.01
        n_steps = int(T_sync_dim / dt_k)
        t_dim = np.arange(n_steps) * dt_k
        order_hist = r_inf * (1 - np.exp(-t_dim / tau_converge))
        order_hist = np.maximum(order_hist, 1.0 / np.sqrt(N))

        # Synthetic final phases from von Mises distribution
        if r_inf >= 0.99:
            kappa = 100.0
        elif r_inf > 0.01:
            kappa = r_inf * (2 - r_inf**2) / (1 - r_inf**2)
        else:
            kappa = 0.0

        if kappa > 0:
            theta = np.random.vonmises(0, kappa, N)
            theta = theta % (2 * np.pi)
        else:
            theta = np.random.uniform(0, 2 * np.pi, N)

        # Physical time mapping (same as ODE path)
        time_scale = self.tau_sync / 5.0
        times = t_dim * time_scale

        if verbose:
            print(f"\n  Kuramoto (analytical): N={N}, K_eff={K_eff:.4f}, "
                  f"K_c={K_c:.4f}, \u03c3_\u03c9={omega_spread_dim:.4f}")
            print(f"    r_\u221e = {r_inf:.6f}  "
                  f"({'synchronized' if r_inf > 0.5 else 'unsynchronized'})")
            if self.p_scatter > 0.001:
                print(f"    p_scatter correction: {self.p_scatter:.4f}")

        return {
            'theta':       theta,
            'order_hist':  order_hist,
            'times':       times,
            'r_final':     float(r_inf),
            'beta_eff':    self.beta * K_scale,
            'N':           N,
            'T_max':       times[-1],
            'p_scatter':   self.p_scatter,
            'n_scattered': int(self.p_scatter * N),
            'mode':        'analytic',
        }

    def validate_kuramoto_modes(self, N=500, tol=0.15, T_sync_dim=200,
                                seed=42, verbose=True):
        """Validation gate (fable5 T9): analytic vs ODE agreement.

        Runs both modes at a size where the ODE is affordable and checks
        |r_ode − r_analytic| < tol.  Since both modes now implement the
        same theory (first-order Kuramoto, Lorentzian frequencies), any
        disagreement beyond finite-N/finite-T effects is a bug.

        Returns (ok, r_ode, r_analytic).
        """
        res_ode = self.synchronize(N_particles=N, mode='ode', seed=seed,
                                   T_sync_dim=T_sync_dim, verbose=False)
        res_ana = self.synchronize(N_particles=N, mode='analytic', seed=seed,
                                   T_sync_dim=T_sync_dim, verbose=False)
        gap = abs(res_ode['r_final'] - res_ana['r_final'])
        ok = gap < tol
        if verbose:
            print(f"  Kuramoto mode gate (N={N}): "
                  f"r_ode = {res_ode['r_final']:.4f}, "
                  f"r_analytic = {res_ana['r_final']:.4f}, "
                  f"|Δ| = {gap:.4f} {'< ' if ok else '>= '}{tol}  "
                  f"[{'OK' if ok else 'FAIL'}]")
        return ok, res_ode['r_final'], res_ana['r_final']

    @staticmethod
    def analytical_r(K_dim, omega_spread_dim, N=None, p_scatter=0.0):
        """Compute steady-state Kuramoto order parameter analytically.

        Exact for first-order Kuramoto with Lorentzian(γ=omega_spread_dim)
        frequencies: K_c = 2γ, r = √(1 − K_c/K).

        Parameters
        ----------
        K_dim : float
            Dimensionless coupling strength.
        omega_spread_dim : float
            Half-width γ of the Lorentzian frequency distribution
            (dimensionless).
        N : int or None
            Ensemble size for finite-size correction. None = infinite.
        p_scatter : float
            Single-transit scattering probability.

        Returns
        -------
        r : float
            Steady-state order parameter.
        """
        K_c = 2 * omega_spread_dim
        if K_dim > K_c:
            r = np.sqrt(1 - K_c / K_dim)
        else:
            r = 0.0

        if N is not None and N > 0:
            noise_floor = 1.0 / np.sqrt(N)
            r = max(r, noise_floor)
            r = min(r, 1.0 - noise_floor)

        if p_scatter > 0:
            N_eff = N if N is not None else int(1e12)
            r = r * (1 - p_scatter) + p_scatter / np.sqrt(N_eff)

        return float(np.clip(r, 0.0, 1.0))

    # -------------------------------------------------------------------
    # AB phase imprinting (single-cavity demonstration)
    # -------------------------------------------------------------------
    def ab_phase_screen(self, N_grid=256, L=1e-6):
        """Generate a sample AB phase screen from the cavity B-field.

        Returns a 2D phase map that a particle beam traversing the cavity
        would acquire from the vector potential.  This connects the beam
        generator to the downstream holographic pipeline.
        """
        x = np.linspace(-L/2, L/2, N_grid)
        y = np.linspace(-L/2, L/2, N_grid)
        X, Y = np.meshgrid(x, y, indexing='ij')

        # For uniform B along z, A = (B/2)(-y, x, 0) in symmetric gauge
        # Phase = (q/ℏ) ∫ A·ds along beam direction (x)
        # For beam along x: ∫ A_x dx = -(B/2) y · L_traverse
        L_traverse = self.cavity.AK_gap
        phase = (self.charge / hbar) * (self.B / 2) * (-Y) * L_traverse

        return phase, x, y

    # -------------------------------------------------------------------
    # Multi-species comparison
    # -------------------------------------------------------------------
    @staticmethod
    def compare_species(species_list=None, E_kinetic_eV=1.0,
                        B_field=0.01, N_particles=None, seed=42,
                        dE_frac=0.01, beam_current_A=1e-6,
                        verbose=True):
        """Run synchronization for multiple species and compare.

        Returns list of (simulator, result) tuples.
        """
        if species_list is None:
            species_list = ['electron', 'He+', 'Na+', 'Rb+']

        runs = []
        for sp in species_list:
            sim = CoherentMatterwaveBeam(species=sp, E_kinetic_eV=E_kinetic_eV,
                                         B_field=B_field,
                                         dE_frac=dE_frac,
                                         beam_current_A=beam_current_A)
            if verbose:
                sim.info()
            result = sim.synchronize(N_particles=N_particles, seed=seed,
                                     verbose=verbose)
            runs.append((sim, result))
        return runs

    # -------------------------------------------------------------------
    # Beam wavefunction construction
    # -------------------------------------------------------------------
    def build_beam(self, sync_result, N_grid=256, L=None, sigma_frac=0.35,
                   xi_perp=None, rng=None):
        """Construct one 2D beam-wavefunction realization from sync result.

        Gaussian transverse envelope with one realization of correlated
        phase noise of RMS σ_θ = √(−2 ln r) and transverse correlation
        length ξ⊥ (fable5 E2/T6).  Partial coherence means the deposited
        intensity is the *ensemble average* over such realizations —
        callers modelling r < 1 should average |ψ|² over many draws.

        Parameters
        ----------
        sync_result : dict
            Output of self.synchronize().
        N_grid : int
            Spatial grid resolution.
        L : float or None
            Physical domain size [m].  Defaults to 1000 × λ_dB.
        sigma_frac : float
            Beam waist as fraction of L.
        xi_perp : float or None
            Transverse correlation length of the phase noise [m].
            Defaults to L/8.
        rng : np.random.Generator or None
            Source of randomness (fresh default_rng() if None).

        Returns
        -------
        psi : ndarray, complex128, shape (N_grid, N_grid)
            Normalized wavefunction.
        x : ndarray
            Coordinate array [m].
        """
        if L is None:
            L = 1000 * self.lam_dB
        if xi_perp is None:
            xi_perp = L / 8
        if rng is None:
            rng = np.random.default_rng()

        dx = L / N_grid
        x  = np.linspace(-L/2, L/2, N_grid)
        X, Y = np.meshgrid(x, x, indexing='ij')

        sigma = sigma_frac * L
        r = sync_result['r_final']

        # Gaussian transverse envelope.  The forward momentum hbar*k0 is
        # carried by the propagation axis (kz in the angular-spectrum
        # propagator), not a transverse phase ramp (fable5 E1).
        psi = np.exp(-(X**2 + Y**2) / (2*sigma**2)).astype(np.complex128)

        # One phase-noise realization at the physical RMS
        sigma_theta = coherence_sigma_theta(r)
        noise = sample_phase_noise(N_grid, sigma_theta, xi_perp / dx, rng)
        psi = psi * np.exp(1j * noise)

        # Normalize
        psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * dx**2)

        return psi, x

    # -------------------------------------------------------------------
    # Visualization
    # -------------------------------------------------------------------
    def plot_synchronization(self, result, save_path=None):
        """Plot Kuramoto synchronization dynamics."""
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(figsize=(14, 8))
        gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

        t_ns = result['times'] * 1e9

        # --- Order parameter vs time ---
        ax0 = fig.add_subplot(gs[0, 0:2])
        ax0.plot(t_ns, result['order_hist'], 'b-', linewidth=1.5)
        ax0.axhline(0.95, color='r', ls='--', alpha=0.5, label='r = 0.95')
        ax0.axvline(self.tau_sync*1e9, color='g', ls=':', alpha=0.5,
                    label=f'τ_sync = {self.tau_sync*1e9:.2f} ns')
        ax0.set_xlabel('Time [ns]')
        ax0.set_ylabel('Order parameter r')
        ax0.set_title(f'{self.label} ({self.symbol}) — Kuramoto Synchronization')
        ax0.set_ylim(-0.05, 1.05)
        ax0.legend(loc='lower right')
        ax0.grid(True, alpha=0.3)

        # --- Final phase distribution ---
        ax1 = fig.add_subplot(gs[0, 2], projection='polar')
        theta = result['theta']
        ax1.hist(theta, bins=36, density=True, alpha=0.7, color='steelblue',
                 edgecolor='navy', linewidth=0.5)
        ax1.set_title(f'Final phases\nr = {result["r_final"]:.4f}', pad=15)

        # --- Phase evolution snapshots (unit circle) ---
        ax2 = fig.add_subplot(gs[1, 0])
        n_steps = len(result['order_hist'])
        # Re-run in dimensionless time to capture snapshots
        np.random.seed(42)
        N = result['N']
        th = np.random.uniform(0, 2*np.pi, N)
        om = self.omega_spread_dim_default * np.random.standard_cauchy(N)
        if self.p_scatter > 0:
            scattered = np.random.random(N) < self.p_scatter
            if scattered.sum() > 0:
                th[scattered] = np.random.uniform(0, 2*np.pi, scattered.sum())
                om[scattered] = (max(self.omega_spread_dim_default, 0.5)
                                 * np.random.standard_cauchy(scattered.sum()))
        dt_k = 0.01
        K_snap = self.K_dim
        snap_indices = [0, n_steps//10, n_steps//4, n_steps//2]
        snap_thetas = []
        for i in range(n_steps):
            if i in snap_indices:
                snap_thetas.append(th.copy())
            z = np.mean(np.exp(1j * th))
            coupling = np.imag(z * np.exp(-1j * th))
            th += (om + K_snap * coupling) * dt_k

        colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']
        labels = ['t=0', 't=T/10', 't=T/4', 't=T/2']
        for th_snap, col, lab in zip(snap_thetas, colors, labels):
            x_c = np.cos(th_snap)
            y_c = np.sin(th_snap)
            ax2.scatter(x_c, y_c, s=4, alpha=0.5, color=col, label=lab)
        circle = plt.Circle((0, 0), 1, fill=False, color='gray', ls='--', alpha=0.4)
        ax2.add_patch(circle)
        ax2.set_xlim(-1.3, 1.3)
        ax2.set_ylim(-1.3, 1.3)
        ax2.set_aspect('equal')
        ax2.set_title('Phase snapshots (unit circle)')
        ax2.legend(fontsize=7, loc='upper right')
        ax2.grid(True, alpha=0.3)

        # --- Species parameters ---
        ax3 = fig.add_subplot(gs[1, 1:])
        ax3.axis('off')
        params = [
            ['Species',      f'{self.label} ({self.symbol})'],
            ['Mass',         f'{self.mass:.4e} kg'],
            ['Charge',       f'{self.charge_e}e'],
            ['E_kinetic',    f'{self.E_k/e_C:.2f} eV'],
            ['Velocity',     f'{self.v:.3e} m/s'],
            ['λ_dB',         f'{self.lam_dB*1e12:.3f} pm'],
            ['B field',      f'{self.B*1e4:.1f} gauss ({self.B:.4f} T)'],
            ['β',            f'{self.beta:.3e} s⁻¹'],
            ['τ_sync',       f'{self.tau_sync*1e9:.2f} ns'],
            ['ΔE/E',         f'{self.dE_frac:.2%}'],
            ['σ_ω (dim)',    f'{self.omega_spread_dim_default:.4f}'],
            ['φ/transit',    f'{self.phi_single:.4f} rad'],
            ['n_eff',        f'{self.n_eff:.1f} ({self.K_limit})'],
            ['K_dim',        f'{self.K_dim:.3f}'],
            ['N_per_vol',    f'{self.N_per_volume}'],
            ['p_scatter',    f'{self.p_scatter:.4f}'],
            ['r_final',      f'{result["r_final"]:.6f}'],
            ['N_particles',  f'{result["N"]}'],
            ['Ballistic?',   f'{"Yes" if self.cavity.is_ballistic()[0] else "NO"}'],
        ]
        table = ax3.table(cellText=params, colLabels=['Parameter', 'Value'],
                          loc='center', cellLoc='left')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.0, 1.4)
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor('#4472C4')
                cell.set_text_props(color='white', weight='bold')

        fig.suptitle('Coherent Matterwave Beam — AB-Effect Kuramoto Synchronization',
                     fontsize=13, weight='bold', y=0.98)

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  Saved: {save_path}")
        plt.close(fig)
        return fig

    def plot_species_comparison(self, runs, save_path=None):
        """Plot synchronization comparison across species."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # --- Order parameter vs time ---
        ax = axes[0]
        for sim, result in runs:
            t_ns = result['times'] * 1e9
            ax.plot(t_ns, result['order_hist'], linewidth=1.5,
                    label=f'{sim.symbol} (τ={sim.tau_sync*1e9:.1f} ns)')
        ax.axhline(0.95, color='k', ls='--', alpha=0.3)
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Order parameter r')
        ax.set_title('Synchronization Dynamics')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # --- Summary table ---
        ax2 = axes[1]
        ax2.axis('off')
        rows = []
        for sim, result in runs:
            rows.append([
                sim.symbol,
                f'{sim.mass:.3e}',
                f'{sim.v:.3e}',
                f'{sim.lam_dB*1e12:.2f}',
                f'{sim.beta:.3e}',
                f'{sim.tau_sync*1e9:.2f}',
                f'{sim.omega_spread_dim_default:.4f}',
                f'{sim.phi_single:.4f}',
                f'{sim.n_eff:.0f}',
                f'{sim.K_dim:.3f}',
                f'{sim.N_per_volume}',
                f'{sim.p_scatter:.4f}',
                f'{result["r_final"]:.4f}',
            ])
        col_labels = ['Species', 'Mass [kg]', 'v [m/s]', 'λ_dB [pm]',
                      'β [s⁻¹]', 'τ_sync [ns]', 'σ_ω (dim)',
                      'φ/pass', 'n_eff', 'K_dim', 'N_vol', 'p_scat', 'r_final']
        table = ax2.table(cellText=rows, colLabels=col_labels,
                          loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1.0, 1.5)
        for (row, col), cell in table.get_celld().items():
            if row == 0:
                cell.set_facecolor('#4472C4')
                cell.set_text_props(color='white', weight='bold')
        ax2.set_title('Species Comparison', pad=20)

        fig.suptitle('Multi-Species Coherent Matterwave Beam',
                     fontsize=13, weight='bold')
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"  Saved: {save_path}")
        plt.close(fig)
        return fig


# ===========================================================================
# Standalone execution
# ===========================================================================
if __name__ == '__main__':
    import os
    import matplotlib.pyplot as plt

    os.makedirs('results', exist_ok=True)

    # --- Single-species demo: electron beam ---
    print("\n" + "="*65)
    print("DEMO 1: Electron beam synchronization")
    print("="*65)
    sim_e = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0,
                                   B_field=0.01)
    sim_e.info()
    res_e = sim_e.synchronize(seed=42)
    sim_e.plot_synchronization(res_e, save_path='results/cmb_electron_sync.png')

    # --- Single-species demo: He+ beam ---
    print("\n" + "="*65)
    print("DEMO 2: He+ ion beam synchronization")
    print("="*65)
    sim_he = CoherentMatterwaveBeam(species='He+', E_kinetic_eV=1.0,
                                    B_field=0.01)
    sim_he.info()
    res_he = sim_he.synchronize(seed=42)
    sim_he.plot_synchronization(res_he, save_path='results/cmb_he_ion_sync.png')

    # --- Multi-species comparison ---
    print("\n" + "="*65)
    print("DEMO 3: Multi-species comparison")
    print("="*65)
    runs = CoherentMatterwaveBeam.compare_species(
        species_list=['electron', 'He+', 'Na+', 'Rb+'],
        E_kinetic_eV=1.0, B_field=0.01, seed=42)
    # Use the first simulator for the plot method
    runs[0][0].plot_species_comparison(runs,
        save_path='results/cmb_species_comparison.png')

    # --- Energy spread comparison ---
    print("\n" + "="*65)
    print("DEMO 4: Effect of energy spread on synchronization")
    print("="*65)
    dE_fracs = [0.001, 0.01, 0.05, 0.20]
    fig, ax = plt.subplots(figsize=(8, 5))
    for dE_frac in dE_fracs:
        sim = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0,
                                     B_field=0.01, dE_frac=dE_frac)
        res = sim.synchronize(seed=42, verbose=False)
        t_ns = res['times'] * 1e9
        ax.plot(t_ns, res['order_hist'], linewidth=1.5,
                label=f'ΔE/E = {dE_frac:.1%}')
    ax.axhline(0.95, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Time [ns]')
    ax.set_ylabel('Order parameter r')
    ax.set_title('Electron Beam: Energy Spread vs Synchronization')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('results/cmb_energy_spread.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: results/cmb_energy_spread.png")

    # --- Beam wavefunction output ---
    print("\n" + "="*65)
    print("DEMO 5: Beam wavefunction (He+ for holographic pipeline)")
    print("="*65)
    sim_he2 = CoherentMatterwaveBeam(species='He+', E_kinetic_eV=1.0,
                                      B_field=0.01)
    res_he2 = sim_he2.synchronize(seed=42)
    psi, x = sim_he2.build_beam(res_he2, N_grid=256)
    density = np.abs(psi)**2
    phase   = np.angle(psi)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    extent = [x[0]*1e9, x[-1]*1e9, x[0]*1e9, x[-1]*1e9]
    im0 = axes[0].imshow(density.T, origin='lower', extent=extent, cmap='inferno')
    axes[0].set_title(f'He⁺ beam density (r={res_he2["r_final"]:.4f})')
    axes[0].set_xlabel('x [nm]')
    axes[0].set_ylabel('y [nm]')
    plt.colorbar(im0, ax=axes[0], label='|ψ|²')
    im1 = axes[1].imshow(phase.T, origin='lower', extent=extent, cmap='twilight')
    axes[1].set_title('Phase')
    axes[1].set_xlabel('x [nm]')
    axes[1].set_ylabel('y [nm]')
    plt.colorbar(im1, ax=axes[1], label='arg(ψ) [rad]')
    fig.suptitle('Coherent He⁺ Beam — Ready for Holographic Pipeline',
                 fontsize=12, weight='bold')
    fig.tight_layout()
    fig.savefig('results/cmb_he_beam_wavefunction.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: results/cmb_he_beam_wavefunction.png")

    # --- B-field comparison (vacuum cavity to see multi-pass effect) ---
    print("\n" + "="*65)
    print("DEMO 6: Effect of B-field on coupling and synchronization")
    print("="*65)
    B_fields = [0.001, 0.005, 0.01, 0.05, 0.1]
    cav_vac = CavityGeometry(pressure_Pa=100)
    fig, ax = plt.subplots(figsize=(8, 5))
    for B in B_fields:
        sim = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0,
                                     B_field=B, cavity=cav_vac)
        res = sim.synchronize(seed=42, verbose=False)
        t_ns = res['times'] * 1e9
        ax.plot(t_ns, res['order_hist'], linewidth=1.5,
                label=f'B={B*1e4:.0f} G, K={sim.K_dim:.2f}, n={sim.n_eff:.0f}')
    ax.axhline(0.95, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Time [ns]')
    ax.set_ylabel('Order parameter r')
    ax.set_title('Electron Beam: B-field → Coupling Strength (vacuum cavity)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('results/cmb_B_field_sweep.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: results/cmb_B_field_sweep.png")

    # --- Pressure sweep ---
    print("\n" + "="*65)
    print("DEMO 7: Pressure sweep — collision-limited vs Q-limited")
    print("="*65)
    pressures = [1e-3, 1.0, 100, 10000, 101325]
    fig, ax = plt.subplots(figsize=(8, 5))
    for P in pressures:
        cav = CavityGeometry(pressure_Pa=P)
        sim = CoherentMatterwaveBeam(species='electron', E_kinetic_eV=1.0,
                                     B_field=0.01, cavity=cav)
        res = sim.synchronize(seed=42, verbose=False)
        t_ns = res['times'] * 1e9
        label = f'P={P:.0e} Pa, n={sim.n_eff:.0f} ({sim.K_limit}), K={sim.K_dim:.2f}'
        ax.plot(t_ns, res['order_hist'], linewidth=1.5, label=label)
    ax.axhline(0.95, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Time [ns]')
    ax.set_ylabel('Order parameter r')
    ax.set_title('Electron Beam: Pressure → Effective Passes → Synchronization')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig('results/cmb_pressure_sweep.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved: results/cmb_pressure_sweep.png")

    print("\n  All demos complete.")
