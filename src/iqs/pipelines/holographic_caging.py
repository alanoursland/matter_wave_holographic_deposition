"""
Integrated Quantum Substrate v10
=================================

End-to-end pipeline for controllable matter deposition:

  He⁺ source → ideal phase-plate holography → deposition model

Physical concept:
  1. A coherent He⁺ ion beam is generated in a diode cavity array via
     Aharonov-Bohm phase synchronization (US Patent 9,502,202 B2).
     The Kuramoto order parameter r quantifies beam coherence and
     depends on B-field, cavity pressure, and energy spread.
  2. The charged He⁺ beam passes through an ideal programmable phase
     plate. The legacy SQUID control grid parameterizes that ideal plate;
     a physical coplanar-loop realization fails the axial phase-authority
     gate and remains an implementation hypothesis.
  3. Loop phases are computed by inverse holography to produce a desired
     target deposition pattern after Fresnel propagation.
  4. The deposited pattern is stabilized by Aharonov-Bohm caging on a
     2D diamond lattice (flat bands at Φ=π).

The corrected pipeline treats transport bandwidth, transverse coherence,
and longitudinal wavelength spread as ensemble properties. The v10 stage is
not wavelength-matched; see reports/v11_report.md.

Dependencies: coherent_matterwave_beam.py, inverse_holography.py,
              sim_v9.py, diamond_caging.py
"""

import os
import sys
import time
import numpy as np
import torch

from iqs.sources.coherent_matterwave import (
    CoherentMatterwaveBeam, CavityGeometry, SPECIES,
    sample_phase_noise,
    hbar, k_B, e_C, m_u,
)
from iqs.holography import (
    SQUIDArray, InverseHolographySolver,
    target_single_spot, target_line, target_grid_of_dots,
    target_ring, target_letter, smooth_target, compute_metrics,
    validate_roundtrip,
)
from iqs.pipelines.patterned_substrate import IntegratedQuantumSubstrate
from iqs.constants import m_He
from iqs.numerics.metrics import michelson_contrast, ssim_score, min_feature_size
from iqs.lattices.diamond import DiamondNetwork
from iqs.lattices.density_mapping import DensityPeakMapper
from iqs.numerics.device import get_device
from iqs.sources import (
    SourceParams, DirectSource, KuramotoPatentSource,
    gaussian_wavelength_samples, fractional_rms,
)

_device = get_device()
# ===========================================================================
# INTEGRATED v10 PIPELINE
# ===========================================================================

class IntegratedPipelineV10:
    """
    End-to-end: He⁺ CMWB source → SQUID holography
    → (optional Floquet) → diamond caging.

    The CMWB generates a coherent He⁺ beam whose order parameter r
    depends on cavity pressure, B-field, and energy spread.  The
    charged beam passes through the SQUID array, where each loop
    imprints an AB phase shift φ = (q/ℏ)∮A·dl.  Charge is required
    for the AB effect.

    Parameters
    ----------
    B_field : float
        Magnetic field for AB synchronization in the CMWB cavity [T].
    pressure_Pa : float or None
        CMWB cavity pressure [Pa].  Determines effective passes.
        None (default): use the supplied cavity's pressure, or 100 Pa
        for a freshly created cavity.  Explicit values override the
        supplied cavity's pressure.
    dE_frac : float
        Beam energy spread ΔE/E.
    beam_current_A : float
        CMWB beam current [A].
    cavity : CavityGeometry or None
        Cavity parameters.  Pressure is overridden by pressure_Pa.
    T_beam : float
        Beam temperature [K].  Sets λ_dB for the He⁺ beam.
    N : int
        Spatial grid resolution.
    L : float
        Physical domain size [m].
    N_loops : int
        SQUID array resolution per axis.
    prop_distance_lam : float
        Propagation distance in units of de Broglie wavelength.
    use_floquet : bool
        Whether to apply Floquet dressing after holography.
    V_max_floquet : float
        Floquet drive strength (if enabled).
    n_resonant : int
        Floquet sideband for binding filter (if enabled).
    N_diamond_x, N_diamond_y : int
        Diamond lattice dimensions for caging.
    """

    def __init__(self, B_field=0.01, pressure_Pa=None,
                 dE_frac=0.01, beam_current_A=1e-6,
                 cavity=None, T_beam=1e-3,
                 N=256, L=400e-9, N_loops=32,
                 prop_distance_lam=20.0,
                 use_floquet=False, V_max_floquet=0.9,
                 n_resonant=0,
                 N_diamond_x=16, N_diamond_y=16,
                 coherence_xi=50e-9, n_noise_realizations=100,
                 n_wavelength_samples=5,
                 phase_actuator='achromatic',
                 electrostatic_plate=None,
                 source=None):
        self.N = N
        self.L = L
        self.dx = L / N
        self.N_loops = N_loops
        # Partial-coherence model (fable5 E2/T6): transverse correlation
        # length ξ⊥ of the source phase noise, and ensemble size M for
        # the intensity average that models incoherent ion arrivals.
        self.coherence_xi = coherence_xi
        self.n_noise_realizations = n_noise_realizations
        self.n_wavelength_samples = n_wavelength_samples
        self.phase_actuator = phase_actuator
        self.electrostatic_plate = electrostatic_plate
        self.prop_distance_lam = prop_distance_lam
        self.use_floquet = use_floquet
        self.V_max_floquet = V_max_floquet
        self.n_resonant = n_resonant
        self.Lx = N_diamond_x
        self.Ly = N_diamond_y
        self.T_beam = T_beam

        # --- CMWB source (He⁺) ---
        # pressure_Pa=None means: use the supplied cavity's own pressure
        # (or the 100 Pa default for a fresh cavity).  An explicit
        # pressure_Pa always wins (fable5 E7: the old logic silently
        # clobbered a user-supplied cavity's pressure with the default).
        if cavity is None:
            cav = CavityGeometry(
                pressure_Pa=100 if pressure_Pa is None else pressure_Pa)
        else:
            cav = cavity
            if pressure_Pa is not None:
                cav.pressure = pressure_Pa
        self.cmwb = CoherentMatterwaveBeam(
            species='He+', E_kinetic_eV=k_B * T_beam / e_C,
            B_field=B_field, cavity=cav,
            dE_frac=dE_frac, beam_current_A=beam_current_A,
        )

        # --- Source-parameter provider (fable5 T12/M2) ---
        # Downstream stages consume only SourceParams (Δλ/λ, ξ⊥, current,
        # σ_θ).  Default provider: the patent's AB-Kuramoto model wrapping
        # the CMWB above.  Pass source=DirectSource(...) to bypass the
        # Kuramoto abstraction entirely.
        if source is None:
            source = KuramotoPatentSource(self.cmwb, xi_perp=coherence_xi)
        self.source = source
        self._source_params = None
        self._wavelength_samples = None

        # --- He⁺ beam parameters ---
        self.mass = m_He
        self.v = np.sqrt(2 * k_B * T_beam / self.mass)
        self.k0 = self.mass * self.v / hbar
        self.lam = 2 * np.pi / self.k0

        # --- SQUID holography solver ---
        self.squid = SQUIDArray(N_loops=N_loops, N_grid=N, L_grid=L)
        self.holo_solver = InverseHolographySolver(
            squid_array=self.squid, N=N, L=L, T_beam=T_beam,
            prop_distance_lam=prop_distance_lam,
            phase_response=phase_actuator,
        )
        # Use one exact central wavelength for reporting, conditioning, and
        # every spectral component in the propagation ensemble.
        self.mass = self.holo_solver.mass
        self.v = self.holo_solver.v
        self.k0 = self.holo_solver.k0
        self.lam = self.holo_solver.lam
        if (electrostatic_plate is not None
                and self.holo_solver.phase_response.name != 'electrostatic'):
            raise ValueError(
                "electrostatic_plate requires phase_actuator='electrostatic'")
        if (electrostatic_plate is not None
                and not callable(getattr(electrostatic_plate, 'evaluate', None))):
            raise TypeError("electrostatic_plate must provide evaluate()")

        # --- v9 substrate sim (for optional Floquet) ---
        if self.use_floquet:
            self.substrate = IntegratedQuantumSubstrate(
                N=N, L=L, T_beam=T_beam,
            )

        # --- Diamond caging ---
        self.diamond = DiamondNetwork(Lx=N_diamond_x, Ly=N_diamond_y)
        self.mapper  = DensityPeakMapper(N_grid=N)

    def info(self):
        print("=" * 65)
        print("INTEGRATED QUANTUM SUBSTRATE  v10")
        print("  He⁺ CMWB → SQUID holography → caging")
        print("=" * 65)
        print(f"  CMWB source: He⁺ at {self.T_beam*1e3:.1f} mK")
        print(f"    B = {self.cmwb.B*1e4:.0f} G, "
              f"P = {self.cmwb.cavity.pressure:.0f} Pa, "
              f"ΔE/E = {self.cmwb.dE_frac:.1%}")
        print(f"    K_dim = {self.cmwb.K_dim:.3f}, "
              f"n_eff = {self.cmwb.n_eff:.0f} "
              f"({self.cmwb.K_limit}-limited), "
              f"N_vol = {self.cmwb.N_per_volume}")
        print(f"  He⁺ beam:")
        print(f"    λ_dB = {self.lam*1e9:.2f} nm, "
              f"v = {self.v:.3f} m/s")
        print(f"  Grid: {self.N}×{self.N}, L = {self.L*1e9:.0f} nm, "
              f"dx = {self.dx*1e9:.2f} nm")
        print(f"  SQUID: {self.N_loops}×{self.N_loops} "
              f"(pitch = {self.squid.pitch*1e9:.1f} nm)")
        print(f"  Phase actuator: {self.holo_solver.phase_response.name}")
        print(f"  Propagation: {self.prop_distance_lam}λ "
              f"= {self.prop_distance_lam * self.lam * 1e9:.1f} nm")
        print(f"  Floquet: {'ON' if self.use_floquet else 'OFF'}")
        print(f"  Diamond: {self.Lx}×{self.Ly}")
        print(f"  Device: {_device}")

    # -------------------------------------------------------------------
    # Stage 1: CMWB He⁺ source
    # -------------------------------------------------------------------
    def generate_source(self, seed=42, verbose=True):
        """Generate one He⁺ beam-wavefunction realization.

        The source provider (fable5 T12) supplies the physical parameters
        (Δλ/λ, ξ⊥, current, σ_θ); this stage consumes only those.

        Returns
        -------
        psi : ndarray (N, N) complex128, normalized
        sync_result : dict (Kuramoto sync result, or a direct-provider
            stand-in with 'r_final'/'mode')
        """
        if verbose:
            print("\n  STAGE 1: He⁺ source "
                  f"[provider: {type(self.source).__name__}]")

        params, sync = self.source.source_params(seed=seed, verbose=verbose)
        self._source_params = params
        self._wavelength_samples = gaussian_wavelength_samples(
            self.lam, params.dlam_frac, self.n_wavelength_samples
        )
        r = params.r_equivalent

        # Build the He⁺ beam wavefunction.
        x = np.linspace(-self.L / 2, self.L / 2, self.N)
        X, Y = np.meshgrid(x, x, indexing='ij')
        sigma = 0.35 * self.L

        # Transverse envelope only: the forward momentum hbar*k0 lives in the
        # propagator's kz, so a normally incident beam has no transverse
        # phase ramp (fable5 E1: exp(i*k0*X) here put the whole spectrum at
        # the evanescent cutoff).
        psi = np.exp(-(X**2 + Y**2) / (2 * sigma**2)).astype(np.complex128)

        # One phase-noise realization at the physical RMS σ_θ, correlated
        # over ξ⊥ and renormalized so the applied RMS equals the nominal
        # (fable5 E2/T6; the old (1−r)π + smoothing applied ~10× less
        # noise than stated).
        rng = np.random.default_rng(seed)
        noise = sample_phase_noise(self.N, params.sigma_theta,
                                   params.xi_perp / self.dx, rng)
        psi = psi * np.exp(1j * noise)

        # Normalize
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * self.dx**2)

        if verbose:
            print(f"    He⁺ beam: r_eq = {r:.4f}, "
                  f"λ = {self.lam*1e9:.2f} nm  "
                  f"[sync mode: {sync.get('mode', 'unknown')}]")
            print(f"    Phase noise: σ_θ nominal = {params.sigma_theta:.4f} rad, "
                  f"applied RMS = {noise.std():.4f} rad, "
                  f"ξ⊥ = {params.xi_perp*1e9:.0f} nm, "
                  f"Δλ/λ = {params.dlam_frac:.2%}")
            print(f"    Wavelength ensemble: Q = "
                  f"{len(self._wavelength_samples)}, "
                  f"applied RMS Δλ/λ = "
                  f"{fractional_rms(self._wavelength_samples):.2%}")

        return psi, sync

    # -------------------------------------------------------------------
    # Stage 2: Inverse holography
    # -------------------------------------------------------------------
    def solve_holography(self, target_pattern, method='gd',
                         n_iter_gs=300, n_restarts_gs=3,
                         n_iter_gd=500, lr_gd=0.05,
                         verbose=True):
        """Solve the inverse holography problem for a target pattern."""
        if verbose:
            print(f"\n  STAGE 2: Inverse holography ({method.upper()})")

        # Condition the target to the *deliverable* bandwidth: array
        # Nyquist AND the geometric transport aperture (fable5 T14) —
        # the optimizer must not chase evanescent/out-of-frame content.
        wavelengths = self._wavelength_samples
        if wavelengths is None:
            wavelengths = np.array([self.lam])
        self.holo_solver.set_wavelength_ensemble(wavelengths)
        if method == 'gs' and self.holo_solver.is_polychromatic:
            raise ValueError(
                "Gerchberg-Saxton cannot optimize an incoherent wavelength "
                "mixture; use method='gd' or set dlam_frac=0"
            )
        # Condition to the longest sampled wavelength, which has the
        # narrowest propagating band. Every source component can therefore
        # reach the target bandwidth presented to the optimizer.
        k0_condition = 2 * np.pi / float(np.max(wavelengths))
        target_smooth = smooth_target(
            target_pattern, N_loops=self.N_loops,
            corner_radius=0.03, sigma=2,
            k0=k0_condition, L=self.L, z=self.holo_solver.z,
        )

        t0 = time.time()
        if method == 'gs':
            raw = self.holo_solver.solve_gerchberg_saxton(
                target_smooth, n_iter=n_iter_gs,
                n_restarts=n_restarts_gs, verbose=verbose)
        elif method == 'gd':
            raw = self.holo_solver.solve_gradient_descent(
                target_smooth, n_iter=n_iter_gd,
                lr=lr_gd, verbose=verbose)
        else:
            raise ValueError(f"Unknown method: {method}")
        elapsed = time.time() - t0

        metrics = compute_metrics(raw['achieved'], target_smooth)
        screen_t = torch.tensor(
            raw['phase_screen'], dtype=torch.float64, device=_device)
        with torch.no_grad():
            central = self.holo_solver.forward_central(screen_t)
        central_np = central.detach().cpu().numpy()
        metrics_central = compute_metrics(central_np, target_smooth)
        electrostatic_report = None
        if (self.holo_solver.phase_response.name == 'electrostatic'
                and self.electrostatic_plate is not None):
            controls = np.asarray(raw['phi_loops_control']).reshape(
                self.N_loops, self.N_loops)
            electrostatic_report = self.electrostatic_plate.evaluate(controls)

        if verbose:
            label = "polychromatic" if self.holo_solver.is_polychromatic \
                else "central"
            print(f"    SSIM = {metrics['ssim']:.4f} ({label}), "
                  f"central = {metrics_central['ssim']:.4f}, "
                  f"eff = {metrics['efficiency']:.4f}, "
                  f"time = {elapsed:.1f}s")
            if electrostatic_report is not None:
                print(f"    Electrostatic gate: "
                      f"{'PASS' if electrostatic_report.ok else 'FAIL'}, "
                      f"|V|max = {electrostatic_report.max_voltage_V:.3e} V, "
                      f"theta_max = "
                      f"{electrostatic_report.max_deflection_rad:.3e} rad")

        return {
            'phase_screen':   raw['phase_screen'],
            'P_achieved':     raw['achieved'],
            'P_achieved_central': central_np,
            'P_target':       target_smooth,
            'P_target_raw':   target_pattern,
            'convergence':    raw['convergence'],
            'metrics':        metrics,
            'metrics_central': metrics_central,
            'polychromatic':  self.holo_solver.is_polychromatic,
            'phase_response': self.holo_solver.phase_response.name,
            'phase_scales': self.holo_solver.ensemble_phase_scales.copy(),
            'electrostatic_validity': electrostatic_report,
            'method':         method,
            'time_s':         elapsed,
            'phi_loops':      raw.get('phi_loops', None),
            'phi_loops_control': raw.get('phi_loops_control', None),
            'phi_loops_wrapped': raw.get('phi_loops_wrapped', None),
        }

    # -------------------------------------------------------------------
    # Stage 3 (optional): Floquet dressing + binding filter
    # -------------------------------------------------------------------
    def floquet_filter(self, density_2d, phase_screen, verbose=True):
        """Apply Floquet dressing and sideband-selective binding."""
        if not self.use_floquet:
            if verbose:
                print("\n  STAGE 3: Floquet — SKIPPED (disabled)")
            return density_2d, {}

        if verbose:
            print("\n  STAGE 3: Floquet dressing + binding filter")

        psi = np.sqrt(np.maximum(density_2d, 0)) * np.exp(
            1j * phase_screen)

        psi_d, avg_pops, entropy, V_map = \
            self.substrate.stage2_floquet_dress_spatial(
                psi, phase_screen, V_max=self.V_max_floquet,
                verbose=verbose)

        psi_ads, ads_frac, weights = \
            self.substrate.stage3_binding_filter(
                psi_d, avg_pops, n_resonant=self.n_resonant,
                verbose=verbose)

        return np.abs(psi_ads)**2, {
            'avg_pops': avg_pops,
            'entropy': entropy,
            'ads_frac': ads_frac,
            'V_map': V_map,
        }

    # -------------------------------------------------------------------
    # Stage 4: Diamond lattice AB caging
    # -------------------------------------------------------------------
    def cage_pattern(self, density_2d, phi_cage=np.pi,
                     T_evolve=40.0, verbose=True):
        """Evolve the deposited density on a 2D diamond lattice."""
        if verbose:
            print(f"\n  STAGE 4: Diamond caging "
                  f"(Φ = {phi_cage/np.pi:.1f}π)")

        psi0 = self.mapper.map_to_a_sites(density_2d, self.diamond)
        return self.diamond.evolve(
            psi0, phi=phi_cage, T=T_evolve, verbose=verbose)

    # -------------------------------------------------------------------
    # Full pipeline
    # -------------------------------------------------------------------
    def run(self, target_pattern, method='gd',
            run_caging=True, T_evolve=40.0,
            seed=42, verbose=True):
        """Execute the full v10 pipeline."""
        if verbose:
            self.info()

        # Stage 1: CMWB He⁺ source
        psi_source, sync_result = self.generate_source(
            seed=seed, verbose=verbose)

        # Build deterministic beam profile (transverse envelope, no noise).
        # No exp(i*k0*X) factor: forward momentum is carried by the
        # propagator's kz (fable5 E1).
        x = np.linspace(-self.L / 2, self.L / 2, self.N)
        X, Y = np.meshgrid(x, x, indexing='ij')
        sigma = 0.35 * self.L
        psi_profile = np.exp(-(X**2 + Y**2)
                             / (2 * sigma**2)).astype(np.complex128)
        psi_profile /= np.sqrt(np.sum(np.abs(psi_profile)**2) * self.dx**2)

        # Solver optimizes against the deterministic profile
        psi_profile_t = torch.tensor(psi_profile, dtype=torch.complex128,
                                      device=_device)
        self.holo_solver.psi_in_t = psi_profile_t
        self.holo_solver.A_in_t = torch.abs(psi_profile_t)
        self.holo_solver.psi_in_np = psi_profile
        self.holo_solver.A_in_np = np.abs(psi_profile)

        # Stage 2: Inverse holography
        holo = self.solve_holography(
            target_pattern, method=method, verbose=verbose)

        density_holo = holo['P_achieved']
        density_longitudinal = density_holo.copy()
        metrics_longitudinal = holo['metrics']
        # The batched kernels are no longer needed after optimization. Free
        # them before the stochastic evaluation, which streams one
        # wavelength propagator at a time to bound device memory.
        self.holo_solver.release_wavelength_ensemble()

        # Evaluate actual deposition (fable5 E2/T6/T23): partial coherence
        # is an ensemble property. Sequential ions sample both the source
        # wavelength distribution and independent transverse phase noise;
        # their intensities, not amplitudes, add at the target.
        screen_t = torch.tensor(holo['phase_screen'], dtype=torch.float64,
                                 device=_device)

        params = self._source_params
        sigma_theta = params.sigma_theta
        wavelengths = self._wavelength_samples
        if wavelengths is None:
            wavelengths = np.array([self.lam])
        Q = len(wavelengths)
        phase_scales = self.holo_solver.phase_scales(wavelengths)

        # Balance independent transverse-noise realizations across the
        # equal-weight wavelength quadrature. This evaluates the joint
        # incoherent source distribution with approximately the old total
        # propagation count, rather than multiplying M by Q.
        M_requested = self.n_noise_realizations if sigma_theta > 0 else 1
        phase_per_wavelength = (
            max(1, int(np.ceil(M_requested / Q))) if sigma_theta > 0 else 1
        )
        M = Q * phase_per_wavelength
        rng = np.random.default_rng(seed)
        xi_pix = params.xi_perp / self.dx

        density_actual = np.zeros((self.N, self.N))
        rms_applied = []

        for wavelength, phase_scale in zip(wavelengths, phase_scales):
            # All spectral components cross the same physical source-to-
            # target distance. Only k0 changes with wavelength.
            propagator = self.holo_solver.make_propagator(wavelength)
            T_screen = SQUIDArray.phase_screen_to_transmission(
                screen_t * float(phase_scale))

            for _ in range(phase_per_wavelength):
                noise = sample_phase_noise(
                    self.N, sigma_theta, xi_pix, rng)
                rms_applied.append(noise.std())
                psi_m = psi_profile * np.exp(1j * noise)
                psi_m /= np.sqrt(
                    np.sum(np.abs(psi_m)**2) * self.dx**2)
                psi_m_t = torch.tensor(
                    psi_m, dtype=torch.complex128, device=_device)
                psi_out = self.holo_solver._propagate_torch(
                    psi_m_t * T_screen, forward=True,
                    propagator=propagator)
                density_actual += torch.abs(psi_out).cpu().numpy()**2

            del propagator

        density_actual /= M
        rms_applied = float(np.mean(rms_applied))
        metrics_actual = compute_metrics(density_actual, holo['P_target'])

        # Transport attenuation (fable5 E6/T8): survival of the beam over
        # the cavity-exit → substrate leg at the ambient pressure, using
        # the Langevin mfp.  Reported separately — it does not rescale
        # the (peak-normalized) pattern metrics.
        L_transport = self.prop_distance_lam * self.lam
        transport_surv = self.cmwb.transport_survival(L_transport)

        # Space-charge gate (fable5 T10/M1): the single-particle picture
        # requires ≤ 1 ion in flight at a time.
        space_charge = self.cmwb.space_charge_check(
            path_length=L_transport, beam_radius=0.35 * self.L,
            current_A=params.current_A, verbose=verbose)

        if verbose:
            print(f"\n  Ensemble deposition: Q = {Q} wavelengths, "
                  f"M_θ = {phase_per_wavelength} per wavelength, "
                  f"M_total = {M}")
            print(f"  Longitudinal: Δλ/λ nominal = "
                  f"{params.dlam_frac:.3%}, applied RMS = "
                  f"{fractional_rms(wavelengths):.3%}, "
                  f"range = {wavelengths.min()*1e9:.3f}–"
                  f"{wavelengths.max()*1e9:.3f} nm")
            print(f"  Transverse: σ_θ nominal = {sigma_theta:.4f} rad, "
                  f"applied RMS = {rms_applied:.4f} rad")
            print(f"  Transport survival over {L_transport*1e9:.0f} nm at "
                  f"P = {self.cmwb.cavity.pressure:.3g} Pa: "
                  f"{transport_surv:.3e} "
                  f"(mfp = {self.cmwb.mfp_beam:.3e} m)")

        # Stage 3: Optional Floquet
        density_post_floquet, floquet_info = self.floquet_filter(
            density_actual, holo['phase_screen'], verbose=verbose)

        # Stage 4: Caging
        cage_pi, cage_0 = None, None
        if run_caging:
            cage_pi = self.cage_pattern(
                density_post_floquet, phi_cage=np.pi,
                T_evolve=T_evolve, verbose=verbose)
            cage_0 = self.cage_pattern(
                density_post_floquet, phi_cage=0.0,
                T_evolve=T_evolve, verbose=verbose)

        results = {
            'psi_source':       psi_source,
            'sync_result':      sync_result,
            'r_source':         sync_result['r_final'],
            'holo':             holo,
            'density_holo':     density_holo,
            'density_longitudinal': density_longitudinal,
            'metrics_longitudinal': metrics_longitudinal,
            'density_actual':   density_actual,
            'metrics_actual':   metrics_actual,
            'density_filtered': density_post_floquet,
            'floquet_info':     floquet_info,
            'cage_pi':          cage_pi,
            'cage_0':           cage_0,
            'sigma_theta':      sigma_theta,
            'noise_rms_applied': rms_applied,
            'M_realizations':   M,
            'M_phase_per_wavelength': phase_per_wavelength,
            'n_wavelength_samples': Q,
            'wavelength_samples_nm': wavelengths * 1e9,
            'dlam_frac_nominal': params.dlam_frac,
            'dlam_frac_applied': fractional_rms(wavelengths),
            'phase_actuator': self.holo_solver.phase_response.name,
            'phase_scales': phase_scales,
            'electrostatic_validity': holo.get('electrostatic_validity'),
            'source_params':    params,
            'space_charge':     space_charge,
            'transport_survival': transport_surv,
            'mfp_beam_m':       self.cmwb.mfp_beam,
            'lam_nm':           self.lam * 1e9,
            'pressure_Pa':      self.cmwb.cavity.pressure,
            'B_gauss':          self.cmwb.B * 1e4,
            'K_dim':            self.cmwb.K_dim,
            'n_eff':            self.cmwb.n_eff,
            'sync_mode':        sync_result.get('mode', 'unknown'),
            'N_per_volume':     self.cmwb.N_per_volume,
        }

        if verbose:
            self._print_summary(results)

        return results

    def _print_summary(self, r):
        print("\n" + "=" * 65)
        print("v10 PIPELINE SUMMARY")
        print("=" * 65)
        print(f"  Source:  He⁺, r = {r['r_source']:.4f}")
        sync_mode = r.get('sync_mode', 'unknown')
        print(f"  CMWB:   B = {r['B_gauss']:.0f} G, "
              f"P = {r['pressure_Pa']:.0f} Pa, "
              f"K_dim = {r['K_dim']:.3f}, "
              f"n_eff = {r['n_eff']:.0f}, "
              f"mode = {sync_mode}")
        m_clean = r['holo'].get('metrics_central', r['holo']['metrics'])
        m_long = r.get('metrics_longitudinal', r['holo']['metrics'])
        m_actual = r.get('metrics_actual', m_long)
        print(f"  Holo:   SSIM = {m_clean['ssim']:.4f} (clean), "
              f"{m_long['ssim']:.4f} (longitudinal), "
              f"{m_actual['ssim']:.4f} "
              f"(joint, M={r.get('M_realizations', 1)}), "
              f"gap = {m_clean['ssim'] - m_actual['ssim']:.4f}")
        print(f"  Noise:  σ_θ = {r.get('sigma_theta', 0):.4f} rad nominal, "
              f"{r.get('noise_rms_applied', 0):.4f} rad applied")
        print(f"  Spectrum: Δλ/λ = "
              f"{r.get('dlam_frac_applied', 0):.3%} RMS, "
              f"Q = {r.get('n_wavelength_samples', 1)}")
        scales = np.asarray(r.get('phase_scales', [1.0]))
        print(f"  Actuator: {r.get('phase_actuator', 'achromatic')}, "
              f"phase scale = {scales.min():.6f}–{scales.max():.6f}")
        electrostatic = r.get('electrostatic_validity')
        if electrostatic is not None:
            print(f"  Electrostatic thin screen: "
                  f"{'PASS' if electrostatic.ok else 'FAIL'}, "
                  f"|qV|/E = {electrostatic.max_energy_ratio:.3e}, "
                  f"walkoff/pitch = "
                  f"{electrostatic.max_walkoff_pitch_ratio:.3e}")
        print(f"  Transport survival: {r.get('transport_survival', 1):.3e}")
        sc = r.get('space_charge')
        if sc is not None:
            print(f"  Space charge: {'PASS' if sc['ok'] else 'FAIL'} "
                  f"({sc['N_in_flight']:.3g} ions in flight; "
                  f"single-ion limit {sc['I_max_single_A']:.3e} A)")
        if r['floquet_info']:
            fi = r['floquet_info']
            print(f"  Floquet: S = {fi.get('entropy', 'N/A')}, "
                  f"ads = {fi.get('ads_frac', 'N/A')}")
        if r['cage_pi'] is not None:
            loc_pi = r['cage_pi']['localization'][-1]
            loc_0 = r['cage_0']['localization'][-1]
            print(f"  Caging: Φ=π loc = {loc_pi:.4f}, "
                  f"Φ=0 loc = {loc_0:.4f}, "
                  f"gap = {loc_pi - loc_0:.4f}")


# ===========================================================================
# SINGLE-ION STATISTICAL ACCUMULATION (fable5 T11/M1)
# ===========================================================================

HolographicCagingPipeline = IntegratedPipelineV10


def dose_fidelity_curve(density, target, doses=None, seed=0, n_repeats=3):
    """Dose-vs-fidelity for one-ion-at-a-time deposition.

    With space charge forbidding co-resident ions (T10), the pattern
    accumulates statistically: each ion arrives independently and lands
    with probability density ∝ the ensemble-averaged |ψ|² (which already
    folds in the T6 noise ensemble — for sequential ions the intensity
    average is exact, not an approximation).  Arrivals are Poissonian:
    counts per pixel ~ Poisson(D · p).

    Parameters
    ----------
    density : ndarray (N, N)
        Ensemble-averaged arrival intensity (e.g. `density_actual`).
    target : ndarray (N, N)
        Target pattern the SSIM is scored against.
    doses : array-like or None
        Mean total ion numbers D to evaluate.  Default: log-spaced
        10²…10⁸.
    seed : int
        RNG seed.
    n_repeats : int
        Independent shot-noise realizations averaged per dose.

    Returns
    -------
    dict with 'doses', 'ssim_mean', 'ssim_std', and 'ssim_ceiling'
    (the infinite-dose ensemble value the curve saturates to).
    """
    if doses is None:
        doses = np.logspace(2, 8, 13)
    doses = np.asarray(doses, dtype=float)
    rng = np.random.default_rng(seed)
    p = density / density.sum()

    ssim_mean, ssim_std = [], []
    for D in doses:
        vals = [ssim_score(target, rng.poisson(D * p).astype(float))
                for _ in range(n_repeats)]
        ssim_mean.append(np.mean(vals))
        ssim_std.append(np.std(vals))

    return {
        'doses':        doses,
        'ssim_mean':    np.array(ssim_mean),
        'ssim_std':     np.array(ssim_std),
        'ssim_ceiling': float(ssim_score(target, density)),
    }


def dose_to_ssim(curve, ssim_target):
    """Smallest dose whose mean SSIM reaches ssim_target.

    Log-interpolates the dose–fidelity curve; returns None if the target
    is never reached (i.e. it exceeds the shot-noise-free ceiling).
    """
    d, s = curve['doses'], curve['ssim_mean']
    above = np.where(s >= ssim_target)[0]
    if len(above) == 0:
        return None
    i = above[0]
    if i == 0:
        return float(d[0])
    # interpolate in log-dose between the bracketing points
    f = (ssim_target - s[i - 1]) / (s[i] - s[i - 1] + 1e-30)
    return float(10 ** (np.log10(d[i - 1])
                        + f * (np.log10(d[i]) - np.log10(d[i - 1]))))


def plot_dose_fidelity(curve, marks=(), fname='results/v10_dose_fidelity.png'):
    """Plot SSIM vs ion dose with the ensemble ceiling."""
    import matplotlib.pyplot as plt

    print(f"\n  Plotting → {fname}")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(curve['doses'], curve['ssim_mean'], yerr=curve['ssim_std'],
                fmt='bo-', lw=2, capsize=3, label='Poisson accumulation')
    ax.axhline(curve['ssim_ceiling'], color='k', ls='--', alpha=0.6,
               label=f"ensemble ceiling = {curve['ssim_ceiling']:.3f}")
    for label, dose in marks:
        if dose is not None:
            ax.axvline(dose, color='g', ls=':', alpha=0.7)
            ax.annotate(label, (dose, 0.05), fontsize=8, rotation=90,
                        textcoords='offset points', xytext=(4, 0))
    ax.set_xscale('log')
    ax.set_xlabel('Ion dose (total arrivals)')
    ax.set_ylabel('SSIM vs target')
    ax.set_title('Single-ion statistical accumulation (T11)',
                 fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close(fig)


# ===========================================================================
# PRESSURE SWEEP
# ===========================================================================

def pressure_sweep(target_pattern, pressures=None,
                   B_field=0.01, seed=42, verbose=True):
    """Sweep CMWB cavity pressure to find the coherence threshold
    for usable holographic deposition."""
    if pressures is None:
        pressures = [1e-3, 1.0, 100, 1000, 10000, 101325]

    results = []
    for P in pressures:
        if verbose:
            print(f"\n{'#' * 65}")
            print(f"# PRESSURE = {P:.1e} Pa")
            print(f"{'#' * 65}")

        pipe = IntegratedPipelineV10(
            B_field=B_field, pressure_Pa=P,
        )
        r = pipe.run(target_pattern, method='gd',
                     run_caging=False, seed=seed, verbose=verbose)
        r['pressure'] = P
        results.append(r)

    return results


# ===========================================================================
# PLOTTING
# ===========================================================================

def plot_v10_pipeline(pipeline, results, fname='results/v10_pipeline.png'):
    """Comprehensive pipeline figure."""
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    print(f"\n  Plotting → {fname}")

    fig = plt.figure(figsize=(24, 20))
    gs = GridSpec(4, 4, figure=fig, hspace=0.45, wspace=0.35,
                  left=0.05, right=0.97, top=0.95, bottom=0.03)

    x_nm = np.linspace(-pipeline.L / 2, pipeline.L / 2,
                       pipeline.N) * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    # Row 0: Banner
    ax = fig.add_subplot(gs[0, :])
    ax.axis('off')
    m = results['holo']['metrics']
    cage_str = ''
    if results['cage_pi'] is not None:
        loc_pi = results['cage_pi']['localization'][-1]
        loc_0 = results['cage_0']['localization'][-1]
        cage_str = (f"\nCaging: Φ=π loc={loc_pi:.3f}  "
                    f"Φ=0 loc={loc_0:.3f}")
    txt = (
        "INTEGRATED QUANTUM SUBSTRATE  v10\n"
        "He⁺ CMWB → SQUID holography → diamond cage\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Source: He⁺ r={results['r_source']:.3f}  "
        f"B={results['B_gauss']:.0f}G  "
        f"P={results['pressure_Pa']:.0f}Pa  "
        f"K={results['K_dim']:.3f}\n"
        f"λ={pipeline.lam*1e9:.1f}nm  "
        f"Holo SSIM={m['ssim']:.3f}  "
        f"eff={m['efficiency']:.3f}"
        f"{cage_str}"
    )
    ax.text(0.5, 0.5, txt, transform=ax.transAxes,
            ha='center', va='center',
            fontsize=11, fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#e3f2fd',
                      edgecolor='#1565c0', linewidth=2))

    # Row 1: Target, achieved, phase screen, source density
    panels = [
        ('Target (smoothed)', results['holo']['P_target']),
        ('Holographic output', results['density_holo']),
        ('Phase screen', None),
        ('Source |ψ|²', np.abs(results['psi_source'])**2),
    ]
    for idx, (title, data) in enumerate(panels):
        ax = fig.add_subplot(gs[1, idx])
        if title == 'Phase screen':
            ph = results['holo']['phase_screen']
            im = ax.imshow(ph.T, extent=ext, cmap='twilight_shifted',
                           origin='lower', vmin=-np.pi, vmax=np.pi)
            ax.set_title(title, fontsize=10, fontweight='bold')
            plt.colorbar(im, ax=ax, shrink=0.8, label='rad')
        else:
            d = data / (data.max() + 1e-30)
            im = ax.imshow(d.T, extent=ext, cmap='inferno',
                           origin='lower')
            C = michelson_contrast(data)
            ax.set_title(f'{title}\nC={C:.3f}',
                         fontsize=10, fontweight='bold')
            plt.colorbar(im, ax=ax, shrink=0.8)
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)

    # Row 2: Convergence and cross-sections
    ax = fig.add_subplot(gs[2, 0:2])
    conv = results['holo']['convergence']
    ax.semilogy(conv, 'r-', lw=1.5)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Loss')
    ax.set_title('Holography convergence', fontsize=10,
                 fontweight='bold')
    ax.grid(True, alpha=0.3)

    ax = fig.add_subplot(gs[2, 2:4])
    mid = pipeline.N // 2
    target_line = results['holo']['P_target'][:, mid]
    achieved_line = results['density_holo'][:, mid]
    target_line /= (target_line.max() + 1e-30)
    achieved_line /= (achieved_line.max() + 1e-30)
    ax.plot(x_nm, target_line, 'k--', lw=2, label='Target')
    ax.plot(x_nm, achieved_line, 'r-', lw=2, label='Achieved')
    if results['floquet_info']:
        filt_line = results['density_filtered'][:, mid]
        filt_line /= (filt_line.max() + 1e-30)
        ax.plot(x_nm, filt_line, 'b-', lw=1.5,
                label='Post-Floquet', alpha=0.7)
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('Normalised intensity')
    ax.set_title('Cross-sections', fontsize=10, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Row 3: Caging
    if results['cage_pi'] is not None:
        cage_pi = results['cage_pi']
        cage_0 = results['cage_0']

        ax = fig.add_subplot(gs[3, 0])
        ax.plot(cage_pi['times'], cage_pi['localization'],
                'b-', lw=2,
                label=f"Φ=π loc={cage_pi['localization'][-1]:.3f}")
        ax.plot(cage_0['times'], cage_0['localization'],
                'r-', lw=2,
                label=f"Φ=0 loc={cage_0['localization'][-1]:.3f}")
        ax.set_xlabel('Time (ℏ/J)')
        ax.set_ylabel('Localization')
        ax.set_title('Diamond AB caging', fontsize=10,
                     fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)

        ax = fig.add_subplot(gs[3, 1])
        ax.plot(cage_pi['times'], cage_pi['fidelity'],
                'b-', lw=2, label='Φ=π')
        ax.plot(cage_0['times'], cage_0['fidelity'],
                'r-', lw=2, label='Φ=0')
        ax.set_xlabel('Time (ℏ/J)')
        ax.set_ylabel('Fidelity')
        ax.set_title('Fidelity evolution', fontsize=10,
                     fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        for si, (lbl, snaps) in enumerate([
                ('Caged t=0', cage_pi['snapshots'][0]),
                ('Caged t=T', cage_pi['snapshots'][-1])]):
            ax = fig.add_subplot(gs[3, 2 + si])
            _, d_sn = snaps
            im = ax.imshow(d_sn.T, cmap='hot', origin='lower',
                           aspect='auto')
            ax.set_title(f'Φ=π {lbl}', fontsize=10,
                         fontweight='bold')
            ax.set_xlabel('Cell x')
            ax.set_ylabel('Cell y')
            plt.colorbar(im, ax=ax, shrink=0.8)
    else:
        ax = fig.add_subplot(gs[3, :])
        ax.axis('off')
        ax.text(0.5, 0.5, 'Caging stage not run',
                ha='center', va='center', fontsize=14,
                transform=ax.transAxes, color='gray')

    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_multi_target(pipeline, all_results,
                      fname='results/v10_multi_target.png'):
    """Summary figure comparing multiple targets."""
    import matplotlib.pyplot as plt

    print(f"\n  Plotting → {fname}")
    names = list(all_results.keys())
    n = len(names)

    fig, axes = plt.subplots(3, n, figsize=(5 * n, 14))
    if n == 1:
        axes = axes.reshape(-1, 1)

    x_nm = np.linspace(-pipeline.L / 2, pipeline.L / 2,
                       pipeline.N) * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    for j, name in enumerate(names):
        r = all_results[name]
        m = r['holo']['metrics']

        ax = axes[0, j]
        P = r['holo']['P_target']
        ax.imshow(P.T / (P.max() + 1e-30), extent=ext,
                  cmap='inferno', origin='lower')
        ax.set_title(f'{name}\n(target)', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)

        ax = axes[1, j]
        I = r['density_holo']
        ax.imshow(I.T / (I.max() + 1e-30), extent=ext,
                  cmap='inferno', origin='lower')
        ax.set_title(f"SSIM={m['ssim']:.3f}\n"
                     f"eff={m['efficiency']:.3f}",
                     fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)

        ax = axes[2, j]
        ph = r['holo']['phase_screen']
        ax.imshow(ph.T, extent=ext, cmap='twilight_shifted',
                  origin='lower', vmin=-np.pi, vmax=np.pi)
        ax.set_title('Phase screen', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)

    fig.suptitle(f'v10 Multi-Target — He⁺ CMWB '
                 f'(r={all_results[names[0]]["r_source"]:.3f})',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


def plot_pressure_sweep(sweep_results,
                        fname='results/v10_pressure_sweep.png'):
    """Plot source coherence and holographic SSIM vs cavity pressure."""
    import matplotlib.pyplot as plt

    print(f"\n  Plotting → {fname}")

    pressures = [r['pressure'] for r in sweep_results]
    rs = [r['r_source'] for r in sweep_results]
    ssims = [r['holo']['metrics']['ssim'] for r in sweep_results]
    ssims_act = [r.get('metrics_actual', r['holo']['metrics'])['ssim']
                 for r in sweep_results]
    effs = [r['holo']['metrics']['efficiency'] for r in sweep_results]
    K_dims = [r['K_dim'] for r in sweep_results]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    ax = axes[0]
    ax.semilogx(pressures, rs, 'bo-', lw=2, markersize=8)
    ax.set_xlabel('Cavity pressure (Pa)')
    ax.set_ylabel('Source coherence r')
    ax.set_title('CMWB He⁺ synchronization', fontweight='bold')
    ax.axhline(0.95, color='g', ls='--', alpha=0.5, label='r = 0.95')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_ylim(-0.05, 1.05)
    for p, r, k in zip(pressures, rs, K_dims):
        ax.annotate(f'K={k:.2f}', (p, r), fontsize=7,
                    textcoords='offset points', xytext=(5, 5))

    ax = axes[1]
    ax.semilogx(pressures, ssims, 'rs-', lw=2, markersize=8,
                label='clean')
    ax.semilogx(pressures, ssims_act, 'b^--', lw=2, markersize=8,
                label='actual (ensemble)')
    ax.set_xlabel('Cavity pressure (Pa)')
    ax.set_ylabel('Holographic SSIM')
    ax.set_title('Deposition fidelity vs pressure', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.legend()

    ax = axes[2]
    ax.plot(rs, ssims_act, 'ko-', lw=2, markersize=8)
    ax.set_xlabel('Source coherence r')
    ax.set_ylabel('Holographic SSIM')
    ax.set_title('SSIM vs source coherence', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylabel('Holographic SSIM (actual)')
    for p, r, s in zip(pressures, rs, ssims_act):
        label = f'{p:.0e}Pa' if p >= 1 else f'{p:.0e}'
        ax.annotate(label, (r, s), fontsize=7,
                    textcoords='offset points', xytext=(5, 5))

    fig.suptitle('v10: Cavity Pressure → Source Coherence → '
                 'Holographic Fidelity',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    os.makedirs('results', exist_ok=True)

    print("╔═══════════════════════════════════════════════════════════╗")
    print("║  INTEGRATED QUANTUM SUBSTRATE  v10                       ║")
    print("║  He⁺ CMWB → SQUID holography → diamond caging           ║")
    print("╚═══════════════════════════════════════════════════════════╝")

    np.random.seed(0)
    N, L = 256, 400e-9

    # ==================================================================
    # VALIDATION GATES
    # ==================================================================
    print("\n" + "=" * 65)
    print("STEP 0: VALIDATION GATES")
    print("=" * 65)

    DiamondNetwork(Lx=4, Ly=4).validate_spectrum(
        phi=np.pi, abort_on_fail=True)

    # Kuramoto mode gate (fable5 T9): analytic and ODE branches must
    # agree at a size where both are affordable, above and below
    # threshold.  mode='auto' may then switch numerical method without
    # switching theory.
    cmwb_gate = CoherentMatterwaveBeam(
        species='He+', E_kinetic_eV=k_B * 1e-3 / e_C, B_field=0.01,
        cavity=CavityGeometry(pressure_Pa=1e-3), dE_frac=0.01)
    ok_sync, _, _ = cmwb_gate.validate_kuramoto_modes(N=500, tol=0.15)
    cmwb_gate_weak = CoherentMatterwaveBeam(
        species='He+', E_kinetic_eV=k_B * 1e-3 / e_C, B_field=0.0001,
        cavity=CavityGeometry(pressure_Pa=1e-3), dE_frac=0.01)
    ok_weak, _, _ = cmwb_gate_weak.validate_kuramoto_modes(N=500, tol=0.15)
    if not (ok_sync and ok_weak):
        sys.exit("ABORT: Kuramoto analytic/ODE validation gate failed.")

    squid_val = SQUIDArray(N_loops=32, N_grid=N, L_grid=L)
    solver_val = InverseHolographySolver(
        squid_array=squid_val, N=N, L=L, T_beam=1e-3,
        prop_distance_lam=20.0)
    rt_ok, rt_ssim = validate_roundtrip(solver_val)
    if not rt_ok:
        sys.exit("ABORT: Holography roundtrip validation failed.")

    # ==================================================================
    # DEMO 1: Full pipeline — vacuum He⁺ source, dots target
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 1: Full pipeline — vacuum He⁺, dots target")
    print("=" * 65)

    pipeline = IntegratedPipelineV10(
        B_field=0.01, pressure_Pa=100,  # rough vacuum
        T_beam=1e-3,
        N=N, L=L, N_loops=32,
        prop_distance_lam=20.0,
        use_floquet=False,
        N_diamond_x=16, N_diamond_y=16,
    )

    target_dots = target_grid_of_dots(N, L, n_dots=3)
    results_dots = pipeline.run(
        target_dots, method='gd', run_caging=True, seed=42)
    plot_v10_pipeline(pipeline, results_dots,
                      fname='results/v10_pipeline_dots.png')

    # ==================================================================
    # DEMO 2: Multi-target benchmark (same vacuum source)
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 2: Multi-target benchmark")
    print("=" * 65)

    targets = {
        'spot':   target_single_spot(N, L),
        'line':   target_line(N, L),
        'dots':   target_grid_of_dots(N, L, n_dots=3),
        'ring':   target_ring(N, L),
        'letter': target_letter(N, L, letter='H'),
    }

    all_results = {}
    for name, P in targets.items():
        print(f"\n{'#' * 65}")
        print(f"# TARGET: {name}")
        print(f"{'#' * 65}")
        all_results[name] = pipeline.run(
            P, method='gd', run_caging=True, seed=42)
    plot_multi_target(pipeline, all_results)

    # ==================================================================
    # DEMO 3: Pressure sweep — coherence threshold
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 3: Pressure sweep — source coherence threshold")
    print("=" * 65)

    sweep = pressure_sweep(
        target_dots,
        pressures=[1e-3, 1.0, 100, 1000, 10000, 101325],
        B_field=0.01, seed=42,
    )
    plot_pressure_sweep(sweep)

    # ==================================================================
    # DEMO 4: Single-ion accumulation with a direct-parameterized source
    # ==================================================================
    print("\n" + "=" * 65)
    print("DEMO 4: Single-ion statistical accumulation (T11)")
    print("        Direct source parameterization (T12), 0.1 pA (T10-safe)")
    print("=" * 65)

    # DirectSource bypasses the Kuramoto model: state the physical beam
    # parameters outright.  σ_θ = 0.079 rad is the Q-limited coherence
    # ceiling from the Phase 2 sweep; 0.1 pA keeps occupancy < 1.
    source_direct = DirectSource(
        dlam_frac=0.005, xi_perp=50e-9, current_A=0.1e-12,
        sigma_theta=0.079)
    pipeline_si = IntegratedPipelineV10(
        B_field=0.01, pressure_Pa=1e-3,  # transport-safe vacuum
        T_beam=1e-3, N=N, L=L, N_loops=32,
        prop_distance_lam=20.0, source=source_direct,
    )
    results_si = pipeline_si.run(
        target_dots, method='gd', run_caging=False, seed=42)

    curve = dose_fidelity_curve(
        results_si['density_actual'], results_si['holo']['P_target'],
        seed=42)
    d80 = dose_to_ssim(curve, 0.8)
    ceiling = curve['ssim_ceiling']
    d95c = dose_to_ssim(curve, 0.95 * ceiling)
    plot_dose_fidelity(curve, marks=[('95% of ceiling', d95c)])

    sc = results_si['space_charge']
    rate = sc['current_A'] / e_C  # singly charged
    print(f"\n  Dose–fidelity (dots target, σ_θ = 0.079 rad):")
    print(f"    Ensemble SSIM ceiling: {ceiling:.4f}")
    if d80 is not None:
        print(f"    Dose to SSIM 0.8: {d80:.3g} ions "
              f"({d80/rate:.3g} s at {sc['current_A']:.1e} A)")
    else:
        print(f"    Dose to SSIM 0.8: NOT REACHABLE — the shot-noise-free "
              f"ceiling is {ceiling:.3f} (aperture/coherence limited)")
    if d95c is not None:
        print(f"    Dose to 95% of ceiling ({0.95*ceiling:.3f}): "
              f"{d95c:.3g} ions ({d95c/rate:.3g} s at "
              f"{sc['current_A']:.1e} A)")

    # ==================================================================
    # SUMMARY
    # ==================================================================
    print("\n" + "=" * 65)
    print("ALL SIMULATIONS COMPLETE")
    print("=" * 65)
    print("\nOutput files:")
    print("  results/v10_pipeline_dots.png     — Full pipeline demo")
    print("  results/v10_multi_target.png      — 5 target benchmark")
    print("  results/v10_pressure_sweep.png    — Coherence threshold")
    print("  results/v10_dose_fidelity.png     — Single-ion dose curve")

    print(f"\nKey results (vacuum source, P = 100 Pa):")
    m = results_dots['holo']['metrics']
    print(f"  Source r:     {results_dots['r_source']:.4f}")
    print(f"  Holo SSIM:    {m['ssim']:.4f}")
    print(f"  Holo eff:     {m['efficiency']:.4f}")
    if results_dots['cage_pi'] is not None:
        loc = results_dots['cage_pi']['localization'][-1]
        print(f"  Cage loc:     {loc:.4f}")

    print(f"\nMulti-target SSIMs:")
    for name, res in all_results.items():
        m = res['holo']['metrics']
        print(f"  {name:<10}: SSIM={m['ssim']:.3f}  "
              f"eff={m['efficiency']:.3f}")

    print(f"\nPressure sweep:")
    print(f"  {'P (Pa)':<12} {'r':>6} {'σ_θ':>8} {'SSIM_cl':>8} "
          f"{'SSIM_act':>9} {'gap':>7} {'K_dim':>8} {'survival':>10}")
    print(f"  {'-'*74}")
    for r in sweep:
        m = r['holo']['metrics']
        ma = r.get('metrics_actual', m)
        print(f"  {r['pressure']:<12.1e} {r['r_source']:6.4f} "
              f"{r.get('sigma_theta', 0):8.4f} {m['ssim']:8.4f} "
              f"{ma['ssim']:9.4f} {m['ssim'] - ma['ssim']:7.4f} "
              f"{r['K_dim']:8.3f} {r.get('transport_survival', 1):10.3e}")


if __name__ == "__main__":
    main()
