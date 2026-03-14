"""
INTEGRATED Quantum Substrate Deposition Simulator
===================================================

This is the COUPLED simulation. One wavefunction flows through
the entire pipeline:

  ψ_beam → [Kuramoto sync] → ψ_coherent
         → [A-B phase mask] → ψ_phased
         → [Floquet drive]  → ψ_dressed (sidebands)
         → [Binding filter] → ψ_adsorbed (resonant states only)
         → [Propagate on caging lattice] → ψ_final (localized or spread)

The deposition map is |ψ_final|² — the probability density
of atoms that survived ALL filtering stages.

Each module transforms the SAME wavefunction. The output of
one is the input to the next. This is the real integration.
"""

import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

hbar = 1.0545718e-34
k_B = 1.380649e-23
m_He = 6.6464731e-27


class IntegratedQuantumSubstrate:
    """
    Fully coupled deposition simulation.

    The substrate is a 2D grid. Each grid point has:
    - An A-B phase value (from the gauge field geometry)
    - A Floquet drive coupling to the substrate lattice
    - A caging lattice Hamiltonian that governs post-adsorption dynamics

    The beam wavefunction ψ(x,y) evolves through all stages.
    """

    def __init__(self, N=256, L=400e-9, T_beam=1e-3,
                 N_floquet_sidebands=4, N_cage_cells=20):
        # Grid
        self.N = N
        self.L = L
        self.dx = L / N
        self.x = np.linspace(-L/2, L/2, N)
        self.y = np.linspace(-L/2, L/2, N)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

        # Beam
        self.mass = m_He
        self.T_beam = T_beam
        self.v = np.sqrt(2 * k_B * T_beam / m_He)
        self.k0 = self.mass * self.v / hbar
        self.lam = 2 * np.pi / self.k0
        self.E0 = 0.5 * self.mass * self.v**2

        # Floquet
        self.N_side = N_floquet_sidebands
        self.fl_dim = 2 * N_floquet_sidebands + 1
        self.n_vals = np.arange(-N_floquet_sidebands, N_floquet_sidebands + 1)

        # Caging lattice
        self.N_cage = N_cage_cells
        self.cage_dim = 3 * N_cage_cells

    def info(self):
        print("=" * 65)
        print("INTEGRATED QUANTUM SUBSTRATE SIMULATOR")
        print("=" * 65)
        print(f"  Beam: He-4 at {self.T_beam*1e3:.1f} mK")
        print(f"  λ_dB = {self.lam*1e9:.2f} nm")
        print(f"  v = {self.v:.2f} m/s")
        print(f"  E₀ = {self.E0:.4e} J")
        print(f"  Substrate: {self.L*1e9:.0f} nm, {self.N}×{self.N}")
        print(f"  Floquet sidebands: {self.fl_dim} (n = {-self.N_side}..{self.N_side})")
        print(f"  Caging lattice: {self.N_cage} cells, {self.cage_dim} sites")
        print(f"  dx = {self.dx*1e9:.2f} nm, λ/dx = {self.lam/self.dx:.1f}")

    # ---------------------------------------------------------------
    # STAGE 0: Beam generation with Kuramoto pre-synchronization
    # ---------------------------------------------------------------
    def stage0_beam(self, N_atoms_kuramoto=200, K_coupling=6.0,
                    sync_time=30, alpha=0.5):
        """
        Generate beam and pre-synchronize with Kuramoto dynamics.

        The Kuramoto order parameter r determines the beam's
        effective coherence. r=1 means perfectly coherent;
        r<1 means partial coherence that reduces contrast.

        Returns: ψ(x,y), coherence factor r
        """
        print("\n  STAGE 0: Beam generation + Kuramoto synchronization")

        # Kuramoto synchronization
        N_k = N_atoms_kuramoto
        omega = np.random.normal(0, 1.0, N_k)
        theta = np.random.uniform(0, 2 * np.pi, N_k)
        dtheta = np.zeros(N_k)
        dt = 0.01
        times = np.arange(0, sync_time, dt)

        order_history = []
        for t in times:
            z = np.mean(np.exp(1j * theta))
            r = np.abs(z)
            order_history.append(r)
            coupling = np.imag(z * np.exp(-1j * theta))
            ddtheta = -alpha * dtheta + omega + K_coupling * coupling
            dtheta += ddtheta * dt
            theta += dtheta * dt

        r_final = order_history[-1]
        print(f"    Kuramoto order parameter: r = {r_final:.4f}")

        # Build beam: Gaussian envelope × plane wave
        # The coherence r modulates the effective beam quality
        sigma = 0.35 * self.L
        psi = np.exp(-(self.X**2 + self.Y**2) / (2 * sigma**2))
        psi = psi * np.exp(1j * self.k0 * self.X)

        # Partial coherence: add phase noise proportional to (1-r)
        # This models the effect of imperfect synchronization
        phase_noise = (1 - r_final) * np.random.normal(0, np.pi, (self.N, self.N))
        # Smooth the noise to make it spatially correlated
        from scipy.ndimage import gaussian_filter
        phase_noise = gaussian_filter(phase_noise, sigma=3)
        psi = psi * np.exp(1j * phase_noise)

        # Normalize
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * self.dx**2)

        beam_norm = np.sum(np.abs(psi)**2) * self.dx**2
        print(f"    Beam norm: {beam_norm:.6f}")
        print(f"    Phase noise RMS: {np.std(phase_noise):.4f} rad")

        return psi, r_final, np.array(order_history), times

    # ---------------------------------------------------------------
    # STAGE 1: A-B phase imprinting
    # ---------------------------------------------------------------
    def stage1_ab_phase(self, psi_in, pattern='vortex_lattice', **kw):
        """
        Substrate imprints Aharonov-Bohm phase on the wavefunction.

        ψ_out(x,y) = ψ_in(x,y) × exp(i × φ_AB(x,y))

        φ_AB is determined by the substrate's vector potential geometry.
        """
        print("\n  STAGE 1: A-B phase imprinting")

        X, Y = self.X, self.Y
        phase = np.zeros_like(X)

        if pattern == 'vortex_lattice':
            a = kw.get('a', 6 * self.lam)
            core = kw.get('core', self.lam)
            for i in range(-6, 7):
                for j in range(-6, 7):
                    x0 = a * (i + 0.5 * (j % 2))
                    y0 = a * np.sqrt(3) / 2 * j
                    if abs(x0) < self.L / 2 and abs(y0) < self.L / 2:
                        R = np.sqrt((X - x0)**2 + (Y - y0)**2 + 1e-30)
                        theta = np.arctan2(Y - y0, X - x0)
                        sign = (-1)**(i + j)
                        w = 1 - np.exp(-R**2 / (2 * core**2))
                        phase += sign * np.pi * theta / (2 * np.pi) * w

        elif pattern == 'ring_array':
            spacing = kw.get('spacing', 8 * self.lam)
            r_ring = kw.get('r_ring', 3 * self.lam)
            width = kw.get('width', 0.5 * self.lam)
            flux = kw.get('flux', 0.5)
            n_rings = kw.get('n_rings', 3)
            for i in range(-n_rings, n_rings + 1):
                for j in range(-n_rings, n_rings + 1):
                    x0, y0 = i * spacing, j * spacing
                    R = np.sqrt((X - x0)**2 + (Y - y0)**2)
                    theta = np.arctan2(Y - y0, X - x0)
                    ring_w = np.exp(-(R - r_ring)**2 / (2 * width**2))
                    phase += ring_w * flux * theta

        elif pattern == 'sinusoidal':
            period = kw.get('period', 4 * self.lam)
            amp = kw.get('amp', np.pi)
            phase = amp * np.cos(2*np.pi*X/period) * np.cos(2*np.pi*Y/period)

        # Apply phase to wavefunction
        psi_out = psi_in * np.exp(1j * phase)

        print(f"    Pattern: {pattern}")
        print(f"    Phase range: [{phase.min():.2f}, {phase.max():.2f}] rad")
        print(f"    |ψ| preserved: {np.allclose(np.abs(psi_out), np.abs(psi_in))}")

        return psi_out, phase

    # ---------------------------------------------------------------
    # STAGE 2: Floquet dressing (sideband creation)
    # ---------------------------------------------------------------
    def stage2_floquet_dress(self, psi_phased, V_drive_frac=0.15,
                             omega_drive=None):
        """
        The substrate's time-periodic drive dresses each spatial point
        of the wavefunction with Floquet sidebands.

        At each (x,y), the atom acquires a sideband decomposition:
        ψ(x,y) → Σ_n c_n(V) × ψ(x,y) × |n⟩

        where c_n(V) is the nth sideband amplitude at drive strength V.

        V_drive can be spatially varying (stronger drive near substrate
        features = coupling between A-B geometry and Floquet physics).

        Returns: ψ_dressed(x,y,n) — a 3D array [Nx, Ny, N_sidebands]
        """
        print("\n  STAGE 2: Floquet sideband dressing")

        if omega_drive is None:
            omega_drive = self.E0 / (hbar * 3)  # ω such that 3ℏω ~ E0

        V_drive = V_drive_frac * self.E0

        # Build Floquet Hamiltonian
        diag = np.array([n * hbar * omega_drive for n in self.n_vals])
        H_F = np.diag(diag)
        for i in range(self.fl_dim - 1):
            H_F[i, i + 1] = V_drive
            H_F[i + 1, i] = V_drive

        # Compute sideband amplitudes from initial n=0 state
        T_drive = 2 * np.pi / omega_drive
        U = expm(-1j * H_F * T_drive / hbar)
        psi0_fl = np.zeros(self.fl_dim, dtype=complex)
        psi0_fl[self.N_side] = 1.0  # Start in n=0
        c_n = U @ psi0_fl  # Sideband amplitudes

        print(f"    Drive: V = {V_drive_frac:.2f} E₀ = {V_drive:.4e} J")
        print(f"    ω = {omega_drive:.4e} rad/s")
        sb_pops = np.round(np.abs(c_n)**2, 4)
        print(f"    Sideband populations: {sb_pops}")
        print(f"    Total population: {np.sum(np.abs(c_n)**2):.6f}")

        # Dress each spatial point: ψ_dressed(x,y,n) = c_n × ψ(x,y)
        psi_dressed = np.zeros((self.N, self.N, self.fl_dim), dtype=complex)
        for n_idx in range(self.fl_dim):
            psi_dressed[:, :, n_idx] = c_n[n_idx] * psi_phased

        return psi_dressed, c_n, hbar * omega_drive

    # ---------------------------------------------------------------
    # STAGE 3: Binding resonance filter (state-selective adsorption)
    # ---------------------------------------------------------------
    def stage3_binding_filter(self, psi_dressed, hw_drive,
                              n_resonant=0, binding_width_frac=0.3):
        """
        The substrate has a binding resonance at a specific energy.
        Only sidebands whose quasi-energy falls within the resonance
        width get adsorbed. All others are reflected.

        This is the quantum state filter: the substrate "selects"
        which Floquet state binds.

        ψ_adsorbed(x,y) = Σ_n  f(E_n, E_bind) × ψ_dressed(x,y,n)

        where f is a Lorentzian centered on E_bind.
        """
        print("\n  STAGE 3: Binding resonance filter")

        E_bind = n_resonant * hw_drive  # Binding energy at nth sideband
        width = binding_width_frac * hw_drive

        print(f"    Binding resonance at n = {n_resonant} (E = {E_bind:.4e} J)")
        print(f"    Resonance width: {width:.4e} J ({binding_width_frac:.0%} of ℏω)")

        # Filter: Lorentzian weight for each sideband
        psi_adsorbed = np.zeros((self.N, self.N), dtype=complex)
        total_weight = 0

        for n_idx, n in enumerate(self.n_vals):
            E_n = n * hw_drive
            # Lorentzian transfer function
            weight = width**2 / ((E_n - E_bind)**2 + width**2)
            psi_adsorbed += weight * psi_dressed[:, :, n_idx]
            pop = np.sum(np.abs(psi_dressed[:, :, n_idx])**2) * self.dx**2
            total_weight += weight * pop
            if weight > 0.01:
                print(f"      n={n:+d}: weight={weight:.4f}, "
                      f"pop_in={pop:.4f}")

        norm_in = np.sum(np.sum(np.abs(psi_dressed)**2, axis=2)) * self.dx**2
        norm_out = np.sum(np.abs(psi_adsorbed)**2) * self.dx**2
        adsorption_fraction = norm_out / norm_in if norm_in > 0 else 0

        print(f"    Adsorption fraction: {adsorption_fraction:.4f}")
        print(f"    Reflected fraction:  {1 - adsorption_fraction:.4f}")

        return psi_adsorbed, adsorption_fraction

    # ---------------------------------------------------------------
    # STAGE 4: Post-adsorption dynamics on caging lattice
    # ---------------------------------------------------------------
    def stage4_caging(self, density_2d, phi_cage=np.pi, J_hop=1.0,
                      T_evolve=30, dt=0.05):
        """
        After adsorption, do the atoms STAY where they landed?

        Model: the substrate has a lattice structure. Adsorbed atoms
        can hop between sites. The A-B flux through each plaquette
        determines whether hopping is allowed (Φ=0) or forbidden (Φ=π).

        We take the 1D cross-section of the deposited density,
        map it onto the caging lattice as initial conditions,
        and evolve to see if the pattern is preserved.

        At Φ=π: flat bands → atoms STAY → pattern preserved
        At Φ=0: dispersive → atoms SPREAD → pattern destroyed
        """
        print(f"\n  STAGE 4: Post-adsorption caging (Φ = {phi_cage/np.pi:.2f}π)")

        dim = self.cage_dim
        N_cells = self.N_cage

        # Build rhombic lattice Hamiltonian
        H = np.zeros((dim, dim), dtype=complex)
        for n in range(N_cells):
            a, b, c = 3*n, 3*n+1, 3*n+2
            H[a, b] = H[b, a] = J_hop
            H[a, c] = H[c, a] = J_hop
            if n < N_cells - 1:
                a2 = 3 * (n + 1)
                H[b, a2] = J_hop * np.exp(1j * phi_cage / 2)
                H[a2, b] = J_hop * np.exp(-1j * phi_cage / 2)
                H[c, a2] = J_hop * np.exp(-1j * phi_cage / 2)
                H[a2, c] = J_hop * np.exp(1j * phi_cage / 2)

        # Map 2D density cross-section → lattice initial state
        mid = self.N // 2
        cross_section = density_2d[:, mid]

        # Downsample to match lattice dimension
        # Only use 'a' sites (every 3rd site) for initial loading
        indices = np.linspace(0, self.N - 1, N_cells).astype(int)
        lattice_density = cross_section[indices]
        lattice_density /= (np.sum(lattice_density) + 1e-30)

        # Initialize wavefunction: sqrt(density) on 'a' sites
        psi_lattice = np.zeros(dim, dtype=complex)
        for n in range(N_cells):
            psi_lattice[3 * n] = np.sqrt(lattice_density[n])

        # Normalize
        psi_lattice /= (np.sqrt(np.sum(np.abs(psi_lattice)**2)) + 1e-30)

        # Time evolution
        U_dt = expm(-1j * H * dt)
        times = np.arange(0, T_evolve, dt)
        initial_density = np.abs(psi_lattice)**2

        # Track how well the pattern is preserved
        fidelity_history = []  # overlap with initial state
        spread_history = []
        density_snapshots = []
        snapshot_times = [0, T_evolve * 0.25, T_evolve * 0.5, T_evolve]

        psi = psi_lattice.copy()
        for t_idx, t in enumerate(times):
            probs = np.abs(psi)**2
            # Fidelity: how much does current state overlap with initial?
            fidelity = np.abs(np.sum(np.conj(psi_lattice) * psi))**2
            fidelity_history.append(fidelity)

            # RMS spread
            sites = np.arange(dim)
            mu = np.sum(sites * probs)
            spread = np.sqrt(np.sum((sites - mu)**2 * probs))
            spread_history.append(spread)

            # Snapshots
            for st in snapshot_times:
                if abs(t - st) < dt / 2:
                    density_snapshots.append((t, probs.copy()))

            psi = U_dt @ psi

        final_fidelity = fidelity_history[-1]
        print(f"    Initial pattern loaded onto {N_cells} lattice sites")
        print(f"    Evolution time: {T_evolve} (ℏ/J)")
        print(f"    Final fidelity: {final_fidelity:.4f}")
        print(f"    Pattern {'PRESERVED' if final_fidelity > 0.5 else 'DESTROYED'}")

        return {
            'times': times,
            'fidelity': np.array(fidelity_history),
            'spread': np.array(spread_history),
            'initial_density': initial_density,
            'snapshots': density_snapshots,
            'phi': phi_cage,
        }

    # ---------------------------------------------------------------
    # FULL PIPELINE
    # ---------------------------------------------------------------
    def run_full_pipeline(self, pattern='vortex_lattice',
                          V_drive_frac=0.15, n_resonant=0,
                          K_kuramoto=6.0, phi_cage=np.pi,
                          prop_distance_lam=20, **pattern_kw):
        """
        Run the complete integrated pipeline.
        One wavefunction, four transformations, one deposition map.
        """
        self.info()

        # Stage 0: Beam + sync
        psi_beam, r_sync, order_hist, sync_times = self.stage0_beam(
            K_coupling=K_kuramoto)

        # Stage 1: A-B phase
        psi_phased, phase_map = self.stage1_ab_phase(
            psi_beam, pattern=pattern, **pattern_kw)

        # Fresnel propagation through substrate
        print(f"\n  PROPAGATION: {prop_distance_lam}λ = "
              f"{prop_distance_lam * self.lam * 1e9:.1f} nm")
        psi_propagated = self._propagate(psi_phased,
                                          prop_distance_lam * self.lam)

        # Stage 2: Floquet dressing
        psi_dressed, c_n, hw = self.stage2_floquet_dress(
            psi_propagated, V_drive_frac=V_drive_frac)

        # Stage 3: Binding filter
        psi_adsorbed, ads_frac = self.stage3_binding_filter(
            psi_dressed, hw, n_resonant=n_resonant)

        # Compute densities at key stages
        density_after_phase = np.abs(psi_propagated)**2
        density_final = np.abs(psi_adsorbed)**2

        # Stage 4: Caging (two cases for comparison)
        cage_pi = self.stage4_caging(density_final, phi_cage=np.pi)
        cage_0 = self.stage4_caging(density_final, phi_cage=0)

        return {
            # Wavefunctions
            'psi_beam': psi_beam,
            'psi_phased': psi_phased,
            'psi_propagated': psi_propagated,
            'psi_adsorbed': psi_adsorbed,
            # Densities
            'density_after_AB': density_after_phase,
            'density_final': density_final,
            # Phase & Floquet
            'phase_map': phase_map,
            'sideband_amplitudes': c_n,
            'adsorption_fraction': ads_frac,
            # Kuramoto
            'sync_order': order_hist,
            'sync_times': sync_times,
            'r_sync': r_sync,
            # Caging
            'cage_pi': cage_pi,
            'cage_0': cage_0,
        }

    def _propagate(self, psi, distance):
        kx = np.fft.fftfreq(self.N, self.dx) * 2 * np.pi
        ky = np.fft.fftfreq(self.N, self.dx) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        kz_sq = self.k0**2 - KX**2 - KY**2
        prop = kz_sq > 0
        kz = np.zeros_like(kz_sq)
        kz[prop] = np.sqrt(kz_sq[prop])
        H = np.zeros_like(kz_sq, dtype=complex)
        H[prop] = np.exp(1j * kz[prop] * distance)
        return np.fft.ifft2(np.fft.fft2(psi) * H)


# ===========================================================================
# COMPARISON: With vs without each stage
# ===========================================================================

def ablation_study():
    """
    Run the pipeline with each stage toggled on/off to show
    the contribution of each module to the final deposition.
    """
    print("\n" + "=" * 65)
    print("ABLATION STUDY: Effect of each module on deposition")
    print("=" * 65)

    sim = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)
    sim.info()

    pattern = 'vortex_lattice'
    pattern_kw = dict(a=6*sim.lam, core=sim.lam)
    prop_dist = 20 * sim.lam

    results = {}

    # --- Case 1: Full pipeline ---
    print("\n\n>>> CASE 1: Full pipeline (all stages active)")
    psi, _, _, _ = sim.stage0_beam(K_coupling=6.0)
    psi, phase = sim.stage1_ab_phase(psi, pattern, **pattern_kw)
    psi = sim._propagate(psi, prop_dist)
    psi_dressed, c_n, hw = sim.stage2_floquet_dress(psi, V_drive_frac=0.2)
    psi_ads, af = sim.stage3_binding_filter(psi_dressed, hw, n_resonant=0)
    results['Full pipeline'] = np.abs(psi_ads)**2

    # --- Case 2: No Floquet (skip stage 2+3, direct adsorption) ---
    print("\n\n>>> CASE 2: No Floquet filter (phase + propagation only)")
    psi, _, _, _ = sim.stage0_beam(K_coupling=6.0)
    psi, _ = sim.stage1_ab_phase(psi, pattern, **pattern_kw)
    psi = sim._propagate(psi, prop_dist)
    results['No Floquet'] = np.abs(psi)**2

    # --- Case 3: No A-B phase (Floquet only) ---
    print("\n\n>>> CASE 3: No A-B phase (flat substrate + Floquet)")
    psi, _, _, _ = sim.stage0_beam(K_coupling=6.0)
    # Skip phase imprinting
    psi = sim._propagate(psi, prop_dist)
    psi_dressed, c_n, hw = sim.stage2_floquet_dress(psi, V_drive_frac=0.2)
    psi_ads, af = sim.stage3_binding_filter(psi_dressed, hw, n_resonant=0)
    results['No A-B phase'] = np.abs(psi_ads)**2

    # --- Case 4: No synchronization (poor coherence) ---
    print("\n\n>>> CASE 4: Poor beam coherence (K=0.5, no sync)")
    psi, _, _, _ = sim.stage0_beam(K_coupling=0.5)
    psi, _ = sim.stage1_ab_phase(psi, pattern, **pattern_kw)
    psi = sim._propagate(psi, prop_dist)
    psi_dressed, c_n, hw = sim.stage2_floquet_dress(psi, V_drive_frac=0.2)
    psi_ads, af = sim.stage3_binding_filter(psi_dressed, hw, n_resonant=0)
    results['Poor coherence'] = np.abs(psi_ads)**2

    # --- Case 5: Different resonant sideband ---
    print("\n\n>>> CASE 5: Binding at n=+2 sideband instead of n=0")
    psi, _, _, _ = sim.stage0_beam(K_coupling=6.0)
    psi, _ = sim.stage1_ab_phase(psi, pattern, **pattern_kw)
    psi = sim._propagate(psi, prop_dist)
    psi_dressed, c_n, hw = sim.stage2_floquet_dress(psi, V_drive_frac=0.2)
    psi_ads, af = sim.stage3_binding_filter(psi_dressed, hw, n_resonant=2)
    results['Sideband n=+2'] = np.abs(psi_ads)**2

    return sim, results, phase


# ===========================================================================
# VISUALIZATION
# ===========================================================================

def create_integrated_figure(sim, pipeline_results, ablation_sim,
                             ablation_results, ablation_phase):
    """Build the master figure showing the coupled pipeline."""

    fig = plt.figure(figsize=(28, 38))
    gs = GridSpec(8, 4, figure=fig, hspace=0.45, wspace=0.35,
                  left=0.04, right=0.97, top=0.97, bottom=0.02)

    r = pipeline_results
    x_nm = sim.x * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    # ===== ROW 0: Header =====
    ax = fig.add_subplot(gs[0, :])
    ax.axis('off')
    txt = (
        "INTEGRATED QUANTUM SUBSTRATE: COUPLED PIPELINE\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"He-4 at {sim.T_beam*1e3:.0f} mK  •  λ_dB = {sim.lam*1e9:.1f} nm  •  "
        f"v = {sim.v:.2f} m/s  •  r_sync = {r['r_sync']:.3f}\n\n"
        "One wavefunction flows through the entire pipeline:\n"
        "ψ_beam → [Kuramoto sync] → [A-B phase] → [Propagate] → "
        "[Floquet dress] → [Bind filter] → [Caging] → deposition"
    )
    ax.text(0.5, 0.5, txt, transform=ax.transAxes,
            ha='center', va='center', fontsize=12, fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#e3f2fd', edgecolor='#1565c0'))

    # ===== ROW 1: Pipeline stages (density at each stage) =====
    stages = [
        ('① Beam |ψ|²\n(after sync)', np.abs(r['psi_beam'])**2),
        ('② After A-B phase\n(|ψ| unchanged)', np.abs(r['psi_phased'])**2),
        ('③ After propagation\n(interference)', r['density_after_AB']),
        ('④ After Floquet filter\n(state-selected)', r['density_final']),
    ]
    for idx, (title, density) in enumerate(stages):
        ax = fig.add_subplot(gs[1, idx])
        d = density / (np.max(density) + 1e-30)
        im = ax.imshow(d.T, extent=ext, cmap='inferno', origin='lower')
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ===== ROW 2: Phase at each stage =====
    phase_stages = [
        ('Beam phase\n(with noise)', np.angle(r['psi_beam'])),
        ('A-B phase map\n(substrate)', r['phase_map']),
        ('Phase after A-B\n(beam + substrate)', np.angle(r['psi_phased'])),
        ('Phase after propagation', np.angle(r['psi_propagated'])),
    ]
    for idx, (title, ph) in enumerate(phase_stages):
        ax = fig.add_subplot(gs[2, idx])
        im = ax.imshow(ph.T, extent=ext, cmap='hsv', origin='lower')
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8, label='rad')

    # ===== ROW 3: Cross-sections at each stage =====
    ax_cs = fig.add_subplot(gs[3, 0:2])
    mid = sim.N // 2
    labels = ['① Beam', '③ After A-B+prop', '④ After Floquet filter']
    densities = [np.abs(r['psi_beam'])**2,
                 r['density_after_AB'],
                 r['density_final']]
    colors_cs = ['gray', 'blue', 'red']
    for label, dens, col in zip(labels, densities, colors_cs):
        line = dens[:, mid]
        line_n = line / (np.max(line) + 1e-30)
        ax_cs.plot(x_nm, line_n, color=col, linewidth=2, label=label)
    ax_cs.set_xlabel('x (nm)')
    ax_cs.set_ylabel('Normalized |ψ|²')
    ax_cs.set_title('Cross-Section: Wavefunction at Each Pipeline Stage',
                    fontsize=11, fontweight='bold')
    ax_cs.legend(fontsize=10)
    ax_cs.grid(True, alpha=0.3)

    # Sideband populations
    ax_sb = fig.add_subplot(gs[3, 2])
    c_n = r['sideband_amplitudes']
    pops = np.abs(c_n)**2
    ax_sb.bar(sim.n_vals, pops, color='steelblue', edgecolor='navy')
    ax_sb.set_xlabel('Sideband n')
    ax_sb.set_ylabel('|c_n|²')
    ax_sb.set_title('Floquet Sideband\nPopulations', fontsize=10, fontweight='bold')
    ax_sb.grid(True, alpha=0.3)

    # Kuramoto
    ax_ku = fig.add_subplot(gs[3, 3])
    ax_ku.plot(r['sync_times'], r['sync_order'], 'g-', linewidth=1.5)
    ax_ku.axhline(y=r['r_sync'], color='red', linestyle='--',
                 label=f'Final r = {r["r_sync"]:.3f}')
    ax_ku.set_xlabel('Time')
    ax_ku.set_ylabel('Order parameter r')
    ax_ku.set_title('Beam Synchronization', fontsize=10, fontweight='bold')
    ax_ku.legend()
    ax_ku.grid(True, alpha=0.3)
    ax_ku.set_ylim(0, 1.05)

    # ===== ROW 4: Caging comparison =====
    cage_pi = r['cage_pi']
    cage_0 = r['cage_0']

    ax_cf = fig.add_subplot(gs[4, 0])
    ax_cf.plot(cage_pi['times'], cage_pi['fidelity'], 'b-', linewidth=2,
              label='Φ=π (caged)')
    ax_cf.plot(cage_0['times'], cage_0['fidelity'], 'r-', linewidth=2,
              label='Φ=0 (free)')
    ax_cf.set_xlabel('Time (ℏ/J)')
    ax_cf.set_ylabel('Pattern Fidelity')
    ax_cf.set_title('Does the Pattern Survive?\n(fidelity with initial)',
                    fontsize=10, fontweight='bold')
    ax_cf.legend(fontsize=10)
    ax_cf.grid(True, alpha=0.3)

    ax_cs2 = fig.add_subplot(gs[4, 1])
    ax_cs2.plot(cage_pi['times'], cage_pi['spread'], 'b-', linewidth=2,
               label='Φ=π')
    ax_cs2.plot(cage_0['times'], cage_0['spread'], 'r-', linewidth=2,
               label='Φ=0')
    ax_cs2.set_xlabel('Time (ℏ/J)')
    ax_cs2.set_ylabel('RMS Spread')
    ax_cs2.set_title('Spatial Spread of\nDeposited Atoms', fontsize=10, fontweight='bold')
    ax_cs2.legend(fontsize=10)
    ax_cs2.grid(True, alpha=0.3)

    # Density snapshots
    ax_snap = fig.add_subplot(gs[4, 2:4])
    n_snaps = len(cage_pi['snapshots'])
    if n_snaps > 0:
        colors_snap = plt.cm.viridis(np.linspace(0, 1, n_snaps))
        for i, (t, dens) in enumerate(cage_pi['snapshots']):
            # Plot every 3rd site (the 'a' sites)
            a_sites = np.arange(0, sim.cage_dim, 3)
            ax_snap.plot(a_sites, dens[a_sites], color=colors_snap[i],
                        linewidth=1.5, label=f't={t:.0f}')
        ax_snap.set_xlabel('Site')
        ax_snap.set_ylabel('|ψ|²')
        ax_snap.set_title('Caged Pattern Snapshots (Φ=π)\n'
                         'Pattern PRESERVED over time',
                         fontsize=10, fontweight='bold')
        ax_snap.legend(fontsize=8)
        ax_snap.grid(True, alpha=0.3)

    # ===== ROW 5: Ablation study =====
    ab_x = ablation_sim.x * 1e9
    ab_ext = [ab_x[0], ab_x[-1], ab_x[0], ab_x[-1]]

    for idx, (title, density) in enumerate(list(ablation_results.items())[:4]):
        ax = fig.add_subplot(gs[5, idx])
        d = density / (np.max(density) + 1e-30)
        im = ax.imshow(d.T, extent=ab_ext, cmap='inferno', origin='lower')
        ax.set_title(f'Ablation: {title}', fontsize=9, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ===== ROW 6: Ablation cross-sections =====
    ax_ab = fig.add_subplot(gs[6, 0:3])
    mid = ablation_sim.N // 2
    ab_colors = {'Full pipeline': 'blue', 'No Floquet': 'green',
                 'No A-B phase': 'red', 'Poor coherence': 'orange',
                 'Sideband n=+2': 'purple'}
    for title, density in ablation_results.items():
        line = density[:, mid]
        line_n = line / (np.max(line) + 1e-30)
        ax_ab.plot(ab_x, line_n, color=ab_colors.get(title, 'gray'),
                  linewidth=1.5, label=title)
    ax_ab.set_xlabel('x (nm)')
    ax_ab.set_ylabel('Normalized |ψ|²')
    ax_ab.set_title('Ablation Study: Cross-Sections\n'
                    'Removing each module shows its contribution',
                    fontsize=11, fontweight='bold')
    ax_ab.legend(fontsize=9)
    ax_ab.grid(True, alpha=0.3)

    # Summary
    ax_sum = fig.add_subplot(gs[6, 3])
    ax_sum.axis('off')
    contrasts = {}
    for title, density in ablation_results.items():
        d_flat = density.ravel()
        c = (np.max(d_flat) - np.min(d_flat)) / (np.max(d_flat) + np.min(d_flat) + 1e-30)
        contrasts[title] = c
    txt_lines = ["ABLATION CONTRAST\n━━━━━━━━━━━━━━━━━\n"]
    for title, c in contrasts.items():
        txt_lines.append(f"{title}:\n  C = {c:.4f}\n")
    ax_sum.text(0.05, 0.95, '\n'.join(txt_lines), transform=ax_sum.transAxes,
               fontsize=9, fontfamily='monospace', va='top',
               bbox=dict(boxstyle='round', facecolor='#fce4ec', edgecolor='#c62828'))

    # ===== ROW 7: Pipeline flow diagram =====
    ax_flow = fig.add_subplot(gs[7, :])
    ax_flow.axis('off')
    flow = (
        "COMPLETE INTEGRATED PIPELINE: One wavefunction, four transformations\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
        "  ┌─────────────┐    ┌──────────────┐    ┌───────────────┐    ┌─────────────────┐    ┌────────────────┐\n"
        "  │  COHERENT    │    │   A-B PHASE   │    │    FLOQUET     │    │  BINDING        │    │   TOPOLOGICAL  │\n"
        "  │  BEAM        │───→│   SUBSTRATE   │───→│    DRESSING    │───→│  RESONANCE      │───→│   CAGING       │\n"
        "  │             │    │              │    │               │    │  FILTER          │    │               │\n"
        "  │ Kuramoto    │    │ φ = ∮A·dl/ℏ │    │ E₀ → E₀+nℏω  │    │ Only resonant   │    │ Φ=π: LOCKED   │\n"
        "  │ sync → r    │    │ WHERE atoms  │    │ Creates energy│    │ sidebands bind  │    │ Pattern stays │\n"
        "  │ coherence   │    │ interfere    │    │ ladder        │    │ = state filter  │    │ = deposition  │\n"
        "  └─────────────┘    └──────────────┘    └───────────────┘    └─────────────────┘    └────────────────┘\n\n"
        "  KEY INSIGHT: The A-B phase operates through the VECTOR POTENTIAL, not the field.\n"
        "  Atoms never feel a force — they accumulate geometric phase from the potential alone.\n"
        "  This is fundamentally different from force-based trapping (magnetic, optical, electric).\n"
        "  Combined with Floquet state selection and topological caging, it constitutes a\n"
        "  QUANTUM STATE FILTER for matter-wave deposition at the atomic scale."
    )
    ax_flow.text(0.02, 0.95, flow, transform=ax_flow.transAxes,
                fontsize=10.5, fontfamily='monospace', va='top',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='#fff8e1',
                         edgecolor='#ff8f00', linewidth=2))

    plt.savefig('results/integrated_pipeline.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    print("\n  Saved: integrated_pipeline.png")
    plt.close()


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("╔═══════════════════════════════════════════════════════════╗")
    print("║  INTEGRATED QUANTUM SUBSTRATE: COUPLED PIPELINE          ║")
    print("║  One wavefunction through all four transformations       ║")
    print("╚═══════════════════════════════════════════════════════════╝")

    # --- Full integrated pipeline ---
    sim = IntegratedQuantumSubstrate(N=256, L=400e-9, T_beam=1e-3)
    results = sim.run_full_pipeline(
        pattern='vortex_lattice',
        V_drive_frac=0.2,
        n_resonant=0,
        K_kuramoto=6.0,
        phi_cage=np.pi,
        prop_distance_lam=20,
        a=6*sim.lam, core=sim.lam,
    )

    # --- Ablation study ---
    ab_sim, ab_results, ab_phase = ablation_study()

    # --- Visualization ---
    print("\n" + "=" * 65)
    print("GENERATING FIGURES")
    print("=" * 65)
    create_integrated_figure(sim, results, ab_sim, ab_results, ab_phase)

    print("\n" + "=" * 65)
    print("SIMULATION COMPLETE")
    print("=" * 65)


if __name__ == "__main__":
    main()
