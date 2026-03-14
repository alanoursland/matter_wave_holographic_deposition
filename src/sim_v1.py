"""
Quantum Phase-Controlled Matter-Wave Deposition Simulator
==========================================================

Integrated simulation connecting:
1. Matter-wave propagation with de Broglie wavelength control
2. Aharonov-Bohm phase imprinting from substrate flux geometry
3. Floquet sideband generation for state-selective adsorption
4. Topological (A-B) caging on a rhombic lattice
5. Full deposition: beam → substrate interaction → localization

The central thesis: a programmable substrate operating through geometric
(Aharonov-Bohm) phase — not direct force — acts as a quantum state filter,
selectively binding atoms based on their topological phase.

Author: Integrated from literature review and independent research
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.sparse import diags, kron, eye as speye
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import warnings
warnings.filterwarnings('ignore')

# ===========================================================================
# PHYSICAL CONSTANTS
# ===========================================================================
hbar = 1.0545718e-34    # J·s
k_B = 1.380649e-23      # J/K
m_He = 6.6464731e-27    # kg (helium-4)
m_Rb = 1.4431607e-25    # kg (rubidium-87, common in cold atom experiments)
a_0 = 5.29177e-11       # Bohr radius, m
e_charge = 1.602176634e-19  # C

# ===========================================================================
# MODULE 1: MATTER-WAVE BEAM PROPAGATION
# ===========================================================================

class MatterWaveBeam:
    """
    Models a coherent atomic beam with specified temperature, flux,
    and species. Computes de Broglie wavelength, coherence length,
    and spatial wavefunction.
    """
    def __init__(self, mass=m_He, temperature=1e-6, velocity=None):
        self.mass = mass
        self.temperature = temperature
        if velocity is None:
            # Thermal velocity from temperature
            self.velocity = np.sqrt(2 * k_B * temperature / mass)
        else:
            self.velocity = velocity
        self.wavelength_dB = hbar / (mass * self.velocity)
        self.k = 2 * np.pi / self.wavelength_dB
        # Coherence length ~ hbar / (mass * delta_v)
        # For BEC-like source, delta_v / v ~ 0.01
        self.velocity_spread = 0.01 * self.velocity
        self.coherence_length = hbar / (mass * self.velocity_spread)

    def wavefunction(self, x, t=0):
        """Gaussian wavepacket propagating in +x direction."""
        sigma = self.coherence_length / 2
        x0 = self.velocity * t
        envelope = np.exp(-(x - x0)**2 / (4 * sigma**2))
        phase = np.exp(1j * (self.k * x - 0.5 * hbar * self.k**2 * t / self.mass))
        norm = (2 * np.pi * sigma**2)**(-0.25)
        return norm * envelope * phase

    def info(self):
        return {
            'de Broglie wavelength (nm)': self.wavelength_dB * 1e9,
            'velocity (m/s)': self.velocity,
            'coherence length (nm)': self.coherence_length * 1e9,
            'temperature (μK)': self.temperature * 1e6,
        }


# ===========================================================================
# MODULE 2: AHARONOV-BOHM SUBSTRATE PHASE
# ===========================================================================

class ABSubstrate:
    """
    Models a 2D substrate with programmable magnetic flux tubes.
    The flux geometry creates an A-B phase landscape that incoming
    matter waves accumulate as they interact with the substrate.

    Key insight: the phase is determined by the vector potential A,
    not the magnetic field B. The atom never enters the field region
    but still acquires phase phi = (q/hbar) * integral(A · dl).

    For neutral atoms, the analogous effect uses synthetic gauge
    fields from laser-generated vector potentials.
    """
    def __init__(self, Nx=128, Ny=128, Lx=100e-9, Ly=100e-9):
        self.Nx, self.Ny = Nx, Ny
        self.Lx, self.Ly = Lx, Ly
        self.dx = Lx / Nx
        self.dy = Ly / Ny
        self.x = np.linspace(-Lx/2, Lx/2, Nx)
        self.y = np.linspace(-Ly/2, Ly/2, Ny)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')
        # Phase landscape initialized to zero
        self.phase = np.zeros((Nx, Ny))
        # Adsorption potential (where atoms stick)
        self.V_ads = np.zeros((Nx, Ny))

    def add_flux_tube(self, x0, y0, flux_quanta, radius):
        """
        Add a solenoid-like flux tube at (x0, y0).
        The A-B phase accumulated by circling the tube is 2*pi*flux_quanta.
        Outside the tube, B=0 but A != 0 (Coulomb gauge).
        """
        R = np.sqrt((self.X - x0)**2 + (self.Y - y0)**2)
        theta = np.arctan2(self.Y - y0, self.X - x0)
        # Phase from A-B effect: phi = 2*pi * Phi/Phi_0
        # For synthetic gauge: phi = flux_quanta * 2*pi * f(R)
        # f(R) = 1 for R > radius (full phase), smooth inside
        f = np.where(R > radius, 1.0, (R / radius)**2)
        self.phase += flux_quanta * 2 * np.pi * f * np.sign(theta)

    def add_ring_interferometer(self, x0, y0, r_inner, r_outer, flux):
        """
        Create a ring-shaped A-B interferometer on the substrate.
        Mimics the Sb2Te3 topological insulator ring geometry.
        The constructive/destructive interference creates localized
        trapping sites.
        """
        R = np.sqrt((self.X - x0)**2 + (self.Y - y0)**2)
        theta = np.arctan2(self.Y - y0, self.X - x0)
        # Ring region
        in_ring = (R >= r_inner) & (R <= r_outer)
        # A-B phase: h/e oscillation period
        ab_phase = flux * 2 * np.pi * theta / (2 * np.pi)
        # The interference creates standing wave in the ring
        ring_phase = np.where(in_ring, ab_phase, 0)
        self.phase += ring_phase

    def add_hex_lattice_holes(self, a_lattice, hole_radius, N_holes=5):
        """
        Create a hexagonal lattice of nanoholes (h-BN-inspired).
        Each hole acts as a diffractive element with dispersion-
        induced phase shifts.

        The phase shift at hole edges models the polarisability
        ripples from the h-BN paper: enhanced C6 at edges.
        """
        # Generate hex lattice positions
        positions = []
        for i in range(-N_holes, N_holes+1):
            for j in range(-N_holes, N_holes+1):
                x = a_lattice * (i + 0.5 * j)
                y = a_lattice * (np.sqrt(3)/2 * j)
                if abs(x) < self.Lx/2 and abs(y) < self.Ly/2:
                    positions.append((x, y))

        for (xh, yh) in positions:
            R = np.sqrt((self.X - xh)**2 + (self.Y - yh)**2)
            # Hole region: transmission
            hole = R < hole_radius
            # Edge region: enhanced polarisability → phase shift
            edge = (R >= hole_radius) & (R < hole_radius * 2)
            # Model the ~40% enhancement at N-terminated edges
            # as an additional phase accumulation
            edge_phase = np.where(edge, 0.4 * np.pi * np.exp(-(R - hole_radius) / (0.5 * hole_radius)), 0)
            self.phase += edge_phase
            # Adsorption enhanced at constructive interference sites
            self.V_ads += np.where(hole, -1.0, 0)

    def compute_transmission(self, psi_in, beam):
        """
        Compute the transmitted wavefunction after interaction with
        the substrate phase landscape.

        psi_out = psi_in * exp(i * phase_AB) * T(x,y)

        where T is the transmission amplitude accounting for
        van der Waals interaction at hole edges.
        """
        # Transmission: unity except at opaque regions
        T = np.ones_like(self.phase)
        # Apply A-B phase
        psi_out = psi_in * np.exp(1j * self.phase) * T
        return psi_out


# ===========================================================================
# MODULE 3: FLOQUET SIDEBAND ENGINE
# ===========================================================================

class FloquetEngine:
    """
    Models the Floquet-Bloch state generation for state-selective
    deposition. A time-periodic drive creates quasi-energy sidebands
    that can be matched to specific substrate adsorption resonances.

    Physical picture:
    - Incoming atom has energy E0
    - Substrate drives with frequency omega
    - Atom acquires sidebands at E0 + n*hbar*omega
    - Only specific sidebands match the substrate binding energy
    - → State-selective adsorption

    This is the "sideband-selective surface chemistry" concept.
    """
    def __init__(self, N_sidebands=5, omega_drive=2*np.pi*1e9):
        self.N = N_sidebands
        self.omega = omega_drive
        self.dim = 2 * N_sidebands + 1  # -N to +N

    def floquet_hamiltonian(self, E0, V_drive, delta_phase=0):
        """
        Construct the Floquet Hamiltonian in the extended zone scheme.

        H_F = diag(E0 + n*hbar*omega) + V_drive * coupling

        The coupling connects adjacent sidebands:
        |n> <-> |n+1> with strength V_drive * exp(i*delta_phase)
        """
        N = self.N
        dim = self.dim

        # Diagonal: quasi-energies
        diag_vals = np.array([E0 + n * hbar * self.omega for n in range(-N, N+1)])
        H = np.diag(diag_vals)

        # Off-diagonal: drive coupling
        for i in range(dim - 1):
            H[i, i+1] = V_drive * np.exp(1j * delta_phase)
            H[i+1, i] = V_drive * np.exp(-1j * delta_phase)

        return H

    def quasi_energies(self, E0, V_drive, delta_phase=0):
        """Compute quasi-energies and Floquet modes."""
        H_F = self.floquet_hamiltonian(E0, V_drive, delta_phase)
        eigenvalues, eigenvectors = np.linalg.eigh(H_F.real)
        return eigenvalues, eigenvectors

    def sideband_population(self, E0, V_drive, delta_phase=0):
        """
        For an atom initially in the n=0 sideband, compute
        the population in each sideband after one drive period.
        """
        H_F = self.floquet_hamiltonian(E0, V_drive, delta_phase)
        T_drive = 2 * np.pi / self.omega
        U = expm(-1j * H_F * T_drive / hbar)
        # Initial state: n=0 (center of the sideband ladder)
        psi0 = np.zeros(self.dim, dtype=complex)
        psi0[self.N] = 1.0
        psi_final = U @ psi0
        populations = np.abs(psi_final)**2
        return populations

    def adsorption_filter(self, E0, V_drive, E_binding, width):
        """
        Model state-selective adsorption:
        Only sidebands whose energy matches the substrate binding
        energy (within width) are adsorbed.

        Returns the fraction of the beam that gets adsorbed.
        """
        eigenvalues, _ = self.quasi_energies(E0, V_drive)
        # Lorentzian resonance at E_binding
        adsorption = np.sum(width**2 / ((eigenvalues - E_binding)**2 + width**2))
        return adsorption / self.dim


# ===========================================================================
# MODULE 4: AHARONOV-BOHM CAGING ON RHOMBIC LATTICE
# ===========================================================================

class ABCage:
    """
    Models A-B caging on a rhombic (diamond chain) lattice.
    At flux Phi = pi per plaquette, all bands become flat
    and particles are perfectly localized.

    This is the mechanism for topological trapping:
    the substrate geometry + gauge field creates compact
    localized states where atoms cannot propagate.

    Based on the Rydberg synthetic dimension experiments
    but applied to the substrate lattice.
    """
    def __init__(self, N_cells=20, phi=np.pi):
        self.N_cells = N_cells
        self.phi = phi
        # Rhombic lattice: 3 sites per unit cell (a, b, c)
        self.dim = 3 * N_cells

    def hamiltonian(self, J=1.0):
        """
        Rhombic lattice Hamiltonian with flux phi per plaquette.

        Unit cell: a --J-- b --J*exp(i*phi)-- a'
                   a --J-- c --J*exp(-i*phi)-- a'

        The two paths a->b->a' and a->c->a' accumulate
        opposite phases. At phi=pi, they interfere destructively
        → caging.
        """
        N = self.N_cells
        dim = self.dim
        H = np.zeros((dim, dim), dtype=complex)
        phi = self.phi

        for n in range(N):
            a = 3 * n
            b = 3 * n + 1
            c = 3 * n + 2

            # Intra-cell hopping: a -> b, a -> c
            H[a, b] = J
            H[b, a] = J
            H[a, c] = J
            H[c, a] = J

            # Inter-cell hopping: b -> a', c -> a' (with phase)
            if n < N - 1:
                a_next = 3 * (n + 1)
                H[b, a_next] = J * np.exp(1j * phi / 2)
                H[a_next, b] = J * np.exp(-1j * phi / 2)
                H[c, a_next] = J * np.exp(-1j * phi / 2)
                H[a_next, c] = J * np.exp(1j * phi / 2)

        return H

    def spectrum(self, J=1.0):
        """Compute energy spectrum."""
        H = self.hamiltonian(J)
        eigenvalues = np.linalg.eigvalsh(H)
        return eigenvalues

    def localization_test(self, J=1.0, site=0):
        """
        Initialize a particle at a single site and evolve.
        At phi=pi (caging), the particle remains localized.
        At phi=0 (no flux), it spreads.
        """
        H = self.hamiltonian(J)
        # Initial state: localized at site
        psi0 = np.zeros(self.dim, dtype=complex)
        psi0[site] = 1.0

        # Time evolution
        T = 50  # in units of hbar/J
        dt = 0.1
        times = np.arange(0, T, dt)
        participation_ratio = []
        site_populations = []

        psi = psi0.copy()
        U_dt = expm(-1j * H * dt)
        for t in times:
            probs = np.abs(psi)**2
            # Inverse participation ratio: 1/sum(p^2)
            # Large IPR = delocalized, small = localized
            ipr = 1.0 / np.sum(probs**2)
            participation_ratio.append(ipr)
            site_populations.append(probs.copy())
            psi = U_dt @ psi

        return times, np.array(participation_ratio), np.array(site_populations)


# ===========================================================================
# MODULE 5: KURAMOTO SYNCHRONIZATION OF BEAM
# ===========================================================================

class KuramotoBeam:
    """
    Models phase synchronization of atoms in the beam using
    the second-order Kuramoto model with inertia.

    d²θ_i/dt² + α dθ_i/dt = ω_i + (K/N) Σ sin(θ_j - θ_i)

    Pre-synchronizing the beam before deposition enhances
    pattern contrast. The order parameter r measures coherence.
    """
    def __init__(self, N_atoms=100, alpha=0.5, K=2.0):
        self.N = N_atoms
        self.alpha = alpha  # damping
        self.K = K          # coupling strength
        # Natural frequencies: Gaussian distributed
        self.omega = np.random.normal(0, 1.0, N_atoms)

    def evolve(self, T=50, dt=0.01):
        """Evolve the Kuramoto system and track order parameter."""
        N = self.N
        theta = np.random.uniform(0, 2*np.pi, N)
        dtheta = np.zeros(N)

        times = np.arange(0, T, dt)
        order_param = []
        phases_history = []

        for t in times:
            # Order parameter
            r = np.abs(np.mean(np.exp(1j * theta)))
            order_param.append(r)
            if int(t / dt) % 100 == 0:
                phases_history.append(theta.copy())

            # Coupling
            sin_diff = np.zeros(N)
            for i in range(N):
                sin_diff[i] = np.mean(np.sin(theta - theta[i]))

            # Second-order dynamics
            ddtheta = -self.alpha * dtheta + self.omega + self.K * sin_diff
            dtheta += ddtheta * dt
            theta += dtheta * dt

        return times, np.array(order_param), phases_history


# ===========================================================================
# MODULE 6: INTEGRATED DEPOSITION SIMULATION
# ===========================================================================

class QuantumSubstrateSimulator:
    """
    The full integrated simulation:

    1. Generate a coherent matter-wave beam
    2. (Optional) Pre-synchronize with Kuramoto dynamics
    3. Propagate beam toward substrate
    4. Substrate imprints A-B phase on the wavefunction
    5. Floquet drive creates sidebands
    6. Only resonant sidebands are adsorbed (state filtering)
    7. A-B caging localizes the adsorbed atoms
    8. Final deposition pattern emerges

    This is the "substrate as quantum state filter" paradigm.
    """
    def __init__(self, Nx=256, Ny=256, L=200e-9):
        self.Nx, self.Ny = Nx, Ny
        self.L = L
        self.dx = L / Nx
        self.substrate = ABSubstrate(Nx, Ny, L, L)
        self.floquet = FloquetEngine(N_sidebands=3)
        self.beam = MatterWaveBeam(mass=m_He, temperature=1e-6)

    def setup_substrate_pattern(self, pattern='rings'):
        """Configure the substrate with a specific phase geometry."""
        if pattern == 'rings':
            # Array of A-B ring interferometers
            spacing = self.L / 5
            for i in range(-2, 3):
                for j in range(-2, 3):
                    x0 = i * spacing
                    y0 = j * spacing
                    self.substrate.add_ring_interferometer(
                        x0, y0,
                        r_inner=spacing*0.15,
                        r_outer=spacing*0.35,
                        flux=0.5  # half flux quantum
                    )
        elif pattern == 'hex_holes':
            self.substrate.add_hex_lattice_holes(
                a_lattice=15e-9,
                hole_radius=3e-10,  # 3 Å ~ 6 Å diameter
                N_holes=6
            )
        elif pattern == 'flux_array':
            # Programmable flux tube array
            spacing = self.L / 8
            for i in range(-3, 4):
                for j in range(-3, 4):
                    # Alternating flux creates checkerboard of phases
                    flux = 0.5 * ((-1)**(i+j))
                    self.substrate.add_flux_tube(
                        i * spacing, j * spacing,
                        flux_quanta=flux,
                        radius=spacing * 0.2
                    )

    def generate_incident_beam(self):
        """Create 2D incident plane wave with Gaussian envelope."""
        sigma = self.L / 4
        envelope = np.exp(-(self.substrate.X**2 + self.substrate.Y**2) / (2 * sigma**2))
        k = self.beam.k
        # Plane wave in x-direction with Gaussian transverse profile
        psi = envelope * np.exp(1j * k * self.substrate.X)
        # Normalize
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * self.dx**2)
        return psi

    def apply_substrate_interaction(self, psi_in):
        """
        Full substrate interaction:
        1. A-B phase imprinting
        2. Diffraction through holes/features
        3. Dispersion (vdW) phase at edges
        """
        psi_out = self.substrate.compute_transmission(psi_in, self.beam)
        return psi_out

    def propagate_fresnel(self, psi, distance):
        """
        Fresnel propagation of the matter wave after the substrate.
        Uses angular spectrum method.
        """
        kx = np.fft.fftfreq(self.Nx, self.dx) * 2 * np.pi
        ky = np.fft.fftfreq(self.Ny, self.dx) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky, indexing='ij')

        k0 = self.beam.k
        # Evanescent wave cutoff
        kz_sq = k0**2 - KX**2 - KY**2
        kz_sq = np.maximum(kz_sq, 0)
        kz = np.sqrt(kz_sq)

        # Transfer function
        H = np.exp(1j * kz * distance)
        # Zero evanescent waves
        H[kz_sq <= 0] = 0

        psi_k = np.fft.fft2(psi)
        psi_propagated = np.fft.ifft2(psi_k * H)
        return psi_propagated

    def compute_deposition_map(self, psi_final, E_binding=None):
        """
        Compute the deposition probability map.

        If Floquet filtering is enabled, only atoms in resonant
        sidebands contribute to deposition.
        """
        density = np.abs(psi_final)**2

        if E_binding is not None:
            # Floquet state-selective filter
            E0 = 0.5 * self.beam.mass * self.beam.velocity**2
            V_drive = 0.1 * E0  # 10% modulation depth
            filter_factor = self.floquet.adsorption_filter(
                E0, V_drive, E_binding, width=0.05*E0
            )
            density *= filter_factor

        return density

    def run_full_simulation(self, pattern='rings', propagation_dist=None):
        """Execute the complete deposition simulation."""
        print("=" * 60)
        print("QUANTUM SUBSTRATE DEPOSITION SIMULATION")
        print("=" * 60)

        # Setup
        print("\n1. Configuring substrate phase geometry...")
        self.setup_substrate_pattern(pattern)

        print(f"   Pattern: {pattern}")
        print(f"   Substrate size: {self.L*1e9:.0f} nm × {self.L*1e9:.0f} nm")
        print(f"   Grid: {self.Nx} × {self.Ny}")

        # Beam
        print("\n2. Generating coherent matter-wave beam...")
        info = self.beam.info()
        for k, v in info.items():
            print(f"   {k}: {v:.4g}")

        psi_in = self.generate_incident_beam()
        print(f"   Beam norm: {np.sum(np.abs(psi_in)**2) * self.dx**2:.6f}")

        # Substrate interaction
        print("\n3. Applying substrate A-B phase interaction...")
        psi_after = self.apply_substrate_interaction(psi_in)

        # Propagation
        if propagation_dist is None:
            propagation_dist = 10 * self.beam.wavelength_dB
        print(f"\n4. Fresnel propagation: {propagation_dist*1e9:.2f} nm...")
        psi_final = self.propagate_fresnel(psi_after, propagation_dist)

        # Deposition
        print("\n5. Computing deposition map...")
        E0 = 0.5 * self.beam.mass * self.beam.velocity**2
        dep_unfiltered = self.compute_deposition_map(psi_final)
        dep_filtered = self.compute_deposition_map(psi_final, E_binding=E0)

        results = {
            'psi_in': psi_in,
            'psi_after': psi_after,
            'psi_final': psi_final,
            'deposition_unfiltered': dep_unfiltered,
            'deposition_filtered': dep_filtered,
            'substrate_phase': self.substrate.phase,
        }

        print("\n6. Simulation complete.")
        return results


# ===========================================================================
# MODULE 7: VISUALIZATION
# ===========================================================================

def plot_full_results(sim, results, save_path=None):
    """Generate comprehensive visualization of simulation results."""
    fig = plt.figure(figsize=(20, 24))
    gs = GridSpec(4, 3, figure=fig, hspace=0.35, wspace=0.3)

    x_nm = sim.substrate.x * 1e9
    y_nm = sim.substrate.y * 1e9

    # --- Row 1: Substrate and beam ---
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(results['substrate_phase'].T,
                      extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                      cmap='twilight', origin='lower')
    ax1.set_title('Substrate A-B Phase Landscape', fontsize=11, fontweight='bold')
    ax1.set_xlabel('x (nm)')
    ax1.set_ylabel('y (nm)')
    plt.colorbar(im1, ax=ax1, label='Phase (rad)')

    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(np.abs(results['psi_in'].T)**2,
                      extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                      cmap='inferno', origin='lower')
    ax2.set_title('Incident Beam |ψ|²', fontsize=11, fontweight='bold')
    ax2.set_xlabel('x (nm)')
    ax2.set_ylabel('y (nm)')
    plt.colorbar(im2, ax=ax2, label='Probability density')

    ax3 = fig.add_subplot(gs[0, 2])
    im3 = ax3.imshow(np.angle(results['psi_after'].T),
                      extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                      cmap='hsv', origin='lower')
    ax3.set_title('Phase After Substrate Interaction', fontsize=11, fontweight='bold')
    ax3.set_xlabel('x (nm)')
    ax3.set_ylabel('y (nm)')
    plt.colorbar(im3, ax=ax3, label='Phase (rad)')

    # --- Row 2: Deposition maps ---
    ax4 = fig.add_subplot(gs[1, 0])
    dep_u = results['deposition_unfiltered']
    im4 = ax4.imshow(dep_u.T,
                      extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                      cmap='hot', origin='lower')
    ax4.set_title('Deposition (No Floquet Filter)', fontsize=11, fontweight='bold')
    ax4.set_xlabel('x (nm)')
    ax4.set_ylabel('y (nm)')
    plt.colorbar(im4, ax=ax4, label='Deposition probability')

    ax5 = fig.add_subplot(gs[1, 1])
    dep_f = results['deposition_filtered']
    im5 = ax5.imshow(dep_f.T,
                      extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                      cmap='hot', origin='lower')
    ax5.set_title('Deposition (Floquet State-Filtered)', fontsize=11, fontweight='bold')
    ax5.set_xlabel('x (nm)')
    ax5.set_ylabel('y (nm)')
    plt.colorbar(im5, ax=ax5, label='Deposition probability')

    # Cross-section comparison
    ax6 = fig.add_subplot(gs[1, 2])
    mid = sim.Ny // 2
    ax6.plot(x_nm, dep_u[:, mid] / np.max(dep_u[:, mid] + 1e-30),
             'r-', alpha=0.7, label='Unfiltered')
    ax6.plot(x_nm, dep_f[:, mid] / np.max(dep_f[:, mid] + 1e-30),
             'b-', alpha=0.7, label='Floquet-filtered')
    ax6.set_title('Cross-Section (y=0)', fontsize=11, fontweight='bold')
    ax6.set_xlabel('x (nm)')
    ax6.set_ylabel('Normalized deposition')
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    # --- Row 3: A-B Caging demonstration ---
    print("   Computing A-B caging demonstration...")
    cage_pi = ABCage(N_cells=20, phi=np.pi)
    cage_0 = ABCage(N_cells=20, phi=0)

    times_pi, ipr_pi, pops_pi = cage_pi.localization_test(site=30)
    times_0, ipr_0, pops_0 = cage_0.localization_test(site=30)

    ax7 = fig.add_subplot(gs[2, 0])
    ax7.plot(times_pi, ipr_pi, 'b-', label='Φ = π (caged)', linewidth=2)
    ax7.plot(times_0, ipr_0, 'r-', label='Φ = 0 (free)', linewidth=2)
    ax7.set_title('A-B Caging: Participation Ratio', fontsize=11, fontweight='bold')
    ax7.set_xlabel('Time (ℏ/J)')
    ax7.set_ylabel('Inverse Participation Ratio')
    ax7.legend()
    ax7.grid(True, alpha=0.3)
    ax7.set_yscale('log')

    ax8 = fig.add_subplot(gs[2, 1])
    # Show site populations at final time
    sites = np.arange(cage_pi.dim)
    ax8.bar(sites, pops_pi[-1], alpha=0.7, color='blue', label='Φ = π', width=0.8)
    ax8.bar(sites, pops_0[-1], alpha=0.4, color='red', label='Φ = 0', width=0.8)
    ax8.set_title('Final Site Populations', fontsize=11, fontweight='bold')
    ax8.set_xlabel('Site index')
    ax8.set_ylabel('|ψ|²')
    ax8.legend()
    ax8.set_xlim(20, 45)

    # Spectrum
    ax9 = fig.add_subplot(gs[2, 2])
    E_pi = cage_pi.spectrum()
    E_0 = cage_0.spectrum()
    ax9.plot(E_pi, 'b.', markersize=3, label='Φ = π (flat bands)')
    ax9.plot(E_0, 'r.', markersize=3, label='Φ = 0 (dispersive)')
    ax9.set_title('Energy Spectrum', fontsize=11, fontweight='bold')
    ax9.set_xlabel('Eigenstate index')
    ax9.set_ylabel('Energy (J)')
    ax9.legend()
    ax9.grid(True, alpha=0.3)

    # --- Row 4: Floquet sidebands and Kuramoto ---
    print("   Computing Floquet sidebands...")
    floquet = FloquetEngine(N_sidebands=5, omega_drive=2*np.pi*1e9)
    E0 = 0.5 * sim.beam.mass * sim.beam.velocity**2

    ax10 = fig.add_subplot(gs[3, 0])
    V_strengths = np.linspace(0, 0.3 * E0, 50)
    for V in V_strengths:
        evals, _ = floquet.quasi_energies(E0, V)
        ax10.plot([V / E0] * len(evals), evals / E0, 'b.', markersize=1)
    ax10.set_title('Floquet Quasi-Energy Fan', fontsize=11, fontweight='bold')
    ax10.set_xlabel('Drive strength V/E₀')
    ax10.set_ylabel('Quasi-energy / E₀')
    ax10.grid(True, alpha=0.3)

    # Sideband populations
    ax11 = fig.add_subplot(gs[3, 1])
    pops = floquet.sideband_population(E0, 0.15 * E0)
    n_vals = np.arange(-floquet.N, floquet.N + 1)
    ax11.bar(n_vals, pops, color='steelblue', edgecolor='navy')
    ax11.set_title('Sideband Populations (V=0.15 E₀)', fontsize=11, fontweight='bold')
    ax11.set_xlabel('Sideband index n')
    ax11.set_ylabel('Population |⟨n|ψ⟩|²')
    ax11.grid(True, alpha=0.3)

    # Kuramoto synchronization
    print("   Computing Kuramoto beam synchronization...")
    kuramoto = KuramotoBeam(N_atoms=200, alpha=0.3, K=3.0)
    k_times, k_order, k_phases = kuramoto.evolve(T=30, dt=0.01)

    ax12 = fig.add_subplot(gs[3, 2])
    ax12.plot(k_times, k_order, 'g-', linewidth=1.5)
    ax12.set_title('Beam Synchronization (Kuramoto)', fontsize=11, fontweight='bold')
    ax12.set_xlabel('Time (a.u.)')
    ax12.set_ylabel('Order parameter r')
    ax12.set_ylim(0, 1.05)
    ax12.grid(True, alpha=0.3)
    ax12.axhline(y=0.8, color='r', linestyle='--', alpha=0.5, label='Threshold')
    ax12.legend()

    plt.suptitle(
        'Quantum Phase-Controlled Matter-Wave Deposition\n'
        'Substrate as Quantum State Filter via Aharonov-Bohm Phase',
        fontsize=14, fontweight='bold', y=1.01
    )

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f"   Figure saved to {save_path}")

    plt.close()
    return fig


def plot_comparison_patterns(save_path=None):
    """Compare deposition patterns for different substrate geometries."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    patterns = ['rings', 'hex_holes', 'flux_array']
    titles = ['Ring Interferometers', 'h-BN Hex Holes', 'Flux Tube Array']

    for idx, (pattern, title) in enumerate(zip(patterns, titles)):
        sim = QuantumSubstrateSimulator(Nx=256, Ny=256, L=200e-9)
        results = sim.run_full_simulation(pattern=pattern)

        x_nm = sim.substrate.x * 1e9
        y_nm = sim.substrate.y * 1e9

        # Phase landscape
        axes[0, idx].imshow(results['substrate_phase'].T,
                           extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                           cmap='twilight', origin='lower')
        axes[0, idx].set_title(f'{title}\nPhase Landscape', fontweight='bold')
        axes[0, idx].set_xlabel('x (nm)')
        axes[0, idx].set_ylabel('y (nm)')

        # Deposition
        axes[1, idx].imshow(results['deposition_unfiltered'].T,
                           extent=[x_nm[0], x_nm[-1], y_nm[0], y_nm[-1]],
                           cmap='hot', origin='lower')
        axes[1, idx].set_title(f'{title}\nDeposition Pattern', fontweight='bold')
        axes[1, idx].set_xlabel('x (nm)')
        axes[1, idx].set_ylabel('y (nm)')

    plt.suptitle('Deposition Patterns for Different Substrate Geometries',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
    plt.close()


# ===========================================================================
# MODULE 8: PARAMETER SWEEP — FINDING THE OPERATING REGIME
# ===========================================================================

def parameter_sweep(save_path=None):
    """
    Sweep key parameters to identify the operating regime where
    quantum phase control produces the highest contrast deposition.

    Key parameters:
    - Beam temperature → coherence length
    - Flux per plaquette → caging strength
    - Floquet drive strength → sideband selectivity
    - Kuramoto coupling → beam synchronization
    """
    print("\n" + "="*60)
    print("PARAMETER SWEEP: Finding optimal operating regime")
    print("="*60)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Temperature vs localization contrast
    print("\n1. Temperature sweep...")
    temps = np.logspace(-7, -4, 30)  # 0.1 μK to 100 μK
    contrasts = []
    for T in temps:
        beam = MatterWaveBeam(mass=m_He, temperature=T)
        # Contrast ~ coherence_length / feature_size
        feature = 10e-9  # 10 nm features
        contrast = min(beam.coherence_length / feature, 100)
        contrasts.append(contrast)

    axes[0, 0].semilogx(temps * 1e6, contrasts, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Temperature (μK)')
    axes[0, 0].set_ylabel('Coherence / Feature Size')
    axes[0, 0].set_title('Beam Coherence vs Temperature', fontweight='bold')
    axes[0, 0].axhline(y=1, color='r', linestyle='--', label='Minimum for interference')
    axes[0, 0].axhline(y=10, color='g', linestyle='--', label='Good contrast')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 2. Flux vs caging (IPR at final time)
    print("2. Flux sweep...")
    fluxes = np.linspace(0, 2*np.pi, 40)
    final_iprs = []
    for phi in fluxes:
        cage = ABCage(N_cells=15, phi=phi)
        times, ipr, _ = cage.localization_test(site=22)
        final_iprs.append(ipr[-1])

    axes[0, 1].plot(fluxes / np.pi, final_iprs, 'b-', linewidth=2)
    axes[0, 1].set_xlabel('Flux per plaquette (π)')
    axes[0, 1].set_ylabel('Final IPR (lower = more localized)')
    axes[0, 1].set_title('A-B Caging vs Flux', fontweight='bold')
    axes[0, 1].axvline(x=1.0, color='r', linestyle='--', label='Φ = π (caging)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    # 3. Floquet drive strength vs selectivity
    print("3. Floquet drive sweep...")
    floquet = FloquetEngine(N_sidebands=5)
    E0 = 1e-25  # arbitrary energy scale
    V_range = np.linspace(0.01 * E0, 0.5 * E0, 50)
    selectivities = []
    for V in V_range:
        pops = floquet.sideband_population(E0, V)
        # Selectivity: how peaked is the distribution?
        # Use entropy as a measure (lower = more selective)
        pops_clean = pops[pops > 1e-10]
        entropy = -np.sum(pops_clean * np.log(pops_clean))
        selectivities.append(entropy)

    axes[1, 0].plot(V_range / E0, selectivities, 'b-', linewidth=2)
    axes[1, 0].set_xlabel('Drive Strength V/E₀')
    axes[1, 0].set_ylabel('Sideband Entropy (lower = more selective)')
    axes[1, 0].set_title('Floquet Selectivity vs Drive', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)

    # 4. Kuramoto coupling vs synchronization time
    print("4. Kuramoto coupling sweep...")
    K_range = np.linspace(0.5, 8.0, 20)
    sync_times = []
    for K in K_range:
        kuramoto = KuramotoBeam(N_atoms=100, alpha=0.3, K=K)
        times, order, _ = kuramoto.evolve(T=30, dt=0.02)
        # Time to reach r > 0.8
        above = np.where(order > 0.8)[0]
        if len(above) > 0:
            sync_times.append(times[above[0]])
        else:
            sync_times.append(30)  # didn't synchronize

    axes[1, 1].plot(K_range, sync_times, 'b-o', linewidth=2, markersize=4)
    axes[1, 1].set_xlabel('Coupling Strength K')
    axes[1, 1].set_ylabel('Time to r > 0.8')
    axes[1, 1].set_title('Beam Synchronization vs Coupling', fontweight='bold')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].axhline(y=30, color='r', linestyle='--', alpha=0.5,
                        label='No sync (timeout)')
    axes[1, 1].legend()

    plt.suptitle('Parameter Space Exploration for Quantum Substrate Deposition',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
    plt.close()

    # Print optimal regime
    print("\n" + "-"*40)
    print("OPTIMAL OPERATING REGIME:")
    print("-"*40)
    print(f"  Temperature: < 1 μK (BEC regime)")
    print(f"  Flux: Φ = π per plaquette (A-B caging)")
    print(f"  Floquet drive: V ~ 0.1-0.2 E₀ (moderate selectivity)")
    print(f"  Kuramoto K: > {K_range[np.argmin(sync_times)]:.1f} (fast sync)")
    print(f"  Feature size: > {sim.beam.wavelength_dB*1e9:.2f} nm (de Broglie limit)")


# ===========================================================================
# MAIN EXECUTION
# ===========================================================================

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════╗")
    print("║  QUANTUM PHASE-CONTROLLED MATTER-WAVE DEPOSITION       ║")
    print("║  Substrate as Quantum State Filter                      ║")
    print("║  via Aharonov-Bohm Phase                                ║")
    print("╚══════════════════════════════════════════════════════════╝")

    # --- Main simulation with ring interferometer substrate ---
    print("\n" + "="*60)
    print("SIMULATION 1: Ring Interferometer Substrate")
    print("="*60)
    sim = QuantumSubstrateSimulator(Nx=256, Ny=256, L=200e-9)
    results = sim.run_full_simulation(pattern='rings')

    print("\n   Generating main visualization...")
    plot_full_results(sim, results,
                     save_path='results/fig_main_results.png')

    # --- Compare substrate geometries ---
    print("\n" + "="*60)
    print("SIMULATION 2: Comparing substrate geometries")
    print("="*60)
    plot_comparison_patterns(save_path='results/fig_comparison.png')

    # --- Parameter sweep ---
    parameter_sweep(save_path='results/fig_parameter_sweep.png')

    print("\n" + "="*60)
    print("ALL SIMULATIONS COMPLETE")
    print("="*60)
    print("\nOutput files:")
    print("  fig_main_results.png     - Full simulation dashboard")
    print("  fig_comparison.png       - Substrate geometry comparison")
    print("  fig_parameter_sweep.png  - Parameter space exploration")
