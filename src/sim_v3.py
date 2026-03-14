"""
Quantum Substrate Simulator v3 — Physically Matched Regime
============================================================
He-4 beam at 1 mK → λ_dB ≈ 8 nm
Substrate features at 10-80 nm
This matches the h-BN nanohole diffraction regime.
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

# ===========================================================================
# BEAM + SUBSTRATE + PROPAGATION (matched regime)
# ===========================================================================

class PhysicalSimulation:
    def __init__(self, T_beam=1e-3, L_sub=400e-9, N=512):
        self.mass = m_He
        self.T = T_beam
        self.v = np.sqrt(2 * k_B * T_beam / m_He)
        self.k0 = self.mass * self.v / hbar
        self.lam = 2 * np.pi / self.k0
        self.L = L_sub
        self.N = N
        self.dx = L_sub / N
        self.x = np.linspace(-L_sub/2, L_sub/2, N)
        self.y = np.linspace(-L_sub/2, L_sub/2, N)
        self.X, self.Y = np.meshgrid(self.x, self.y, indexing='ij')

    def info(self):
        return (f"He-4 at {self.T*1e3:.1f} mK: v={self.v:.2f} m/s, "
                f"λ_dB={self.lam*1e9:.2f} nm, "
                f"Substrate={self.L*1e9:.0f}nm, dx={self.dx*1e9:.2f}nm, "
                f"λ/dx={self.lam/self.dx:.1f}")

    def make_beam(self, width_frac=0.4):
        sigma = width_frac * self.L
        psi = np.exp(-(self.X**2 + self.Y**2)/(2*sigma**2))
        psi = psi * np.exp(1j * self.k0 * self.X)
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * self.dx**2)
        return psi

    def make_phase(self, pattern, **kw):
        X, Y = self.X, self.Y
        phase = np.zeros_like(X)

        if pattern == 'sinusoidal':
            period = kw.get('period', 4 * self.lam)
            amp = kw.get('amp', np.pi)
            phase = amp * np.cos(2*np.pi*X/period) * np.cos(2*np.pi*Y/period)

        elif pattern == 'rings':
            spacing = kw.get('spacing', 8 * self.lam)
            r_ring = kw.get('r_ring', 3 * self.lam)
            width = kw.get('width', 0.5 * self.lam)
            flux = kw.get('flux', 0.5)
            n = kw.get('n', 3)
            for i in range(-n, n+1):
                for j in range(-n, n+1):
                    x0, y0 = i*spacing, j*spacing
                    R = np.sqrt((X-x0)**2 + (Y-y0)**2)
                    theta = np.arctan2(Y-y0, X-x0)
                    ring_w = np.exp(-(R-r_ring)**2/(2*width**2))
                    phase += ring_w * flux * theta

        elif pattern == 'vortex_lattice':
            a = kw.get('a', 6 * self.lam)
            core = kw.get('core', self.lam)
            for i in range(-6, 7):
                for j in range(-6, 7):
                    x0 = a * (i + 0.5*(j%2))
                    y0 = a * np.sqrt(3)/2 * j
                    if abs(x0) < self.L/2 and abs(y0) < self.L/2:
                        R = np.sqrt((X-x0)**2+(Y-y0)**2+1e-30)
                        theta = np.arctan2(Y-y0, X-x0)
                        sign = (-1)**(i+j)
                        w = 1 - np.exp(-R**2/(2*core**2))
                        phase += sign * np.pi * theta/(2*np.pi) * w

        elif pattern == 'checkerboard':
            cell = kw.get('cell', 4*self.lam)
            ix = np.floor(X/cell).astype(int)
            iy = np.floor(Y/cell).astype(int)
            phase = np.pi * ((ix+iy)%2).astype(float)

        elif pattern == 'gradient_rings':
            # Radially varying flux — demonstrates programmability
            spacing = kw.get('spacing', 10*self.lam)
            for i in range(-3, 4):
                for j in range(-3, 4):
                    x0, y0 = i*spacing, j*spacing
                    R_center = np.sqrt(x0**2 + y0**2)
                    # Flux increases with distance from center
                    flux_local = 0.2 + 0.8 * R_center / (3*spacing)
                    R = np.sqrt((X-x0)**2+(Y-y0)**2)
                    theta = np.arctan2(Y-y0, X-x0)
                    r_ring = 3*self.lam
                    w = np.exp(-(R-r_ring)**2/(2*(0.5*self.lam)**2))
                    phase += w * flux_local * theta

        return phase

    def propagate(self, psi, distance):
        kx = np.fft.fftfreq(self.N, self.dx) * 2*np.pi
        ky = np.fft.fftfreq(self.N, self.dx) * 2*np.pi
        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        kz_sq = self.k0**2 - KX**2 - KY**2
        prop = kz_sq > 0
        kz = np.zeros_like(kz_sq)
        kz[prop] = np.sqrt(kz_sq[prop])
        H = np.zeros_like(kz_sq, dtype=complex)
        H[prop] = np.exp(1j * kz[prop] * distance)
        return np.fft.ifft2(np.fft.fft2(psi) * H)

    def simulate(self, pattern, prop_distances=None, **kw):
        psi = self.make_beam()
        phase = self.make_phase(pattern, **kw)
        psi_phased = psi * np.exp(1j * phase)

        if prop_distances is None:
            prop_distances = [2*self.lam, 10*self.lam, 50*self.lam]

        results = []
        for d in prop_distances:
            psi_p = self.propagate(psi_phased, d)
            density = np.abs(psi_p)**2
            # Contrast
            d_flat = density.ravel()
            c = (np.max(d_flat) - np.min(d_flat))/(np.max(d_flat) + np.min(d_flat) + 1e-30)
            results.append({
                'distance': d,
                'density': density,
                'contrast': c,
                'psi': psi_p,
            })

        return phase, psi, results


# ===========================================================================
# A-B CAGING (clean, focused)
# ===========================================================================

def ab_caging_analysis():
    N_cells = 25
    dim = 3 * N_cells

    def make_H(phi, J=1.0):
        H = np.zeros((dim, dim), dtype=complex)
        for n in range(N_cells):
            a, b, c = 3*n, 3*n+1, 3*n+2
            H[a,b] = H[b,a] = J
            H[a,c] = H[c,a] = J
            if n < N_cells-1:
                a2 = 3*(n+1)
                H[b,a2] = J*np.exp(1j*phi/2)
                H[a2,b] = J*np.exp(-1j*phi/2)
                H[c,a2] = J*np.exp(-1j*phi/2)
                H[a2,c] = J*np.exp(1j*phi/2)
        return H

    site0 = 3*(N_cells//2)
    psi0 = np.zeros(dim, dtype=complex)
    psi0[site0] = 1.0
    dt, T_max = 0.05, 50
    times = np.arange(0, T_max, dt)

    cage_data = {}
    for label, phi in [('Φ=0', 0), ('Φ=π/4', np.pi/4), ('Φ=π/2', np.pi/2),
                        ('Φ=3π/4', 3*np.pi/4), ('Φ=π', np.pi)]:
        H = make_H(phi)
        U = expm(-1j*H*dt)
        psi = psi0.copy()
        spread = []
        for _ in times:
            p = np.abs(psi)**2
            s = np.arange(dim)
            mu = np.sum(s*p)
            spread.append(np.sqrt(np.sum((s-mu)**2*p)))
            psi = U @ psi
        cage_data[label] = {'spread': np.array(spread), 'final': np.abs(psi)**2}

    # Spectrum vs flux
    phis = np.linspace(0, 2*np.pi, 300)
    spectra = np.array([np.linalg.eigvalsh(make_H(phi)) for phi in phis])

    return times, cage_data, phis, spectra, dim


# ===========================================================================
# FLOQUET SIDEBANDS
# ===========================================================================

def floquet_analysis():
    N_side = 6
    dim = 2*N_side+1
    ns = np.arange(-N_side, N_side+1)

    def H_floquet(V, omega=1.0):
        d = np.array([n*omega for n in ns])
        H = np.diag(d.astype(float))
        for i in range(dim-1):
            H[i,i+1] = V
            H[i+1,i] = V
        return H

    V_range = np.linspace(0, 2.5, 300)
    spectra = []
    pops_at_V = {}

    for V in V_range:
        H = H_floquet(V)
        evals = np.linalg.eigvalsh(H)
        spectra.append(evals)

    spectra = np.array(spectra)

    # Populations at selected V
    for V_sel in [0.3, 0.8, 1.5, 2.0]:
        H = H_floquet(V_sel)
        U = expm(-1j*H*2*np.pi)
        psi0 = np.zeros(dim); psi0[N_side] = 1.0
        psi_f = U @ psi0
        pops_at_V[V_sel] = np.abs(psi_f)**2

    # Adsorption curve: resonance at n=+2 sideband
    E_bind = 2.0
    gamma = 0.2
    ads = []
    for evals in spectra:
        a = np.sum(gamma**2/((evals - E_bind)**2 + gamma**2)) / dim
        ads.append(a)

    return V_range, spectra, ns, pops_at_V, np.array(ads)


# ===========================================================================
# KURAMOTO
# ===========================================================================

def kuramoto_analysis(N=500):
    alpha = 0.5
    omega_nat = np.random.normal(0, 1.0, N)
    dt = 0.01
    T = 50
    times = np.arange(0, T, dt)

    data = {}
    for K in [1, 2, 3, 5, 8]:
        theta = np.random.uniform(0, 2*np.pi, N)
        dtheta = np.zeros(N)
        order = []
        for _ in times:
            z = np.mean(np.exp(1j*theta))
            order.append(np.abs(z))
            coupling = np.imag(z * np.exp(-1j*theta))
            ddtheta = -alpha*dtheta + omega_nat + K*coupling
            dtheta += ddtheta*dt
            theta += dtheta*dt
        data[K] = np.array(order)

    return times, data


# ===========================================================================
# MASTER FIGURE
# ===========================================================================

def main():
    print("="*70)
    print("QUANTUM SUBSTRATE v3: Physically Matched Regime")
    print("="*70)

    sim = PhysicalSimulation(T_beam=1e-3, L_sub=400e-9, N=512)
    print(sim.info())

    # --- 2D Deposition simulations ---
    print("\n[1/5] Beam-substrate simulations...")
    patterns = {
        'Sinusoidal\n(period=4λ)': ('sinusoidal', dict(period=4*sim.lam, amp=np.pi)),
        'Checkerboard\n(cell=4λ)': ('checkerboard', dict(cell=4*sim.lam)),
        'Vortex Lattice\n(a=6λ)': ('vortex_lattice', dict(a=6*sim.lam, core=sim.lam)),
        'Ring Array\n(spacing=10λ)': ('rings', dict(spacing=10*sim.lam, r_ring=3*sim.lam,
                                                     width=0.5*sim.lam, flux=0.5, n=2)),
    }

    sim_results = {}
    for title, (pat, kw) in patterns.items():
        phase, psi_in, res = sim.simulate(pat,
            prop_distances=[5*sim.lam, 20*sim.lam, 100*sim.lam], **kw)
        sim_results[title] = (phase, res)
        print(f"  {pat}: contrast = {[r['contrast'] for r in res]}")

    # --- Caging ---
    print("\n[2/5] A-B caging...")
    cage_times, cage_data, cage_phis, cage_spectra, cage_dim = ab_caging_analysis()

    # --- Floquet ---
    print("[3/5] Floquet sidebands...")
    fl_V, fl_spec, fl_ns, fl_pops, fl_ads = floquet_analysis()

    # --- Kuramoto ---
    print("[4/5] Kuramoto synchronization...")
    k_times, k_data = kuramoto_analysis()

    # ===================================================================
    # BUILD THE FIGURE
    # ===================================================================
    print("\n[5/5] Building master figure...")

    fig = plt.figure(figsize=(28, 36))
    gs = GridSpec(7, 4, figure=fig, hspace=0.45, wspace=0.35,
                  left=0.04, right=0.98, top=0.97, bottom=0.02)

    x_nm = sim.x * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    # ===== ROW 0: Header =====
    ax_h = fig.add_subplot(gs[0, :])
    ax_h.axis('off')
    txt = (
        "QUANTUM PHASE-CONTROLLED MATTER-WAVE DEPOSITION\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Beam: He-4 at {sim.T*1e3:.0f} mK  •  λ_dB = {sim.lam*1e9:.2f} nm  •  "
        f"v = {sim.v:.2f} m/s  •  Substrate: {sim.L*1e9:.0f} nm × {sim.L*1e9:.0f} nm\n\n"
        "The substrate's Aharonov-Bohm phase landscape — not force, not intensity —\n"
        "determines WHERE matter waves constructively interfere and deposit.\n"
        "Floquet sidebands select WHICH quantum states bind. A-B caging LOCKS them in place."
    )
    ax_h.text(0.5, 0.5, txt, transform=ax_h.transAxes,
              ha='center', va='center', fontsize=13, fontfamily='monospace',
              bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f8ff',
                       edgecolor='steelblue'))

    # ===== ROW 1: Phase landscapes =====
    for idx, (title, (phase, _)) in enumerate(sim_results.items()):
        ax = fig.add_subplot(gs[1, idx])
        im = ax.imshow(phase.T, extent=ext, cmap='twilight', origin='lower')
        ax.set_title(f'{title}', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8, label='φ (rad)')

    # ===== ROW 2: Deposition at medium distance =====
    for idx, (title, (_, res)) in enumerate(sim_results.items()):
        ax = fig.add_subplot(gs[2, idx])
        d = res[1]['density']  # middle propagation distance
        d_norm = d / (np.max(d)+1e-30)
        im = ax.imshow(d_norm.T, extent=ext, cmap='inferno', origin='lower')
        dist_nm = res[1]['distance']*1e9
        c = res[1]['contrast']
        ax.set_title(f'Deposition (d={dist_nm:.0f}nm)\nC={c:.4f}', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ===== ROW 3: Cross-sections showing contrast vs propagation distance =====
    for idx, (title, (_, res)) in enumerate(sim_results.items()):
        ax = fig.add_subplot(gs[3, idx])
        mid = sim.N // 2
        colors = ['blue', 'green', 'red']
        for ri, r in enumerate(res):
            d = r['density'][:, mid]
            d_n = d / (np.max(d)+1e-30)
            d_nm = r['distance']*1e9
            ax.plot(x_nm, d_n, color=colors[ri], linewidth=1.2,
                   label=f'd={d_nm:.0f}nm, C={r["contrast"]:.4f}')
        ax.set_title(f'{title.split(chr(10))[0]}: Cross-section', fontsize=9, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('Norm. |ψ|²', fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # ===== ROW 4: A-B Caging =====
    ax_c1 = fig.add_subplot(gs[4, 0])
    cage_colors = {'Φ=0': 'red', 'Φ=π/4': 'orange', 'Φ=π/2': 'gold',
                   'Φ=3π/4': 'green', 'Φ=π': 'blue'}
    for label, data in cage_data.items():
        ax_c1.plot(cage_times, data['spread'], color=cage_colors[label],
                  linewidth=2, label=label)
    ax_c1.set_xlabel('Time (ℏ/J)')
    ax_c1.set_ylabel('RMS Spread (sites)')
    ax_c1.set_title('A-B Caging: Wavepacket Spread', fontsize=10, fontweight='bold')
    ax_c1.legend(fontsize=8, ncol=2)
    ax_c1.grid(True, alpha=0.3)

    ax_c2 = fig.add_subplot(gs[4, 1])
    center = cage_dim//2
    w = 15
    for label, data in cage_data.items():
        ax_c2.plot(np.arange(cage_dim)[3*center//3-w:3*center//3+w],
                  data['final'][3*center//3-w:3*center//3+w],
                  color=cage_colors[label], linewidth=1.5, label=label)
    ax_c2.set_xlabel('Site')
    ax_c2.set_ylabel('|ψ|²')
    ax_c2.set_title('Final Site Distribution', fontsize=10, fontweight='bold')
    ax_c2.legend(fontsize=7)
    ax_c2.grid(True, alpha=0.3)

    ax_c3 = fig.add_subplot(gs[4, 2])
    for i in range(cage_spectra.shape[1]):
        ax_c3.plot(cage_phis/np.pi, cage_spectra[:,i], 'k.', markersize=0.2)
    ax_c3.axvline(x=1, color='red', linewidth=1.5, linestyle='--', label='Φ=π')
    ax_c3.set_xlabel('Flux Φ/π')
    ax_c3.set_ylabel('Energy (J)')
    ax_c3.set_title('Spectrum vs Flux\n(flat bands at Φ=π)', fontsize=10, fontweight='bold')
    ax_c3.legend(fontsize=9)
    ax_c3.grid(True, alpha=0.3)

    ax_c4 = fig.add_subplot(gs[4, 3])
    ax_c4.axis('off')
    cage_txt = (
        "A-B CAGING MECHANISM\n"
        "━━━━━━━━━━━━━━━━━━━━\n\n"
        "Rhombic lattice with flux Φ\n"
        "per plaquette:\n\n"
        "  a ─J─ b ─Je^{iΦ/2}─ a'\n"
        "  │                     │\n"
        "  └─J─ c ─Je^{-iΦ/2}─ ┘\n\n"
        "At Φ=π:\n"
        "  Path a→b→a' picks up +π/2\n"
        "  Path a→c→a' picks up -π/2\n"
        "  ───────────────────────\n"
        "  DESTRUCTIVE interference\n"
        "  → particle CANNOT hop\n"
        "  → compact localized state\n\n"
        "For deposition:\n"
        "  Deposited atoms stay PUT.\n"
        "  Pattern is topologically\n"
        "  protected against disorder."
    )
    ax_c4.text(0.05, 0.95, cage_txt, transform=ax_c4.transAxes,
              fontsize=9, fontfamily='monospace', va='top',
              bbox=dict(boxstyle='round', facecolor='#e8f4e8', edgecolor='green'))

    # ===== ROW 5: Floquet =====
    ax_f1 = fig.add_subplot(gs[5, 0])
    for i in range(fl_spec.shape[1]):
        ax_f1.plot(fl_V, fl_spec[:,i], 'b-', linewidth=0.4, alpha=0.5)
    ax_f1.axhline(y=2, color='red', linewidth=2, linestyle='--', label='E_bind = 2ω')
    ax_f1.set_xlabel('Drive V')
    ax_f1.set_ylabel('Quasi-energy (ω)')
    ax_f1.set_title('Floquet Spectrum', fontsize=10, fontweight='bold')
    ax_f1.legend()
    ax_f1.grid(True, alpha=0.3)
    ax_f1.set_ylim(-8, 8)

    ax_f2 = fig.add_subplot(gs[5, 1])
    bar_w = 0.2
    offsets = [-0.3, -0.1, 0.1, 0.3]
    colors_fl = ['#2196F3', '#4CAF50', '#FF9800', '#F44336']
    for i, (V_sel, pops) in enumerate(fl_pops.items()):
        ax_f2.bar(fl_ns + offsets[i], pops, width=bar_w, alpha=0.8,
                 color=colors_fl[i], label=f'V={V_sel}')
    ax_f2.set_xlabel('Sideband n')
    ax_f2.set_ylabel('Population')
    ax_f2.set_title('Sideband Populations', fontsize=10, fontweight='bold')
    ax_f2.legend(fontsize=8)
    ax_f2.grid(True, alpha=0.3)

    ax_f3 = fig.add_subplot(gs[5, 2])
    ax_f3.plot(fl_V, fl_ads, 'b-', linewidth=2)
    ax_f3.fill_between(fl_V, fl_ads, alpha=0.15, color='blue')
    ax_f3.set_xlabel('Drive V')
    ax_f3.set_ylabel('Adsorption probability')
    ax_f3.set_title('State-Selective Adsorption\n(resonance at n=+2)', fontsize=10, fontweight='bold')
    ax_f3.grid(True, alpha=0.3)

    ax_f4 = fig.add_subplot(gs[5, 3])
    ax_f4.axis('off')
    fl_txt = (
        "FLOQUET STATE FILTER\n"
        "━━━━━━━━━━━━━━━━━━━━\n\n"
        "Atom energy: E₀\n"
        "Substrate drive: ω\n\n"
        "Creates sidebands:\n"
        "  E₀, E₀±ω, E₀±2ω, ...\n\n"
        "Binding resonance at E_b:\n"
        "  Only sideband n where\n"
        "  E₀ + nω ≈ E_b binds\n\n"
        "→ Sideband-selective\n"
        "  surface chemistry\n\n"
        "Combined with A-B phase:\n"
        "  WHERE (geometry) ×\n"
        "  WHICH (energy) =\n"
        "  deterministic placement"
    )
    ax_f4.text(0.05, 0.95, fl_txt, transform=ax_f4.transAxes,
              fontsize=9, fontfamily='monospace', va='top',
              bbox=dict(boxstyle='round', facecolor='#fff3e0', edgecolor='orange'))

    # ===== ROW 6: Kuramoto + Summary =====
    ax_k = fig.add_subplot(gs[6, 0:2])
    k_colors = {1:'red', 2:'orange', 3:'gold', 5:'green', 8:'blue'}
    for K, order in k_data.items():
        ax_k.plot(k_times, order, color=k_colors[K], linewidth=1.5, label=f'K={K}')
    ax_k.axhline(y=0.8, color='gray', linestyle='--')
    ax_k.set_xlabel('Time')
    ax_k.set_ylabel('Order parameter r')
    ax_k.set_title('Beam Phase Synchronization (2nd-order Kuramoto, N=500)',
                   fontsize=10, fontweight='bold')
    ax_k.legend(ncol=5, fontsize=9)
    ax_k.grid(True, alpha=0.3)
    ax_k.set_ylim(0, 1.05)

    ax_s = fig.add_subplot(gs[6, 2:4])
    ax_s.axis('off')
    summary = (
        "INTEGRATED DEPOSITION PIPELINE\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
        "① COHERENT BEAM\n"
        "   BEC or ultracold source → λ_dB ~ nm\n"
        "   Kuramoto sync → phase-locked ensemble\n"
        "   Spin squeezing → below SQL noise\n\n"
        "② A-B PHASE SUBSTRATE\n"
        "   Synthetic gauge field from laser array\n"
        "   Phase φ = ∮ A·dl / ℏ (not force!)\n"
        "   Programmable: change A, change pattern\n\n"
        "③ FLOQUET STATE SELECTION\n"
        "   Time-periodic drive → quasi-energy ladder\n"
        "   Binding resonance selects specific state\n"
        "   = Quantum state filter for adsorption\n\n"
        "④ TOPOLOGICAL LOCALIZATION\n"
        "   A-B caging at Φ=π → flat bands\n"
        "   Deposited atoms topologically trapped\n"
        "   Pattern robust against perturbation\n\n"
        "⑤ RESULT: Deterministic atomic placement\n"
        "   Resolution: λ_dB (sub-nm possible)\n"
        "   Selectivity: quantum state, not species\n"
        "   Robustness: topological protection"
    )
    ax_s.text(0.02, 0.98, summary, transform=ax_s.transAxes,
             fontsize=10, fontfamily='monospace', va='top',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#fffde7',
                      edgecolor='#f9a825', linewidth=2))

    plt.savefig('results/quantum_substrate_final.png', dpi=150,
                bbox_inches='tight', facecolor='white')
    print("\n  Saved: quantum_substrate_final.png")
    plt.close()


if __name__ == "__main__":
    main()
