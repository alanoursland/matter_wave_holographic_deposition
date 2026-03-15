"""
Inverse Holography via Simulated SQUID Array
==============================================

Given a target 2D deposition pattern P_target(x,y), solve for the SQUID
array phase distribution that produces it via Aharonov-Bohm phase
modulation and Fresnel propagation of He-4 matter waves at 1 mK.

Two solvers:
  1. Gerchberg-Saxton  (alternating projections, phase-only constraint)
  2. Gradient descent   (PyTorch autograd through differentiable propagator)

Forward model:
  psi_in -> [SQUID phase screen] -> [Fresnel propagation] -> |psi_target|^2

Reuses propagator logic, beam generation, metrics, and physical constants
from sim_v9.
"""

import os
import sys
import time
import numpy as np
import torch
import torch.nn.functional as F
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

from sim_v9 import (
    IntegratedQuantumSubstrate,
    hbar, k_B, m_He,
    michelson_contrast, ssim_score, min_feature_size,
)

_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
mu_0 = 4 * np.pi * 1e-7   # vacuum permeability [H/m]
Phi_0 = 2.067833848e-15    # magnetic flux quantum [Wb]

os.makedirs('results', exist_ok=True)


# ===========================================================================
# SQUID ARRAY MODEL
# ===========================================================================

class SQUIDArray:
    """
    2D grid of N_loops x N_loops superconducting loops.

    Each loop has a controllable phase phi_i (the AB phase imprinted on a
    passing matter wave).  The array generates a continuous phase screen by
    bicubic interpolation from the discrete loop phases onto the simulation
    grid.
    """

    def __init__(self, N_loops=32, N_grid=256, L_grid=400e-9):
        self.N_loops = N_loops
        self.N_grid  = N_grid
        self.L_grid  = L_grid
        self.pitch   = L_grid / N_loops
        self.n_total = N_loops * N_loops

        # Loop center coordinates
        half = L_grid / 2
        pos  = np.linspace(-half + self.pitch / 2, half - self.pitch / 2,
                           N_loops)
        self.loop_x = pos
        self.loop_y = pos
        self.loop_X, self.loop_Y = np.meshgrid(pos, pos, indexing='ij')

    def phases_to_screen(self, phi_loops):
        """
        Interpolate discrete loop phases to continuous (N_grid, N_grid) screen.

        Parameters
        ----------
        phi_loops : torch.Tensor
            Shape (N_loops, N_loops) or (N_loops*N_loops,).

        Returns
        -------
        screen : torch.Tensor, shape (N_grid, N_grid), float64
        """
        if phi_loops.dim() == 1:
            phi_loops = phi_loops.reshape(self.N_loops, self.N_loops)
        # F.interpolate expects (B, C, H, W)
        inp = phi_loops.unsqueeze(0).unsqueeze(0).to(dtype=torch.float64)
        out = F.interpolate(inp, size=(self.N_grid, self.N_grid),
                            mode='bicubic', align_corners=True)
        return out.squeeze(0).squeeze(0)

    @staticmethod
    def phase_screen_to_transmission(screen):
        """T(x,y) = exp(i * screen).  Phase-only: |T| = 1."""
        return torch.exp(1j * screen.to(dtype=torch.complex128))

    # ------------------------------------------------------------------
    # MUTUAL INDUCTANCE MODEL
    # ------------------------------------------------------------------
    def build_inductance_matrix(self, wire_width_frac=0.1):
        """
        Build the mutual inductance matrix M for the loop array.

        Self-inductance of each square loop (side a, wire width w):
            L_self ≈ (2 mu_0 a / pi) * [ln(a/w) - 0.774]

        Mutual inductance between coplanar square loops (dipole approx):
            M_ij ≈ (mu_0 / (4 pi)) * a^4 / r_ij^3
            (valid for r >> a; nearest-neighbour coupling uses a numerical
            prefactor from Grover's tables.)

        Parameters
        ----------
        wire_width_frac : float
            Wire width as a fraction of loop side length (default 0.1).

        Returns
        -------
        M : np.ndarray, shape (n_total, n_total)
            Mutual inductance matrix in Henries.  M[i,j] = Phi_i / I_j.
        """
        a = self.pitch                          # loop side length
        w = wire_width_frac * a                  # wire width

        # Self-inductance per loop
        L_self = (2 * mu_0 * a / np.pi) * (np.log(a / w) - 0.774)

        # Flatten loop positions
        xs = self.loop_X.ravel()   # (n_total,)
        ys = self.loop_Y.ravel()

        # Pairwise distance matrix
        dx = xs[:, None] - xs[None, :]
        dy = ys[:, None] - ys[None, :]
        r = np.sqrt(dx**2 + dy**2)

        # Mutual inductance (dipole approximation, coplanar loops)
        with np.errstate(divide='ignore', invalid='ignore'):
            M = (mu_0 / (4 * np.pi)) * a**4 / r**3

        # Replace diagonal (self-inductance) and fix any inf/nan
        np.fill_diagonal(M, L_self)
        M = np.nan_to_num(M, nan=0.0, posinf=0.0, neginf=0.0)

        self._M = M
        self._M_inv = np.linalg.inv(M)
        self._L_self = L_self

        return M

    def flux_to_currents(self, phi_loops):
        """
        Convert AB phases to the drive currents needed to produce them.

        Phi_i = phi_i * Phi_0 / (2 pi)   (flux from phase)
        I = M^{-1} Phi                    (currents from flux)

        Parameters
        ----------
        phi_loops : np.ndarray, shape (N_loops, N_loops) or (n_total,)

        Returns
        -------
        currents : np.ndarray, shape (n_total,)  in Amperes.
        """
        if not hasattr(self, '_M_inv'):
            self.build_inductance_matrix()

        phi_flat = np.asarray(phi_loops).ravel()
        flux = phi_flat * Phi_0 / (2 * np.pi)
        return self._M_inv @ flux

    def currents_to_flux(self, currents):
        """
        Convert drive currents back to AB phases.

        Phi = M I   -> phi = 2 pi Phi / Phi_0

        Parameters
        ----------
        currents : np.ndarray, shape (n_total,)

        Returns
        -------
        phi_loops : np.ndarray, shape (N_loops, N_loops)
        """
        if not hasattr(self, '_M'):
            self.build_inductance_matrix()

        flux = self._M @ np.asarray(currents).ravel()
        phi = flux * (2 * np.pi) / Phi_0
        return phi.reshape(self.N_loops, self.N_loops)

    def current_map_summary(self, phi_loops):
        """
        Compute and return a dict summarising the current requirements
        for a given phase solution.
        """
        I = self.flux_to_currents(phi_loops)
        I_grid = I.reshape(self.N_loops, self.N_loops)

        # Roundtrip check: currents -> flux -> phase should recover input
        phi_rt = self.currents_to_flux(I)
        phi_in = np.asarray(phi_loops).ravel()
        rt_err = np.max(np.abs(phi_rt.ravel() - phi_in))

        return {
            'currents_A':   I,
            'currents_grid': I_grid,
            'I_max_uA':     float(np.max(np.abs(I)) * 1e6),
            'I_rms_uA':     float(np.sqrt(np.mean(I**2)) * 1e6),
            'roundtrip_err': float(rt_err),
        }

    def info(self):
        print(f"  SQUID array: {self.N_loops}x{self.N_loops} = "
              f"{self.n_total} loops")
        print(f"  Pitch: {self.pitch*1e9:.1f} nm")
        print(f"  Grid:  {self.N_grid}x{self.N_grid}, "
              f"L={self.L_grid*1e9:.0f} nm")
        if hasattr(self, '_L_self'):
            print(f"  L_self: {self._L_self*1e12:.3f} pH")
            cond = np.linalg.cond(self._M)
            print(f"  M matrix condition number: {cond:.1f}")


# ===========================================================================
# INVERSE HOLOGRAPHY SOLVER
# ===========================================================================

class InverseHolographySolver:
    """
    Solves the inverse holography problem.

    Given P_target(x,y), find the SQUID loop phases that make the
    forward-propagated intensity match P_target.
    """

    def __init__(self, squid_array, N=256, L=400e-9, T_beam=1e-3,
                 prop_distance_lam=20.0):
        self.squid  = squid_array
        self.N      = N
        self.L      = L
        self.dx     = L / N
        self.T_beam = T_beam
        self.prop_distance_lam = prop_distance_lam

        # Physical parameters (He-4 at T_beam)
        self.mass = m_He
        self.v    = np.sqrt(2 * k_B * T_beam / m_He)
        self.k0   = self.mass * self.v / hbar
        self.lam  = 2 * np.pi / self.k0
        self.z    = prop_distance_lam * self.lam

        # Generate input beam (plane-wave-like Gaussian for clean inverse problem)
        sigma = 0.35 * L
        x = np.linspace(-L / 2, L / 2, N)
        X, Y = np.meshgrid(x, x, indexing='ij')
        psi_in = np.exp(-(X**2 + Y**2) / (2 * sigma**2)).astype(np.complex128)
        psi_in /= np.sqrt(np.sum(np.abs(psi_in)**2) * self.dx**2)
        self.psi_in_np = psi_in
        self.A_in_np   = np.abs(psi_in)
        self.psi_in_t  = torch.tensor(psi_in, dtype=torch.complex128,
                                      device=_device)
        self.A_in_t    = torch.abs(self.psi_in_t)

        # Precompute propagation transfer function
        kx = torch.fft.fftfreq(N, self.dx, device=_device,
                                dtype=torch.float64) * 2 * np.pi
        ky = torch.fft.fftfreq(N, self.dx, device=_device,
                                dtype=torch.float64) * 2 * np.pi
        KX, KY = torch.meshgrid(kx, ky, indexing='ij')
        kz_sq = self.k0**2 - KX**2 - KY**2
        valid = kz_sq > 0
        kz = torch.where(valid, torch.sqrt(torch.clamp(kz_sq, min=0)),
                         torch.zeros_like(kz_sq))
        # Forward and backward transfer functions
        self._H_fwd = torch.where(
            valid, torch.exp(1j * kz * self.z),
            torch.zeros_like(kz, dtype=torch.complex128))
        self._H_bwd = torch.where(
            valid, torch.exp(-1j * kz * self.z),
            torch.zeros_like(kz, dtype=torch.complex128))

    def _propagate_torch(self, psi_t, forward=True):
        """Angular spectrum propagator — pure torch, differentiable."""
        H = self._H_fwd if forward else self._H_bwd
        return torch.fft.ifft2(torch.fft.fft2(psi_t) * H)

    def forward(self, phase_screen):
        """
        Full forward model: phase screen -> intensity at target plane.

        Parameters
        ----------
        phase_screen : torch.Tensor (N, N) float64

        Returns
        -------
        intensity : torch.Tensor (N, N) float64
        """
        T = SQUIDArray.phase_screen_to_transmission(phase_screen)
        psi_mod = self.psi_in_t * T
        psi_target = self._propagate_torch(psi_mod, forward=True)
        return torch.abs(psi_target)**2

    def forward_from_loops(self, phi_loops):
        """Forward model parameterised by discrete loop phases."""
        screen = self.squid.phases_to_screen(phi_loops)
        return self.forward(screen), screen

    # -------------------------------------------------------------------
    # GERCHBERG-SAXTON
    # -------------------------------------------------------------------
    def solve_gerchberg_saxton(self, P_target, n_iter=300, n_restarts=3,
                               verbose=True):
        if verbose:
            print(f"\n  GS solver: {n_iter} iter x {n_restarts} restarts")

        P_t = torch.tensor(P_target, dtype=torch.float64, device=_device)
        A_target = torch.sqrt(P_t / (P_t.sum() * self.dx**2 + 1e-30))

        best_nrmse  = np.inf
        best_screen = None
        best_I      = None
        best_hist   = None

        for restart in range(n_restarts):
            # Random initial phase at target plane
            phi_target = torch.rand(self.N, self.N, dtype=torch.float64,
                                    device=_device) * 2 * np.pi
            history = []

            for it in range(n_iter):
                # Back-propagate to SQUID plane
                psi_target = A_target * torch.exp(
                    1j * phi_target.to(dtype=torch.complex128))
                psi_grid = self._propagate_torch(psi_target, forward=False)

                # Phase-only constraint: keep phase, replace amplitude
                phi_grid = torch.angle(psi_grid)
                psi_grid_hat = self.A_in_t * torch.exp(
                    1j * phi_grid.to(dtype=torch.complex128))

                # Forward-propagate to target plane
                psi_target_new = self._propagate_torch(psi_grid_hat,
                                                       forward=True)

                # Amplitude constraint: keep phase, replace amplitude
                phi_target = torch.angle(psi_target_new).to(
                    dtype=torch.float64)

                # Compute error
                I_achieved = torch.abs(psi_target_new)**2
                I_norm = I_achieved / (I_achieved.sum() + 1e-30)
                P_norm = P_t / (P_t.sum() + 1e-30)
                nrmse = float(torch.sqrt(
                    torch.mean((I_norm - P_norm)**2)
                    / (torch.mean(P_norm**2) + 1e-30)).item())
                history.append(nrmse)

            if verbose:
                print(f"    restart {restart}: NRMSE={history[-1]:.6f}")

            if history[-1] < best_nrmse:
                best_nrmse  = history[-1]
                best_screen = phi_grid.detach().cpu().numpy()
                best_I      = I_achieved.detach().cpu().numpy()
                best_hist   = history

        if verbose:
            print(f"    best NRMSE: {best_nrmse:.6f}")

        return {
            'phase_screen': best_screen,
            'achieved':     best_I,
            'convergence':  best_hist,
            'nrmse':        best_nrmse,
        }

    # -------------------------------------------------------------------
    # GRADIENT DESCENT
    # -------------------------------------------------------------------
    def solve_gradient_descent(self, P_target, n_iter=500, lr=0.1,
                               reg_smooth=1e-4, verbose=True):
        if verbose:
            print(f"\n  GD solver: {n_iter} iter, lr={lr}, "
                  f"reg={reg_smooth}")

        P_t = torch.tensor(P_target, dtype=torch.float64, device=_device)
        P_t_norm = P_t / (P_t.max() + 1e-30)

        # Initialise with random phases (not zeros — flat phase gives no
        # gradient signal because the propagated beam is nearly symmetric)
        phi_loops = (torch.rand(self.squid.n_total, dtype=torch.float64,
                                device=_device) * 0.5
                     ).requires_grad_(True)
        optimizer = torch.optim.Adam([phi_loops], lr=lr)

        history = []
        best_loss = np.inf
        best_phi  = None

        for it in range(n_iter):
            optimizer.zero_grad()

            screen = self.squid.phases_to_screen(phi_loops)
            I = self.forward(screen)
            I_peak_norm = I / (I.max() + 1e-30)

            # Negative correlation loss: maximise overlap with target shape
            loss_corr = -torch.sum(I_peak_norm * P_t_norm) / (
                torch.sqrt(torch.sum(I_peak_norm**2) *
                           torch.sum(P_t_norm**2)) + 1e-30)

            # MSE on peak-normalised patterns for spatial accuracy
            loss_mse = torch.mean((I_peak_norm - P_t_norm)**2)

            # Total-variation smoothness penalty on phase screen
            grad_x = (screen[1:, :] - screen[:-1, :])**2
            grad_y = (screen[:, 1:] - screen[:, :-1])**2
            loss_smooth = reg_smooth * (grad_x.mean() + grad_y.mean())

            loss = loss_corr + loss_mse + loss_smooth
            loss.backward()
            optimizer.step()

            loss_val = float(loss.item())
            history.append(loss_val)

            if loss_val < best_loss:
                best_loss = loss_val
                best_phi  = phi_loops.detach().clone()

            if verbose and (it % 100 == 0 or it == n_iter - 1):
                print(f"    iter {it:4d}: loss={loss_val:.6e}  "
                      f"corr={loss_corr.item():.4f}  "
                      f"mse={loss_mse.item():.6e}")

        # Recover best solution
        phi_loops_best = best_phi
        screen_best = self.squid.phases_to_screen(phi_loops_best)
        I_best = self.forward(screen_best)

        # Wrap phases mod 2*pi
        phi_wrapped = (phi_loops_best % (2 * np.pi)).detach().cpu().numpy()

        return {
            'phi_loops':    phi_wrapped,
            'phase_screen': screen_best.detach().cpu().numpy(),
            'achieved':     I_best.detach().cpu().numpy(),
            'convergence':  history,
            'nrmse':        best_loss,
        }


# ===========================================================================
# TARGET PATTERN GENERATORS
# ===========================================================================

def target_single_spot(N, L, sigma_frac=0.05):
    """Gaussian spot at center."""
    x = np.linspace(-0.5, 0.5, N)
    X, Y = np.meshgrid(x, x, indexing='ij')
    P = np.exp(-(X**2 + Y**2) / (2 * sigma_frac**2))
    return P / (P.sum() + 1e-30)


def target_line(N, L, width_frac=0.02, angle_deg=0):
    """Line through center at given angle."""
    x = np.linspace(-0.5, 0.5, N)
    X, Y = np.meshgrid(x, x, indexing='ij')
    theta = np.radians(angle_deg)
    d = -X * np.sin(theta) + Y * np.cos(theta)
    P = np.exp(-d**2 / (2 * width_frac**2))
    return P / (P.sum() + 1e-30)


def target_grid_of_dots(N, L, n_dots=3, sigma_frac=0.03):
    """n_dots x n_dots grid of Gaussian spots."""
    x = np.linspace(-0.5, 0.5, N)
    X, Y = np.meshgrid(x, x, indexing='ij')
    P = np.zeros((N, N))
    positions = np.linspace(-0.3, 0.3, n_dots)
    for xi in positions:
        for yi in positions:
            P += np.exp(-((X - xi)**2 + (Y - yi)**2) / (2 * sigma_frac**2))
    return P / (P.sum() + 1e-30)


def target_ring(N, L, radius_frac=0.2, width_frac=0.03):
    """Annular ring pattern."""
    x = np.linspace(-0.5, 0.5, N)
    X, Y = np.meshgrid(x, x, indexing='ij')
    R = np.sqrt(X**2 + Y**2)
    P = np.exp(-(R - radius_frac)**2 / (2 * width_frac**2))
    return P / (P.sum() + 1e-30)


def target_letter(N, L, letter='H', line_frac=0.04):
    """Block letter rendered on grid."""
    P = np.zeros((N, N))
    w = max(int(N * line_frac), 1)
    margin = N // 5

    if letter == 'H':
        # Left vertical
        P[margin:N - margin, margin:margin + w] = 1.0
        # Right vertical
        P[margin:N - margin, N - margin - w:N - margin] = 1.0
        # Horizontal bar
        mid = N // 2
        P[mid - w // 2:mid + w // 2 + 1, margin:N - margin] = 1.0
    elif letter == 'L':
        P[margin:N - margin, margin:margin + w] = 1.0
        P[N - margin - w:N - margin, margin:N - margin] = 1.0
    elif letter == 'T':
        P[margin:margin + w, margin:N - margin] = 1.0
        mid = N // 2
        P[margin:N - margin, mid - w // 2:mid + w // 2 + 1] = 1.0
    else:
        # Fallback: cross
        mid = N // 2
        P[mid - w:mid + w, margin:N - margin] = 1.0
        P[margin:N - margin, mid - w:mid + w] = 1.0

    # Smooth slightly for physical plausibility
    P = gaussian_filter(P, sigma=2)
    return P / (P.sum() + 1e-30)


def bandlimit_target(P, N_loops, rolloff=0.8):
    """
    Low-pass filter a target pattern to the spatial bandwidth of the
    SQUID array.

    The phase screen is generated by bicubic interpolation from an
    N_loops x N_loops grid.  It cannot modulate at spatial frequencies
    above its Nyquist limit (N_loops / 2 cycles across the aperture).
    Any target content above that frequency is unreachable — the
    optimizer wastes energy chasing it, driving power into the highest
    modes it *can* produce, which then diffract out of the target region.

    This function zeroes out unreachable frequencies with a smooth
    Butterworth rolloff (no Gibbs ringing from a hard cutoff).

    Parameters
    ----------
    P : np.ndarray (N, N)
        Target pattern (non-negative).
    N_loops : int
        Number of SQUID loops per axis.
    rolloff : float
        Fraction of the Nyquist limit at which the filter reaches -3 dB.
        Default 0.8 (conservative — keeps only content the array can
        comfortably reproduce).

    Returns
    -------
    P_bl : np.ndarray (N, N), normalised.
    """
    N = P.shape[0]

    # Nyquist in grid-frequency units: N_loops/2 cycles per N pixels
    f_nyquist = N_loops / 2.0
    f_cutoff = rolloff * f_nyquist

    # Frequency grid (in cycles per N pixels)
    fx = np.fft.fftfreq(N) * N   # cycles across the grid
    fy = np.fft.fftfreq(N) * N
    FX, FY = np.meshgrid(fx, fy, indexing='ij')
    F_rad = np.sqrt(FX**2 + FY**2)

    # 4th-order Butterworth low-pass (smooth rolloff, no ringing)
    order = 4
    H = 1.0 / (1.0 + (F_rad / f_cutoff)**(2 * order))

    # Filter in Fourier domain
    P_fft = np.fft.fft2(P)
    P_filtered = np.real(np.fft.ifft2(P_fft * H))
    P_filtered = np.clip(P_filtered, 0, None)

    return P_filtered / (P_filtered.sum() + 1e-30)


def smooth_target(P, N_loops=None, corner_radius=None, sigma=None):
    """
    Apply physically motivated smoothing to a target pattern.

    Three complementary stages (all optional, applied in order):

    1. Band-limiting: zero spatial frequencies above the SQUID array's
       Nyquist limit.  This is the principled fix — you cannot reproduce
       what you cannot modulate.
    2. Corner rounding: morphological open+close with a disk element.
       Rounds sharp corners while preserving flat edges.
    3. Gaussian blur: gentle final smoothing.

    Parameters
    ----------
    P : np.ndarray (N, N)
        Raw target pattern (non-negative).
    N_loops : int or None
        SQUID loops per axis.  If provided, band-limits the target to the
        array's spatial bandwidth.
    corner_radius : float or None
        Corner rounding as a fraction of grid side (e.g. 0.03 = 3%).
    sigma : float or None
        Gaussian blur sigma in pixels.

    Returns
    -------
    P_smooth : np.ndarray (N, N), normalised, same shape as input.
    """
    from scipy.ndimage import binary_opening, binary_closing, distance_transform_edt

    P_out = P.copy()
    N = P.shape[0]

    # Stage 1: band-limit to SQUID Nyquist
    if N_loops is not None:
        P_out = bandlimit_target(P_out, N_loops)

    # Stage 2: morphological corner rounding
    if corner_radius is not None and corner_radius > 0:
        r_px = max(int(corner_radius * N), 1)

        # Build a disk structuring element
        y, x = np.ogrid[-r_px:r_px + 1, -r_px:r_px + 1]
        selem = (x**2 + y**2) <= r_px**2

        # Threshold to binary, apply morphological open+close (rounds
        # both convex and concave corners), then restore smooth profile
        thresh = 0.1 * P_out.max()
        mask = P_out > thresh

        mask = binary_opening(mask, structure=selem)
        mask = binary_closing(mask, structure=selem)

        # Feathered edge: use distance transform to create a smooth
        # transition instead of hard binary boundary
        dist = distance_transform_edt(mask)
        dist_inv = distance_transform_edt(~mask)
        # Signed distance (positive inside, negative outside)
        signed_dist = dist - dist_inv
        # Smooth sigmoid transition over ~r_px/2 pixels
        transition = 1.0 / (1.0 + np.exp(-signed_dist / max(r_px * 0.3, 1)))

        P_out = P_out * transition

    # Stage 3: Gaussian blur
    if sigma is not None and sigma > 0:
        P_out = gaussian_filter(P_out, sigma=sigma)

    P_out = np.clip(P_out, 0, None)
    return P_out / (P_out.sum() + 1e-30)


# ===========================================================================
# METRICS
# ===========================================================================

def compute_metrics(I_achieved, P_target):
    """Compute quality metrics between achieved and target intensity."""
    # Normalise
    I_n = I_achieved / (I_achieved.sum() + 1e-30)
    P_n = P_target / (P_target.sum() + 1e-30)

    nrmse = float(np.sqrt(np.mean((I_n - P_n)**2)
                          / (np.mean(P_n**2) + 1e-30)))
    ssim = ssim_score(P_target, I_achieved)
    contrast = michelson_contrast(I_achieved)

    # Diffraction efficiency: fraction of intensity in target region
    mask = P_target > 0.1 * P_target.max()
    eff = float(I_achieved[mask].sum() / (I_achieved.sum() + 1e-30))

    return {
        'nrmse':    nrmse,
        'ssim':     ssim,
        'contrast': contrast,
        'efficiency': eff,
    }


# ===========================================================================
# ROUNDTRIP VALIDATION
# ===========================================================================

def validate_roundtrip(solver, verbose=True):
    """
    Sanity check: apply a known phase, forward-propagate, then run the
    inverse solver.  The forward -> inverse -> forward chain should
    reproduce the intensity.
    """
    if verbose:
        print("\n  ROUNDTRIP VALIDATION")

    # Known sinusoidal phase screen
    x = np.linspace(-np.pi, np.pi, solver.N)
    X, Y = np.meshgrid(x, x, indexing='ij')
    phase_true = 0.8 * np.cos(4 * X) * np.cos(4 * Y)

    screen_t = torch.tensor(phase_true, dtype=torch.float64, device=_device)
    I_true = solver.forward(screen_t).detach().cpu().numpy()

    # Run GS to recover phase from intensity
    result = solver.solve_gerchberg_saxton(I_true, n_iter=200,
                                           n_restarts=2, verbose=False)

    # Forward-propagate recovered phase
    screen_rec = torch.tensor(result['phase_screen'], dtype=torch.float64,
                              device=_device)
    I_rec = solver.forward(screen_rec).detach().cpu().numpy()

    ssim_val = ssim_score(I_true, I_rec)
    ok = ssim_val > 0.5

    if verbose:
        print(f"    SSIM(I_true, I_recovered) = {ssim_val:.4f}  "
              f"{'PASS' if ok else 'FAIL'}")

    return ok, ssim_val


# ===========================================================================
# VERIFICATION PIPELINE
# ===========================================================================

def run_verification(solver, targets, methods=('gs', 'gd'), verbose=True):
    """
    Run both solvers on all targets, collect metrics.
    """
    results = []

    for name, P_target in targets.items():
        for method in methods:
            if verbose:
                print(f"\n{'='*55}")
                print(f"  Target: {name}  |  Method: {method.upper()}")
                print(f"{'='*55}")

            t0 = time.time()
            if method == 'gs':
                r = solver.solve_gerchberg_saxton(P_target, n_iter=300,
                                                  n_restarts=3,
                                                  verbose=verbose)
            else:
                r = solver.solve_gradient_descent(P_target, n_iter=500,
                                                  lr=0.05, reg_smooth=0.01,
                                                  verbose=verbose)
            elapsed = time.time() - t0

            metrics = compute_metrics(r['achieved'], P_target)

            if verbose:
                print(f"    NRMSE={metrics['nrmse']:.4f}  "
                      f"SSIM={metrics['ssim']:.4f}  "
                      f"eff={metrics['efficiency']:.4f}  "
                      f"time={elapsed:.1f}s")

            entry = {
                'target_name': name,
                'method':      method,
                'P_target':    P_target,
                'P_achieved':  r['achieved'],
                'phase_screen': r['phase_screen'],
                'metrics':     metrics,
                'convergence': r['convergence'],
                'time_s':      elapsed,
            }
            if 'phi_loops' in r:
                entry['phi_loops'] = r['phi_loops']
            results.append(entry)

    return results


# ===========================================================================
# PLOTTING
# ===========================================================================

def plot_inverse_holography(results, solver,
                            fname='results/inverse_holography.png'):
    print(f"\n  Plotting inverse holography -> {fname}")

    # Group by target
    target_names = []
    seen = set()
    for r in results:
        if r['target_name'] not in seen:
            target_names.append(r['target_name'])
            seen.add(r['target_name'])

    n_targets = len(target_names)
    fig = plt.figure(figsize=(28, 5 * n_targets + 6))
    gs = GridSpec(n_targets + 2, 5, figure=fig, hspace=0.5, wspace=0.35,
                  left=0.04, right=0.97, top=0.95, bottom=0.03)

    x_nm = np.linspace(-solver.L / 2, solver.L / 2, solver.N) * 1e9
    ext = [x_nm[0], x_nm[-1], x_nm[0], x_nm[-1]]

    # Banner
    ax = fig.add_subplot(gs[0, :])
    ax.axis('off')
    banner = (
        "INVERSE HOLOGRAPHY VIA SQUID ARRAY\n"
        + chr(9473) * 60 + "\n"
        f"He-4 at {solver.T_beam*1e3:.1f} mK  "
        f"lambda={solver.lam*1e9:.1f} nm  "
        f"SQUID: {solver.squid.N_loops}x{solver.squid.N_loops}  "
        f"prop={solver.prop_distance_lam} lambda = {solver.z*1e9:.1f} nm\n"
        f"Device: {_device}\n\n"
        "Columns: [Target] [GS result] [GD result] "
        "[GD phase screen] [Convergence]"
    )
    ax.text(0.5, 0.5, banner, transform=ax.transAxes, ha='center',
            va='center', fontsize=10, fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='#e3f2fd',
                      edgecolor='#1565c0', linewidth=2))

    for ti, tname in enumerate(target_names):
        row = ti + 1
        gs_r = [r for r in results
                if r['target_name'] == tname and r['method'] == 'gs']
        gd_r = [r for r in results
                if r['target_name'] == tname and r['method'] == 'gd']

        # Target
        ax = fig.add_subplot(gs[row, 0])
        P = gs_r[0]['P_target'] if gs_r else gd_r[0]['P_target']
        im = ax.imshow(P.T / (P.max() + 1e-30), extent=ext,
                       cmap='inferno', origin='lower')
        ax.set_title(f'Target: {tname}', fontsize=10, fontweight='bold')
        ax.set_xlabel('x (nm)', fontsize=8)
        ax.set_ylabel('y (nm)', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)

        # GS result
        if gs_r:
            ax = fig.add_subplot(gs[row, 1])
            I = gs_r[0]['P_achieved']
            im = ax.imshow(I.T / (I.max() + 1e-30), extent=ext,
                           cmap='inferno', origin='lower')
            m = gs_r[0]['metrics']
            ax.set_title(f"GS  SSIM={m['ssim']:.3f}\n"
                         f"eff={m['efficiency']:.3f}  "
                         f"t={gs_r[0]['time_s']:.1f}s",
                         fontsize=9, fontweight='bold')
            ax.set_xlabel('x (nm)', fontsize=8)
            ax.set_ylabel('y (nm)', fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.8)

        # GD result
        if gd_r:
            ax = fig.add_subplot(gs[row, 2])
            I = gd_r[0]['P_achieved']
            im = ax.imshow(I.T / (I.max() + 1e-30), extent=ext,
                           cmap='inferno', origin='lower')
            m = gd_r[0]['metrics']
            ax.set_title(f"GD  SSIM={m['ssim']:.3f}\n"
                         f"eff={m['efficiency']:.3f}  "
                         f"t={gd_r[0]['time_s']:.1f}s",
                         fontsize=9, fontweight='bold')
            ax.set_xlabel('x (nm)', fontsize=8)
            ax.set_ylabel('y (nm)', fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.8)

        # Phase screen (GD)
        if gd_r:
            ax = fig.add_subplot(gs[row, 3])
            ph = gd_r[0]['phase_screen']
            im = ax.imshow(ph.T, extent=ext, cmap='twilight_shifted',
                           origin='lower', vmin=-np.pi, vmax=np.pi)
            ax.set_title('GD phase screen', fontsize=9, fontweight='bold')
            ax.set_xlabel('x (nm)', fontsize=8)
            ax.set_ylabel('y (nm)', fontsize=8)
            plt.colorbar(im, ax=ax, shrink=0.8, label='rad')

        # Convergence
        ax = fig.add_subplot(gs[row, 4])
        if gs_r:
            ax.semilogy(gs_r[0]['convergence'], 'b-', lw=1.5,
                        label='GS (NRMSE)', alpha=0.8)
        if gd_r:
            ax.semilogy(gd_r[0]['convergence'], 'r-', lw=1.5,
                        label='GD (loss)', alpha=0.8)
        ax.set_xlabel('Iteration', fontsize=9)
        ax.set_ylabel('Error', fontsize=9)
        ax.set_title(f'{tname} convergence', fontsize=9, fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Summary bar chart
    ax = fig.add_subplot(gs[n_targets + 1, :3])
    gs_ssims = [r['metrics']['ssim'] for r in results if r['method'] == 'gs']
    gd_ssims = [r['metrics']['ssim'] for r in results if r['method'] == 'gd']
    x_pos = np.arange(n_targets)
    w = 0.35
    ax.bar(x_pos - w / 2, gs_ssims, width=w, color='#1565c0',
           edgecolor='k', alpha=0.8, label='GS')
    ax.bar(x_pos + w / 2, gd_ssims, width=w, color='#c62828',
           edgecolor='k', alpha=0.8, label='GD')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(target_names, fontsize=9)
    ax.set_ylabel('SSIM', fontsize=10)
    ax.set_title('SSIM Comparison: GS vs Gradient Descent',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 1.05)
    for i, (s1, s2) in enumerate(zip(gs_ssims, gd_ssims)):
        ax.text(i - w / 2, s1 + 0.02, f'{s1:.3f}', ha='center', fontsize=7)
        ax.text(i + w / 2, s2 + 0.02, f'{s2:.3f}', ha='center', fontsize=7)

    # Summary text
    ax = fig.add_subplot(gs[n_targets + 1, 3:])
    ax.axis('off')
    lines = ["SUMMARY\n" + chr(9472) * 40]
    for r in results:
        m = r['metrics']
        lines.append(f"{r['target_name']:<8} {r['method'].upper():<3}  "
                     f"SSIM={m['ssim']:.3f}  eff={m['efficiency']:.3f}  "
                     f"{r['time_s']:.1f}s")
    ax.text(0.03, 0.97, '\n'.join(lines), transform=ax.transAxes,
            fontsize=8, fontfamily='monospace', va='top',
            bbox=dict(boxstyle='round', facecolor='#e8f5e9',
                      edgecolor='#2e7d32'))

    plt.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close()


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("=" * 65)
    print("  INVERSE HOLOGRAPHY VIA SQUID ARRAY")
    print("=" * 65)

    squid = SQUIDArray(N_loops=32, N_grid=256, L_grid=400e-9)

    # -- Build mutual inductance matrix --
    print("\n" + "-" * 65)
    print("MUTUAL INDUCTANCE MODEL")
    print("-" * 65)
    M = squid.build_inductance_matrix(wire_width_frac=0.1)
    squid.info()

    solver = InverseHolographySolver(
        squid_array=squid, N=256, L=400e-9, T_beam=1e-3,
        prop_distance_lam=20.0,
    )

    print(f"  He-4 at {solver.T_beam*1e3:.1f} mK  "
          f"lambda={solver.lam*1e9:.1f} nm  v={solver.v:.3f} m/s")
    print(f"  Propagation: {solver.prop_distance_lam} lambda "
          f"= {solver.z*1e9:.1f} nm")
    print(f"  Device: {_device}")

    # -- Roundtrip validation --
    print("\n" + "-" * 65)
    print("ROUNDTRIP VALIDATION")
    print("-" * 65)
    validate_roundtrip(solver)

    # -- Target patterns (with corner rounding) --
    N, L = 256, 400e-9

    # Raw targets
    raw_targets = {
        'spot':   target_single_spot(N, L),
        'line':   target_line(N, L),
        'dots':   target_grid_of_dots(N, L, n_dots=3),
        'ring':   target_ring(N, L),
        'letter': target_letter(N, L, letter='H'),
    }

    # Band-limit targets to the SQUID array's spatial bandwidth, then
    # round corners and apply gentle blur.  The band-limiting is the
    # critical step: it removes spatial frequencies the phase screen
    # cannot modulate, so the optimizer doesn't waste energy chasing
    # unreachable content (which just bleeds into high-order diffraction).
    targets = {}
    for name, P in raw_targets.items():
        targets[name] = smooth_target(P, N_loops=squid.N_loops,
                                      corner_radius=0.03, sigma=2)

    print(f"\n  Target conditioning: band-limited to {squid.N_loops}-loop "
          f"Nyquist + corner rounding + sigma=2 blur")

    # -- Solve --
    print("\n" + "=" * 65)
    print("INVERSE HOLOGRAPHY -- PHASE RETRIEVAL")
    print("=" * 65)
    results = run_verification(solver, targets, methods=('gs', 'gd'))

    # -- Current mapping for GD solutions --
    print("\n" + "=" * 65)
    print("CURRENT REQUIREMENTS (flux -> drive currents via M^{-1})")
    print("=" * 65)
    print(f"  {'Target':<10} {'I_max (uA)':>12} {'I_rms (uA)':>12} "
          f"{'Roundtrip err':>14}")
    print("  " + "-" * 55)
    for r in results:
        if r['method'] != 'gd':
            continue
        if 'phi_loops' not in r:
            # GD results stored in the run_verification wrapper don't carry
            # phi_loops — recompute from phase screen by sampling at loop centers
            continue
        cm = squid.current_map_summary(r.get('phi_loops', r['phase_screen']))
        print(f"  {r['target_name']:<10} {cm['I_max_uA']:12.4f} "
              f"{cm['I_rms_uA']:12.4f} {cm['roundtrip_err']:14.2e}")

    # -- Plot --
    plot_inverse_holography(results, solver)

    # -- Summary table --
    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    print(f"  {'Target':<10} {'Method':<6} {'NRMSE':>8} {'SSIM':>8} "
          f"{'Contrast':>10} {'Efficiency':>12} {'Time(s)':>8}")
    print("  " + "-" * 65)
    for r in results:
        m = r['metrics']
        print(f"  {r['target_name']:<10} {r['method'].upper():<6} "
              f"{m['nrmse']:8.4f} {m['ssim']:8.4f} "
              f"{m['contrast']:10.4f} {m['efficiency']:12.4f} "
              f"{r['time_s']:8.1f}")

    print(f"\nOutput: results/inverse_holography.png")


if __name__ == "__main__":
    main()
