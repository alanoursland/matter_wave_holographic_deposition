"""
T21 — Closed-Loop Stochastic Printing (gen-2 architecture study)
=================================================================

Holographic deposition is inherently stochastic: each ion lands at a
position sampled from the (ensemble-averaged) arrival intensity, so an
open-loop exposure carries Poisson statistics — per-site dose error
~1/√N and feature-edge roughness ~√N.  But unlike photons, every ion
announces its arrival (secondary-electron detection, as practiced in
deterministic single-ion implantation).  This study simulates the
closed loop proposed in notes/design/gen2_architecture.md §6:

    expose a sub-dose batch → count accumulated arrivals per site →
    compute the remaining deficit → re-solve the hologram to aim at
    the deficit (warm-started) → expose the next batch → …

The hologram is the only actuator, so correction is possible only
within the deliverable band (T14): sites are defined at that scale.
Feedback therefore fixes *dose allocation* errors between resolvable
sites; sub-site shot noise remains — exactly the division of labor the
gen-2 note assigns to beam (where) vs surface chemistry (precision).

Metrics
-------
- RMS relative per-site dose error vs cumulative dose, open loop vs
  closed loop vs the 1/√N Poisson floor of a perfect open-loop source.
- Site defect rate (|dose − target| > tol) vs cumulative dose — the
  defect-rate-vs-throughput curve.

Acceptance (gen-2 note): feedback beats the √N scaling; a
defect-rate-vs-throughput curve exists.

Run:  python t21_closed_loop_printing.py
"""

import os
import time
import numpy as np
import torch
import matplotlib.pyplot as plt

from iqs.holography import (
    SQUIDArray, InverseHolographySolver,
    target_grid_of_dots, smooth_target,
)
from iqs.sources import sample_phase_noise
from iqs.numerics.device import get_device

_device = get_device()
os.makedirs('results', exist_ok=True)

# Fixed stage (v10 geometry, transport-aware conditioning)
N_GRID  = 256
L_AP    = 400e-9
N_LOOPS = 32
SIGMA_THETA = 0.079     # Q-limited coherence ceiling (phase2)
XI_PERP     = 50e-9
M_ENSEMBLE  = 50        # noise realizations for the arrival law


# ---------------------------------------------------------------------------
# Forward pieces
# ---------------------------------------------------------------------------

def make_solver():
    squid = SQUIDArray(N_loops=N_LOOPS, N_grid=N_GRID, L_grid=L_AP)
    solver = InverseHolographySolver(
        squid_array=squid, N=N_GRID, L=L_AP, T_beam=1e-3,
        prop_distance_lam=20.0)
    return solver


def arrival_law(solver, phase_screen, rng):
    """Ensemble-averaged arrival density for a hologram.

    Each ion sees one noise realization; the marginal arrival
    distribution over many ions is the intensity average (T6/T11).
    Returns a normalized (N, N) probability array.
    """
    screen_t = torch.tensor(phase_screen, dtype=torch.float64,
                            device=_device)
    T = SQUIDArray.phase_screen_to_transmission(screen_t)
    prof = solver.psi_in_np
    dx = solver.dx
    dens = np.zeros((N_GRID, N_GRID))
    for _ in range(M_ENSEMBLE):
        noise = sample_phase_noise(N_GRID, SIGMA_THETA,
                                   XI_PERP / (L_AP / N_GRID), rng)
        psi = prof * np.exp(1j * noise)
        psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx**2)
        psi_t = torch.tensor(psi, dtype=torch.complex128, device=_device)
        out = solver._propagate_torch(psi_t * T)
        dens += torch.abs(out).cpu().numpy()**2
    dens /= M_ENSEMBLE
    return dens / dens.sum()


def sample_arrivals(p, n_ions, rng):
    """Poisson-sample n_ions (expected) arrivals from probability map p."""
    return rng.poisson(n_ions * p)


# ---------------------------------------------------------------------------
# Sites: the deliverable-band scale at which the hologram can steer dose
# ---------------------------------------------------------------------------

def coarse_grain(A, site_px):
    n = A.shape[0] // site_px
    return A[:n * site_px, :n * site_px].reshape(
        n, site_px, n, site_px).sum(axis=(1, 3))


def site_metrics(counts, target_sites, total_dose, tol=0.2):
    """Per-site relative dose error over target sites + defect rate.

    A target site is in spec when its accumulated dose is within
    ±tol of its target share of the total dose so far.

    Also returns the *scale-calibrated* RMS: the residual after the
    best global dose rescale (least-squares α).  Real lithography
    calibrates exposure ("dose-to-size"), so a common-mode dose offset
    — e.g. the fixed fraction lost to band-limited haze between
    features — is not a pattern error; only the site-to-site shape
    residual is.
    """
    expected = target_sites * total_dose
    mask = target_sites > 0.5 * target_sites.max()   # the pattern sites
    c, e = counts[mask], expected[mask]
    rel_err = (c - e) / (e + 1e-30)
    rms = float(np.sqrt(np.mean(rel_err**2)))
    defect = float(np.mean(np.abs(rel_err) > tol))
    alpha = float(np.sum(c * e) / (np.sum(e * e) + 1e-30))
    rel_cal = (c - alpha * e) / (alpha * e + 1e-30)
    rms_cal = float(np.sqrt(np.mean(rel_cal**2)))
    return rms, defect, rms_cal


# ---------------------------------------------------------------------------
# Open- and closed-loop runs
# ---------------------------------------------------------------------------

def run_open_loop(solver, target, doses, seed, site_px):
    """One hologram solved once; accumulate dose; measure at checkpoints."""
    rng = np.random.default_rng(seed)
    np.random.seed(seed); torch.manual_seed(seed)
    sol = solver.solve_gradient_descent(target, n_iter=500, lr=0.05,
                                        verbose=False)
    p = arrival_law(solver, sol['phase_screen'], rng)
    target_sites = coarse_grain(target / target.sum(), site_px)

    counts = np.zeros_like(coarse_grain(p, site_px))
    rows = []
    prev = 0.0
    for D in doses:
        batch = sample_arrivals(p, D - prev, rng)
        counts += coarse_grain(batch.astype(float), site_px)
        prev = D
        rms, defect, rms_cal = site_metrics(counts, target_sites, D)
        rows.append((D, rms, defect, rms_cal))
    return rows, sol


def _predicted_rms(counts, p_sites, batch, target_sites, D_after):
    """Expected site RMS after adding `batch` ions from allocation
    p_sites (shot noise ignored — it's what we're trying to average
    away; the prediction uses only the deterministic allocation)."""
    expected_counts = counts + batch * p_sites
    exp_target = target_sites * D_after
    mask = target_sites > 0.5 * target_sites.max()
    rel = (expected_counts[mask] - exp_target[mask]) / (exp_target[mask]
                                                        + 1e-30)
    return float(np.sqrt(np.mean(rel**2)))


def run_closed_loop(solver, target, doses, seed, site_px,
                    n_iter_resolve=400, verbose=False):
    """Model-predictive closed loop over a growing hologram library.

    Each round: (1) count accumulated arrivals per site; (2) solve a
    correction hologram aimed at the remaining deficit (warm-started,
    same conditioning as the base target) and add it to the library;
    (3) expose the next batch from whichever library hologram minimizes
    the *predicted* post-batch site error.

    Two properties matter.  The base hologram never leaves the library,
    so the closed loop can only match or beat open loop in expectation
    (open loop is always an available policy).  And alternating between
    library holograms realizes a *convex combination* of their dose
    allocations over time — the mechanism by which feedback beats any
    single hologram's systematic allocation error, not just shot noise.
    The gate uses only each hologram's site allocation, which a real
    machine can calibrate (witness exposures / accumulated statistics).
    """
    rng = np.random.default_rng(seed)
    np.random.seed(seed); torch.manual_seed(seed)
    sol = solver.solve_gradient_descent(target, n_iter=500, lr=0.05,
                                        verbose=False)
    phi = sol['phi_loops']
    p = arrival_law(solver, sol['phase_screen'], rng)

    tnorm = target / target.sum()
    target_sites = coarse_grain(tnorm, site_px)
    D_final = doses[-1]

    # hologram library: (full-grid law, site allocation)
    library = [(p, coarse_grain(p, site_px))]

    counts = np.zeros_like(target_sites)
    rows = []
    prev = 0.0
    for i, D in enumerate(doses):
        batch = sample_arrivals(p, D - prev, rng)
        counts += coarse_grain(batch.astype(float), site_px)
        prev = D
        rms, defect, rms_cal = site_metrics(counts, target_sites, D)
        rows.append((D, rms, defect, rms_cal))
        if verbose:
            print(f"    dose {D:9.3g}: site RMS = {rms:.4f}, "
                  f"defects = {defect:.3f}")

        if D >= D_final:
            break
        # --- solve a correction hologram for the remaining deficit ---
        D_next = doses[i + 1]
        deficit_sites = np.maximum(target_sites * D_next - counts, 0.0)
        if deficit_sites.sum() > 0:
            deficit = np.kron(deficit_sites, np.ones((site_px, site_px)))
            full = np.zeros_like(target)
            full[:deficit.shape[0], :deficit.shape[1]] = deficit
            # same conditioning as the base target (the earlier heavy
            # site-scale blur wrecked the re-solves)
            full = smooth_target(full, N_loops=N_LOOPS, sigma=2,
                                 k0=solver.k0, L=L_AP, z=solver.z)
            # Fresh random init (not warm start): warm-started corrections
            # inherit the base solution's allocation bias, and a library
            # of same-biased holograms cannot cancel it by mixing —
            # actuator diversity is what lets the hull cover the target.
            res = solver.solve_gradient_descent(
                full, n_iter=n_iter_resolve, lr=0.05, verbose=False)
            phi = res['phi_loops']
            p_corr = arrival_law(solver, res['phase_screen'], rng)
            library.append((p_corr, coarse_grain(p_corr, site_px)))

        # --- predictive gate over the whole library, sub-batched ---
        # Splitting the next exposure into sub-batches lets successive
        # gate decisions realize finer convex combinations of the
        # library allocations (no extra hologram solves needed).
        n_sub = 4
        B_sub = (D_next - D) / n_sub
        D_run = D
        sim_counts = counts.copy()
        plan = []
        for _ in range(n_sub):
            D_run += B_sub
            scores = [_predicted_rms(sim_counts, ps, B_sub, target_sites,
                                     D_run)
                      for _, ps in library]
            best = int(np.argmin(scores))
            plan.append(best)
            sim_counts = sim_counts + B_sub * library[best][1]
        # expose the planned sub-batches (arrivals sampled per sub-batch)
        for b in plan[:-1]:
            sub = sample_arrivals(library[b][0], B_sub, rng)
            counts += coarse_grain(sub.astype(float), site_px)
        # the final sub-batch is exposed at the top of the next loop
        # iteration from p, so hand it the last plan entry
        p = library[plan[-1]][0]
        prev = D_next - B_sub
        if verbose:
            print(f"      controller: library of {len(library)}, "
                  f"plan {plan}")
    return rows


# ---------------------------------------------------------------------------
# Study
# ---------------------------------------------------------------------------

def run_study(seeds=(0, 1, 2), verbose=True):
    solver = make_solver()

    raw = target_grid_of_dots(N_GRID, L_AP, n_dots=3)
    target = smooth_target(raw, N_loops=N_LOOPS, corner_radius=0.03,
                           sigma=2, k0=solver.k0, L=L_AP, z=solver.z)

    # Site scale = deliverable band: f_c ≈ k0·sin(atan(L/z))·L/2π ≈ 3.1
    # cycles/aperture → half-period ≈ N/(2 f_c) ≈ 41 px; use 32 px so the
    # 3×3 dots land in distinct sites.
    site_px = 32
    # dose checkpoints double as closed-loop batch boundaries
    doses = np.unique(np.round(np.logspace(2.5, 5.5, 13))).astype(float)

    open_runs, closed_runs = [], []
    for s in seeds:
        if verbose:
            print(f"\n  seed {s}: open loop")
        o, _ = run_open_loop(solver, target, doses, seed=s,
                             site_px=site_px)
        open_runs.append(o)
        if verbose:
            print(f"  seed {s}: closed loop")
        c = run_closed_loop(solver, target, doses, seed=s,
                            site_px=site_px, verbose=verbose)
        closed_runs.append(c)

    def collect(runs, col):
        return np.array([[r[col] for r in run] for run in runs])

    out = {
        'doses':       doses,
        'open_rms':    collect(open_runs, 1),
        'open_def':    collect(open_runs, 2),
        'open_cal':    collect(open_runs, 3),
        'closed_rms':  collect(closed_runs, 1),
        'closed_def':  collect(closed_runs, 2),
        'closed_cal':  collect(closed_runs, 3),
        'site_px':     site_px,
    }
    return out, target


def plot_study(out, fname='results/t21_closed_loop.png'):
    print(f"\n  Plotting → {fname}")
    d = out['doses']
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    for key, color, label in [('open_cal', 'r', 'open loop'),
                              ('closed_cal', 'b', 'closed loop')]:
        m = out[key].mean(axis=0)
        s = out[key].std(axis=0)
        ax.loglog(d, m, f'{color}o-', lw=2, label=label)
        ax.fill_between(d, m - s, m + s, color=color, alpha=0.15)
    # Poisson floor for a perfectly aimed source: rel err = 1/√(D·t_s)
    # with t_s the mean target-site share — use the open-loop first point
    # to anchor the slope.
    anchor = out['open_cal'].mean(axis=0)[0]
    ax.loglog(d, anchor * np.sqrt(d[0] / d), 'k--', alpha=0.6,
              label=r'$1/\sqrt{N}$ scaling')
    ax.set_xlabel('cumulative dose (ions)')
    ax.set_ylabel('RMS site-dose error (dose-calibrated)')
    ax.set_title('Feedback vs shot noise', fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend()

    ax = axes[1]
    for key, color, label in [('open_def', 'r', 'open loop'),
                              ('closed_def', 'b', 'closed loop')]:
        m = out[key].mean(axis=0)
        ax.semilogx(d, m, f'{color}s-', lw=2, label=label)
    ax.set_xlabel('cumulative dose (ions) — throughput')
    ax.set_ylabel('site defect rate (|err| > 20%)')
    ax.set_title('Defect rate vs throughput', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.suptitle('T21 — closed-loop stochastic printing '
                 '(dots target, v10 stage, transport-aware conditioning)',
                 fontweight='bold')
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {fname}")
    plt.close(fig)


def main():
    print("=" * 65)
    print("T21: CLOSED-LOOP STOCHASTIC PRINTING")
    print("=" * 65)
    t0 = time.time()
    out, _ = run_study()
    plot_study(out)

    d = out['doses']
    orms = out['open_rms'].mean(axis=0)
    crms = out['closed_rms'].mean(axis=0)
    ocal = out['open_cal'].mean(axis=0)
    ccal = out['closed_cal'].mean(axis=0)
    odef = out['open_def'].mean(axis=0)
    cdef = out['closed_def'].mean(axis=0)
    print(f"\n  {'dose':>10} {'RMS open':>10} {'RMS closed':>11} "
          f"{'cal open':>9} {'cal closed':>11} "
          f"{'def open':>9} {'def closed':>11}")
    print("  " + "-" * 78)
    for i in range(len(d)):
        print(f"  {d[i]:10.3g} {orms[i]:10.4f} {crms[i]:11.4f} "
              f"{ocal[i]:9.4f} {ccal[i]:11.4f} "
              f"{odef[i]:9.3f} {cdef[i]:11.3f}")

    # headline numbers
    gain = orms[-1] / (crms[-1] + 1e-30)
    gain_cal = ocal[-1] / (ccal[-1] + 1e-30)
    print(f"\n  Final-dose RMS improvement (closed/open): {gain:.2f}× "
          f"absolute, {gain_cal:.2f}× dose-calibrated")
    for rate in (0.1, 0.01):
        do = d[np.argmax(odef <= rate)] if np.any(odef <= rate) else None
        dc = d[np.argmax(cdef <= rate)] if np.any(cdef <= rate) else None
        print(f"  Dose to defect rate ≤ {rate:>4}: "
              f"open = {do if do else 'not reached'}, "
              f"closed = {dc if dc else 'not reached'}")
    print(f"\nDone in {time.time()-t0:.0f}s. "
          f"Output: results/t21_closed_loop.png")


if __name__ == '__main__':
    main()
