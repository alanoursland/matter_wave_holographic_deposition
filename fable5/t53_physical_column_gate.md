# T53 — explicit electrostatic column-family gate

**Verdict: candidate family falsified.**

T50–T52 retained an ideal 80× mapping and supplied the chromatic coefficient
as an input. T53 replaces both assumptions with an explicit five-plate
electrostatic geometry, KinoPulse Laplace solves, and direct charged-particle
ray transfer.

## Falsification contract

The candidate has a 30 keV entrance, one swept focus electrode, and three
plates held at 29.990 kV to establish a nominal 10 eV landing region. The
focus electrode is swept from -15 to +29 kV in 2 kV steps. For every voltage,
the object plane is allowed to move to the position that gives exactly
`|M| = 1/80`; this is more favorable than fixing the package length first.

A member passes only if the exact 80× image:

- lies beyond the final electrode (`z >= 0.55 mm`),
- has kinetic energy `10 ± 1 eV`, and
- is no more than 10 mm from the object.

## Result

No swept voltage passes. On the converged 65×65×169 field grid, 15 voltages
form an inverted 80× first-order image with a real upstream object plane, but
all such images lie only 3.46–57.53 µm after the first landing plate—inside
the landing stack. Their kinetic energies are 240.07–365.27 eV and their
object-to-image lengths are 22.60–58.71 mm.

The member closest in landing energy uses a +9 kV focus electrode:

- image plane: 49.11 µm,
- image energy: 240.07 eV,
- column length: 26.29 mm,
- magnification: exactly -1/80 by the favorable object-distance solve,
- post-final-electrode gate: fail,
- 10 eV gate: fail by 230.07 eV,
- 10 mm length gate: fail by 16.29 mm.

Its local finite-energy focus derivative corresponds to an effective
`C_c = 16.7 µm`, but this is only a diagnostic at the invalid 240 eV image.
It does **not** validate the 1 mm coefficient assumed for a 10 eV image in
T50–T52.

## Numerical checks

The vacuum Laplace solve uses unit coefficient because a constant permittivity
factor does not change the solution and scaling by epsilon-zero prevented the
iterative residual from resolving a 10 eV remainder after 30 kV deceleration.
The project wrapper now exposes KinoPulse's stagnation window. Both voltage
bases converge to relative residual `1e-8`.

- Maximum-principle overshoot is 0.000163 V on the fine grid, far below the
  1 eV landing tolerance.
- The 49- and 65-point transverse grids give 242.42 and 240.07 eV for the
  +9 kV diagnostic image; the 230 eV gate miss is insensitive to that change.
- Direct finite-amplitude ray probes from 1 to 100 nm change the inferred
  column length by 0.28% and the image energy by 0.031 eV.
- The reusable tracer enforces electrostatic energy conservation exactly and
  independently agrees with the on-axis variational transfer calculation.

## Interpretation

This result makes the current implementation less feasible: a single focus
plate followed by a same-potential deceleration stack cannot simultaneously
provide the assumed demagnification and soft landing. It does not rule out a
multistage column. The next credible design must separate high-energy
demagnification, controlled deceleration, and a low-energy final imaging lens,
then derive the complete transfer and chromatic coefficient from that geometry.

Artifacts:

- `results/t53_physical_column_gate.json`
- `results/t53_physical_column_gate.npz`
- `results/t53_physical_column_gate.png`
- `src/t53_physical_column_gate.py`
- `src/iqs/experiments/electrostatic_column.py`
