# T37 - Geometry-resolved SOI surface transport

Implemented 2026-07-12.

## Purpose

T35 estimated trench leakage with `R = Rsheet gap / width`. T36 resolved the
electrostatic field but could only scale that rectangular leakage estimate by
a local field statistic. T37 solves the actual DC surface-current problem and
integrates current across the trench:

```text
div(Gsheet grad(V)) = 0.
```

Study: `src/t37_soi_damage_transport.py`.

Artifacts:

- `results/t37_soi_damage_transport.json`;
- `results/t37_soi_damage_transport.npz`;
- `results/t37_soi_damage_transport.png`.

## Model

The BOX surface is represented as a uniform resistive sheet. The two p+ mesa
footprints are embedded Dirichlet regions at 0 V and 1 V. KinoPulse solves the
two-dimensional sheet potential with open outer boundaries; current is
integrated across three trench cross-sections.

The sheet solve is normalized to unit conductance, producing a dimensionless
geometry factor

```text
K = I Rsheet / V.
```

For any measured post-process sheet resistance,

```text
I = V K / Rsheet.
```

The T35 rectangular approximation gives `K = width / gap = 3.0`.

Conductive sidewall residue is modeled as an equipotential lateral extension
of each mesa by `0`, `2.5`, `5`, or `10 nm`. This is conservative: it assumes
zero vertical resistance between p+ silicon and residue at the BOX surface.
Real finite vertical resistance can only reduce pair leakage.

Plan-view mesa radii `0` and `20 nm` are solved at `5` and `2.5 nm` grid
spacing. Every accepted solve converges below `1e-9` relative residual.

## Results

Fine-grid extraction at 1 V:

```text
radius  residue  remaining gap  geometry K  required Rsheet  I at 1e10 Ohm/sq
 0 nm    0 nm       25 nm          4.158       4.16e9             0.416 nA
20 nm    0 nm       25 nm          3.756       3.76e9             0.376 nA

 0 nm    5 nm       15 nm          6.041       6.04e9             0.604 nA
20 nm    5 nm       15 nm          5.265       5.27e9             0.527 nA

 0 nm   10 nm        5 nm         19.159       1.92e10             1.92 nA
20 nm   10 nm        5 nm         12.982       1.30e10             1.30 nA
```

The clean geometry is already 25% to 39% more conductive than T35's
rectangular approximation because current spreads around the complete mesa
perimeter. The 20 nm radius reduces extracted leakage by about 10% for a clean
surface and 13% at 5 nm residue.

The clean and 5 nm residue cases have coarse-to-fine geometry-factor changes
below 3%. Both pass the `1 nA` limit at the `1e10 ohm/sq` qualification floor.

The 10 nm residue cases fail on the fine grid. Their remaining physical gap is
only 5 nm and coarse-to-fine changes are 33% to 43%, so their exact current is
not grid resolved. The failure conclusion is still conservative: sharp mesas
fail on both grids, and both fine-grid geometries exceed the limit.

The `2.5 nm` residue case is also flagged as grid-sensitive because it is one
fine cell and half a coarse cell. It is not used to set the process limit.

## Solver finding

An attempted 3D conformal-shell solve used a high-conductivity damage layer in
an almost insulating volume. The current KinoPulse Jacobi-preconditioned
elliptic path stagnated at the required conductivity contrast on the fine
grid. No nonconverged result was accepted.

The reduced sheet equation is the physically appropriate model for measured
surface sheet resistance and converges cleanly. The failed 3D attempt is still
useful library feedback: high-contrast coefficients need a stronger
preconditioner, domain elimination, or a mesh that contains only the conducting
manifold.

KinoPulse `0.1.0.dev2026071200` repaired false stagnation during improving PCG
plateaus. T38 revisits the conformal shell with that release and obtains
converged coarse and fine solutions. A 50-iteration plateau window is still
needed by the 1.3-million-point fine case.

## Impact

T35's `3.16e9 ohm/sq` requirement was optimistic. Geometry-resolved transport
raises the clean requirement to as much as `4.16e9 ohm/sq` and the qualified
5 nm residue requirement to `6.04e9 ohm/sq`.

The resulting fabrication specification is:

- retain the 20 nm plan-view corner radius;
- require post-process trench sheet resistance at least `1e10 ohm/sq` at 1 V;
- limit electrically conductive lateral sidewall residue to `5 nm` maximum;
- reject or rework coupons showing a conductive extension near `10 nm`.

The next useful device experiment is now the two-terminal isolated silicon
channel proposed after T35: two Al-1.5%Si contacts on one SOI channel, with an
adjacent isolated monitor mesa. That will separate contact resistance, channel
resistance, and neighbor leakage in one printable coupon.
