# T40 - Holographic device-pattern resolution gate

Implemented 2026-07-12.

## Purpose

T39 demonstrated that ideal 55 nm Al-1.5%Si contact patterns form a working
SOI resistor. It still assumed an 11 nm projection response rather than
passing those contacts through the calibrated holographic printer. T40 closes
that gap and asks a stricter question:

```text
Can the holographic actuator deposit a functional T39 device when the desired
contacts are smaller than its control elements?
```

Study: `src/t40_holographic_device_resolution.py`.

Artifacts:

- `results/t40_holographic_device_resolution.json`;
- `results/t40_holographic_device_resolution.npz`;
- `results/t40_holographic_device_resolution.png`.

## Method

Source, drain, and monitor are optimized as three separate matter-wave
exposures on the T32 operating point:

- 14.4 nm He+ wavelength;
- 400 nm image field;
- 489 nm propagation distance;
- electrostatic phase response;
- 64 by 64 wave grid;
- 32 micrometre physical plate and nominal 80x demagnification.

The raw targets are the actual T39 55 nm contact windows, not targets softened
to the actuator bandwidth. Control grids of 3x3, 5x5, 7x7, and 9x9 are tested.
Only 3x3 is lowered through a field-resolved hardware model: its optimized
phase coordinates are converted to voltages with the measured T31 influence
matrix. Larger grids are information-capacity studies without corresponding
electrode FEM bases.

Each achieved intensity is upsampled to the T39 morphology grid and calibrated
so its target-window 95th percentile equals the nominal metal thickness at
unit dose. The three expected thickness fields are added. A 0.2-5.0 dose
sweep then applies the unchanged T39 gates:

- every contact below 10 kohm;
- total device resistance below 25 kohm;
- monitor leakage below 1 nA;
- no conducting bridge between contacts.

Image similarity alone cannot pass T40.

## Calibrated 3x3 result

The numerical 3x3 interpolation pitch is 133 nm. The actual T31 electrode
pitch is 8 micrometres, or 100 nm after nominal demagnification. A 55 nm
contact is therefore only 0.55 physical electrode pitches wide and 0.41
numerical control pitches wide.

Despite this, the three optimized exposures have a functional dose window:

```text
dose multiplier window                    1.45 to 2.75
minimum contact coverage at lower edge    84.2%
worst contact resistance                   9.80 kohm
total device resistance                   19.83 kohm
device current at 1 V                     50.4 uA
conducting area outside target windows     3.15%
metal bridge                              none
maximum calibrated electrode voltage       1.611 mV
available voltage limit                    2.000 mV
```

At higher dose, sidelobes eventually become conducting and bridge contacts.
At lower dose, the contact area is insufficient. The nonzero interval between
those failures is the relevant manufacturing margin.

Three independent optimizer initializations all recover the same 1.45-2.75
window, 84.2% minimum coverage, and 9.80 kohm worst contact. Their maximum
voltages span only 1.6097-1.6108 mV.

## Fidelity versus function

The 3x3 holograms are not faithful square printers:

```text
contact    target efficiency    raw SSIM    raw NRMSE
source          39.6%             0.545        0.615
drain           39.6%             0.545        0.615
monitor         38.5%             0.534        0.609
```

They pass because the device is tolerant to rounded, nonuniform contacts.
Enough metal lands inside each window to meet contact resistance, while the
remaining dose stays below the conducting threshold between windows. This is
a functional lithography result, not a claim of arbitrary 55 nm image
reproduction.

## Scaling sweep

Every tested ideal control count has a functional dose window:

```text
controls   image control pitch   contact/pitch   functional dose window
  3x3            133 nm              0.41             1.45-2.75
  5x5             80 nm              0.69             1.25-4.75
  7x7             57 nm              0.96             1.25-4.50
  9x9             44 nm              1.24             0.85->=5.00
```

The 9x9 upper limit is truncated by the tested dose range. The 7x7 and 9x9
optimizations each send about 2% of input power beyond the propagating or
accepted transverse band; the angular-spectrum solver removes it. Increasing
control count therefore does not remove the transport-band limit.

## Impact

T40 answers the recursive-printing question in a limited but concrete sense:
**the calibrated array can print a functional feature smaller than the array
elements used to print it.** Interference and thresholded material response,
not one-to-one geometric imaging, make that possible.

The result does not yet show that the printer can manufacture a smaller copy
of its complete electrode array. It demonstrates three isolated sub-element
contacts in one 400 nm field, one exposure at a time, using deterministic
expected thickness. The next defensible experiment is a multi-device field:
tile several T39 resistor coupons, include T39 registration/dose/roughness
statistics on the holographically delivered intensity, and measure array-level
all-device yield and correlated sidelobe bridges.
