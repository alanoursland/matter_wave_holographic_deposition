# T41 - Multi-device field with correlated sidelobes and process variation

Implemented 2026-07-12.

## Purpose

T40 showed that the calibrated 3x3 holographic actuator can print one
functional T39 resistor coupon with contacts smaller than its demagnified
electrode pitch. T41 tests whether that result survives array replication,
overlapping sidelobes, and the qualified T39 process variation.

Study: `src/t41_multidevice_correlated_field.py`.

Artifacts:

- `results/t41_multidevice_correlated_field.json`;
- `results/t41_multidevice_correlated_field.npz`;
- `results/t41_multidevice_correlated_field.png`.

## Field construction

Four T39 resistor/monitor coupons are placed in a 2x2 array on an 800 nm
prepared SOI field. Each device receives the three calibrated 3x3 T40
expected-thickness kernels by stage translation, for twelve exposures total.
Deposited intensities add incoherently between shots. The complete thickness
field is assembled before connected-component and electrical extraction, so
neighboring sidelobes can fill contacts, create residue, or bridge devices.

The reference dose sweep uses 300 nm device pitch. After selecting dose, the
same random realizations are evaluated at 300, 275, and 250 nm pitch. The
220 nm resistor mesa makes 250 nm the tightest tested pitch with a remaining
physical trench between adjacent channel mesas.

## Correlated process model

T41 applies the joint T39 qualification corner:

```text
registration sigma                 2.0 nm
correlated edge-displacement sigma 2.5 nm
dose coefficient of variation     10%
roughness correlation length       5.0 nm
```

For registration, dose, and edge displacement, 70% of variance is common to
corresponding exposures across the array and 30% is device-local. This creates
batch-like failures without increasing the single-device marginal variance.

The isolated control uses each device's identical realized local morphology
but removes contributions from neighboring hologram kernels. Full-field minus
isolated behavior therefore measures sidelobe overlap without changing random
draws.

## Dose result

Forty full fields at each 300 nm-pitch setpoint give:

```text
dose   device yield   all-four yield   intra-device bridge   failure corr.
1.45      10.6%            0.0%               0.0%              0.149
1.75      75.0%           52.5%               0.0%              0.445
2.05      97.5%           95.0%               0.0%              0.457
2.35      96.9%           92.5%               5.0%              0.623
2.65      57.5%           37.5%              62.5%              0.528
```

Low dose produces correlated opens. High dose converts the common sidelobe
structure into within-device bridges. The optimum tested setpoint is 2.05.
Positive failure correlation means array yield does not generally equal the
fourth power of average device yield; failures tend to cluster by realization.

## Packing result

One hundred fifty identical random realizations at dose 2.05 give:

```text
pitch   neighbor dose in contacts   device yield   all-four yield   inter-device bridges
300 nm          0.000%                 99.67%          98.67%               0%
275 nm          0.108%                 99.83%          99.33%               0%
250 nm          0.637%                100.00%         100.00%               0%
```

The 250 nm result is 150 successes in 150 trials. Its exact 95% binomial
confidence interval is 97.57-100%, so the simulation does not establish a
literal perfect yield. Across all 600 device instances, the corresponding
lower confidence bound is 99.39%.

The isolated control at the same 250 nm realizations gives 99.67% device yield
and 98.67% all-four yield. The small neighbor contribution therefore removes
the two isolated-control array failures rather than causing a new failure.
No intra-device or inter-device bridge occurs at the selected pitch and dose.

Qualified electrical statistics are:

```text
median device current at 1 V              55.9 uA
95th-percentile finite contact resistance  9.33 kohm
mean neighbor dose fraction in contacts    0.637%
95th-percentile neighbor dose fraction     0.755%
```

## Residual risk

At the selected 250 nm condition, 21.4% of all conducting pixels lie outside
the nominal contact windows. They remain disconnected from contacts in every
qualified realization, so they do not fail the present DC gate. They are still
real deposited halo material and could obstruct adhesion, alter surface
leakage, or connect after a later layer.

The edge model is a first-order correlated contour displacement, and the
global field repeats the existing monitor-leakage model rather than solving a
new twelve-terminal surface PDE. The reported array yield is simulated, not
measured.

## Impact

T41 advances the result from one tolerant device to a compact four-device
field. At the tested T39 process limits, correlated errors do not destroy the
array, and modest sidelobe overlap is beneficial at 250 nm pitch. The useful
dose window remains bounded on both sides: opens below it and bridges above it.

The next layer-level question is the disconnected halo. A two-layer array
experiment should deposit a subsequent routed or dielectric pattern over this
same morphology and determine whether nominally harmless first-layer residue
creates interlayer shorts or adhesion failures.
