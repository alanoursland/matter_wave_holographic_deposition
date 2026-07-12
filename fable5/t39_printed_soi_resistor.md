# T39 - Complete printed SOI resistor coupon

Implemented 2026-07-12.

## Purpose

T33 demonstrated stochastic multi-material deposition, T34 established a
sourced p+ Si/Al-1.5%Si contact model, and T35-T38 established isolated SOI
geometry and surface-leakage bounds. T39 composes those pieces into the first
complete device-level coupon: a working two-terminal resistor and an adjacent
isolation monitor produced by sequential metal patterns on a prepared
substrate.

Study: `src/t39_printed_soi_resistor.py`.

Reusable circuit extractor: `src/iqs/deposition/electrical.py`.

Artifacts:

- `results/t39_printed_soi_resistor.json`;
- `results/t39_printed_soi_resistor.npz`;
- `results/t39_printed_soi_resistor.png`.

## Coupon

The prepared substrate contains:

- a 220 nm by 75 nm, 20 nm-radius p+ Si(100) SOI resistor mesa;
- a 75 nm square isolated monitor mesa;
- a 25 nm exposed-BOX trench between them;
- the T35 70 nm device-silicon layer and qualified trench surface.

The printer deposits three 55 nm Al-1.5%Si patterns in sequence: source,
drain, and monitor. The source and drain centers are 130 nm apart, leaving a
75 nm silicon channel between their nominal contact-window edges. Each pass
has independent registration, correlated edge roughness, and dose error. The
existing 11 nm FWHM projection response is retained.

The prepared SOI mesas, doping, etched contact windows, and post-metal anneal
are conventional process steps. T39 tests whether the holographic deposition
step can form the device terminals reliably; it does not claim to print the
SOI substrate itself.

## Electrical extraction

Realized metal above half the nominal 700 nm thickness defines each active
contact area. T34's sourced contact resistivity and p+ silicon resistivity give
the interface and spreading terms. The finite channel adds

```text
R_channel = rho_Si L / (W t) = 244 ohm.
```

KinoPulse separately solves sheet conduction between the complete resistor
mesa and monitor mesa. The solved geometry factor is `5.50246`; three trench
cuts agree to `1.57e-9` relative variation. At the `1e10 ohm/sq` surface
qualification floor, worst-case 1 V monitor leakage is `0.550 nA`, below the
`1 nA` criterion.

A realization passes only when all three contacts are below 10 kohm, total
source-to-drain resistance is below 25 kohm, monitor leakage is below 1 nA,
and no conducting metal component bridges nominal contacts.

## Nominal result

Thirty realizations at 1 nm registration sigma, zero imposed dose variation,
and zero imposed edge roughness give:

```text
functional yield                         100%
median contact resistance               8.61 kohm
95th-percentile contact resistance      8.86 kohm
median total device resistance         17.55 kohm
95th-percentile device resistance      17.92 kohm
median device current at 1 V            57.0 uA
5th-percentile device current           55.8 uA
worst-case monitor leakage              0.550 nA
median device/leakage current ratio     1.04e5
metal short rate                         0%
```

The two contacts consume about 98.6% of the applied voltage. The silicon
channel consumes only about 1.4%, so contact quality, not channel conduction,
sets device resistance in this geometry.

## Process window

The broad 30-replicate screen shows the failure boundary:

- 5 nm registration sigma reduces ideal-dose, smooth-edge yield to 50%;
- 10 nm registration sigma reduces it to 7%;
- 15% dose CV at 1 nm registration reduces smooth-edge yield to 90%;
- 5 nm edge roughness at ideal dose and 1 nm registration reduces yield to 73%.

A focused 100-replicate qualification sweep resolves the useful region. The
joint condition

```text
registration sigma <= 2 nm
correlated edge roughness sigma <= 2.5 nm
dose coefficient of variation <= 10%
```

produces 95% functional yield at the tested corner. The corresponding smooth
condition gives 98%. Increasing registration sigma to 3 nm reduces yield to
84-86% across the tested dose and roughness combinations.

Every condition passes the fixed isolation model and no metal-short population
controls the result. Failures are contact opens or contact resistance above
10 kohm. Registration is therefore the first process-control priority,
followed by dose uniformity and contact-edge roughness.

## Impact

This is the first integrated simulation in the project that begins with a
prepared semiconductor substrate, deposits separate terminal patterns, and
ends with a predicted terminal current, resistance, isolation current, and
statistical fabrication yield. Under the stated assumptions it is a working
device, not only a deposited image or disconnected contact array.

The claim remains conditional on nanoscale transfer of the micrometre-derived
specific contact resistivity and on the measured post-process surface sheet
resistance. The direct experimental counterpart is now well defined: print
this resistor/monitor coupon, measure source-drain resistance by Kelvin
methods, and measure monitor leakage at 1 V after the complete clean and
anneal sequence.
