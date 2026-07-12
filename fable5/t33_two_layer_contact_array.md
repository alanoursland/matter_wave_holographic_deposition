# T33 - Two-species, two-layer contact-array coupon

Implemented 2026-07-11.

## Purpose

T20 modeled impact safety, neutralization, and thermal escape after landing,
but explicitly assumed near-unity sticking and had no evolving surface. T33
adds the missing 2.5D surface state and tests the first multi-species,
multi-layer fabrication coupon:

1. deposit nine `Si+` mesas;
2. deposit nine smaller `Al+` contact caps on the mesas.

This is a metal-on-semiconductor contact-array process coupon. It is not yet a
transistor, an electrically extracted diode model, or a complete phase plate.
It is the smallest stack that can fail by poor adherence, stochastic opens,
surface diffusion, registration error, cross-material selectivity, or shorts.

Study: `src/t33_two_layer_contact_array.py`.

Reusable surface model: `src/iqs/deposition/surface.py`.

Artifacts:

- `results/t33_two_layer_device.json`;
- `results/t33_two_layer_device.npz`;
- `results/t33_two_layer_device.png`.

## Surface model

Each material has a thickness field on an advancing surface. A layer update
includes:

- holographic projection blur;
- layer registration offset;
- Poisson-distributed atom count;
- material/interface-dependent sticking and thermal survival;
- conservative post-landing Gaussian diffusion;
- atomic-volume conversion from retained atoms to local thickness.

The model is explicitly not molecular dynamics. It does not predict crystal
structure, chemical reaction pathways, grain boundaries, conductivity, or
interface states.

## Coupon geometry

- field: `400 x 400 nm` on a `256 x 256` process grid;
- array pitch: `100 nm`;
- Si mesa width: `75 nm`;
- Al cap width: `55 nm`;
- layer thickness: `2 nm` each;
- projection blur: `11 nm` FWHM (`4.67 nm` sigma);
- substrate temperature: `300 K`;
- binding barrier: `1.2 eV`;
- hold time: `1000 s`;
- 30 stochastic replicates per process point.

The T20 worst-case `10 eV` landing gates pass for both species. At the assumed
barrier and temperature, post-landing thermal survival is `0.99993068`.

Al sticking on realized Si is swept over `0.5`, `0.7`, and `0.9`. Sticking on
exposed non-Si background is assumed `0.05`, representing selective contact
growth. Al diffusion and second-layer registration sigma are each swept over
`1`, `5`, and `10 nm`.

A replicate passes when all nine Si/Al contacts exceed `80%` layer coverage,
at least `95%` of Al mass lands on realized Si foundation, and no conducting
Al component bridges two contact sites.

## Fixed-dose process window

Commanded Al dose is initially calibrated for `90%` retention.

At `50%` actual Al sticking, every process point fails with nine open contacts.
This is an under-dose result, not an adhesion impossibility.

At `70%` sticking:

```text
diffusion 1 nm, registration 1 nm  -> 100% yield
diffusion 1 nm, registration 5 nm  -> 63% yield
diffusion 1 nm, registration 10 nm -> 13% yield
diffusion >= 5 nm                  -> 0% yield
```

At `90%` sticking:

```text
diffusion 1 nm, registration 1 nm  -> 100% yield
diffusion 1 nm, registration 5 nm  -> 80% yield
diffusion 1 nm, registration 10 nm -> 33% yield
diffusion 5 nm, registration 1 nm  -> 100% yield
diffusion 5 nm, registration 5 nm  -> 63% yield
diffusion 5 nm, registration 10 nm -> 40% yield
diffusion 10 nm                    -> 0% yield
```

No tested point produced a bridge short. The `45 nm` Al gap is generous enough
that failures are opens or insufficient Al-on-Si overlap before lateral
bridging becomes important.

## Measured-sticking dose calibration

At `1 nm` diffusion and `1 nm` registration, recalibrating commanded dose to
the actual sticking probability gives:

```text
Al sticking 0.5 -> 100% yield at 1.80x nominal dose
Al sticking 0.7 -> 100% yield at 1.29x nominal dose
Al sticking 0.9 -> 100% yield at 1.00x nominal dose
```

Thus sticking probability is a measurable throughput parameter as long as it
is nonzero and spatially selective. Surface diffusion and registration cannot
be repaired by scalar dose calibration because they move material to the wrong
place.

## Impact

The program now uses different deposit species and tracks an actual two-layer
surface. T20's landing result was necessary but insufficient: atoms can land
safely and remain thermally bound while the device still fails through
under-retention, diffusion, or misregistration.

The immediate process specifications suggested by this coupon are:

- characterize sticking and close dose around the measured value;
- hold effective Al diffusion below roughly `5 nm` for comfortable margin;
- target layer registration near `1 nm`; `5 nm` is marginal at nine-device
  all-pass yield;
- preserve interface selectivity so Al deposited outside Si does not remain.

These are model-derived requirements, not material measurements. The next
materials-facing work should replace assumed Si/Al sticking, diffusion, and
interface barriers with measured or literature-supported values for one
specific substrate preparation. Electrical extraction should then add contact
resistance, continuity, and leakage criteria.
