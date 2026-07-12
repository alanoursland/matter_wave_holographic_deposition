# T35 - Etched p+ SOI / Al-1.5%Si isolated contact array

Implemented 2026-07-11.

## Purpose

T34 showed that good Al/Si contacts on a continuous p+ silicon wafer still
produce about 60 microamp of pair current at 1 V through the shared substrate.
T35 changes only the substrate architecture: the same contact geometry and
Al-1.5%Si interface are placed on p+ silicon mesas etched through an SOI
device layer to buried oxide.

Study: `src/t35_soi_isolated_contact_array.py`.

Reusable SOI extractor: `src/iqs/deposition/electrical.py`.

Artifacts:

- `results/t35_soi_contact_array.json`;
- `results/t35_soi_contact_array.npz`;
- `results/t35_soi_contact_array.png`.

## Selected substrate

The coupon uses:

- Si(100) SOI;
- 70 nm p+ device silicon;
- 145 nm buried SiO2;
- 75 nm square mesas at 100 nm pitch;
- complete mesa etch through the device layer to BOX;
- the T34 55 nm Al-1.5%Si contacts and post-metal anneal.

The 70 nm silicon / 145 nm BOX dimensions are based on an experimentally
studied partially depleted SOI stack (DOI
`10.1016/j.microrel.2010.01.015`). The device layer must additionally receive
the T34 boron concentration; that dopant profile is a fabrication requirement,
not a property of the cited commercial SOI wafer.

Thermal oxide bulk resistivity is set to the measured order of
`1e17 ohm cm` reported by Dumin, Fossum, and Todorov. Schwarzenbach et al.
(DOI `10.1109/S3S.2015.7333496`) report BOX leakage below `2 nA/cm2` at
`2 MV/cm`, breakdown fields above `10 MV/cm`, and defect density below
`0.2 cm^-2` for studied SOI material.

The Dumin resistivity was measured on much thinner thermal oxides and is used
only as an order-of-magnitude BOX bulk value. The independently measured SOI
leakage ceiling below provides the more direct cross-check.

## Leakage model

Two independent paths are placed in parallel.

The BOX path sends current vertically through each mesa footprint to an
ideally conducting handle wafer:

```text
R_BOX,pair = rho_BOX t_BOX / A_1 + rho_BOX t_BOX / A_2.
```

Omitting handle-wafer spreading resistance makes this an upper-current bound.

The exposed trench path uses a measured-or-qualified sheet resistance:

```text
R_surface,pair = R_sheet gap / facing_width.
```

No literature constant is assigned to the surface. Adsorbed water, residue,
sidewall damage, and vacuum history can move it by orders of magnitude. T35
therefore solves for the minimum sheet resistance required by the electrical
specification and then evaluates a separate qualification floor.

## Result

Thirty contact realizations retain the T34 projection and registration model.

```text
contact yield below 10 kOhm                    100%
median contact + spreading resistance          8.40 kOhm
worst contact + spreading resistance           8.62 kOhm

BOX-only worst pair leakage at 1 V              1.94e-23 A
required trench surface sheet resistance        3.16e9 ohm/sq
qualification floor                             1.00e10 ohm/sq
qualified worst pair leakage                    0.30 nA
isolation limit                                 1.00 nA
```

The published `2 nA/cm2` BOX ceiling, without extrapolating it to the much
lower operating field, scales to only `1.13e-19 A` over one 75 nm mesa. The
published `0.2 cm^-2` defect-density ceiling corresponds to `3.2e-10`
expected BOX defects in the complete 400 nm square simulation field.

The BOX is therefore not the practical leakage limit in an undamaged coupon.
The exposed trench surface controls the result. A measured sheet resistance
above `3.16e9 ohm/sq` is sufficient for the provisional 1 nA criterion; the
chosen `1e10 ohm/sq` qualification floor leaves a factor of about 3.3 margin.

## Impact

Unlike the bulk p+ substrate, this substrate supports nine independently
addressable contacts in the present DC model. The result is conditional, but
the condition is now experimentally sharp: fabricate a Kelvin/contact-chain
coupon plus adjacent isolated mesas, then measure contact resistance and
trench surface leakage after the exact clean, vacuum exposure, and anneal.

The remaining model risk is no longer nominal BOX conductivity. It is process
damage and contamination: incomplete mesa etch, conductive sidewall residue,
BOX pinholes, radiation damage, or field enhancement at etched corners. A
rounded-mesa geometry and a field-resolved corner study are the next useful
simulation refinement.

T36 performs that study. A 20 nm radius modestly reduces the resolved
mid-sidewall field, and all tested radii retain more than 10x sampled BOX
breakdown margin.
