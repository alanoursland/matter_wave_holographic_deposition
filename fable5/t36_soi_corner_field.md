# T36 - Field-resolved SOI corner and sidewall experiment

Implemented 2026-07-11.

## Purpose

T35 represented trench leakage with a scalar sheet resistance and did not
resolve electric-field concentration at etched mesa corners. T36 solves the
three-dimensional heterogeneous-dielectric field around two nearest-neighbor
SOI mesas and tests whether plan-view corner rounding materially changes the
sidewall or buried-oxide field.

Study: `src/t36_soi_corner_field.py`.

Artifacts:

- `results/t36_soi_corner_field.json`;
- `results/t36_soi_corner_field.npz`;
- `results/t36_soi_corner_field.png`.

## Geometry and solve

The T35 geometry is retained:

- 75 nm p+ Si mesas at 100 nm pitch;
- 25 nm nominal trench gap;
- 70 nm device layer;
- 145 nm BOX;
- neighboring mesas held at 0 V and 1 V;
- grounded conducting handle plane below the BOX;
- open side and top boundaries;
- relative permittivity 1.0 above the BOX and 3.9 inside it.

Plan-view corner radius is swept over `0`, `5`, `10`, and `20 nm`. KinoPulse
solves each geometry at `5 nm` and `2.5 nm` spacing with relative residual
below `1e-8`. The fine systems contain `121 x 89 x 121`, or approximately
1.30 million, grid unknowns.

The nominal parallel-edge field is

```text
E_nominal = 1 V / 25 nm = 40 MV/m.
```

Diagnostics use fixed physical regions on both grids:

- mid-sidewall: trench points from 10 to 60 nm above BOX;
- BOX interface: the first 5 nm of oxide beneath the trench;
- near-device: all free points around the mesas and top/bottom edges.

Both sampled maxima and 99th-percentile fields are reported. A maximum at an
ideal conductor edge is mesh dependent and is not treated as a converged
physical observable.

## Fine-grid result

```text
radius   sidewall p99 / E0   BOX p99 / E0   BOX sampled peak   breakdown margin
 0 nm          1.035             1.925          87.1 MV/m            11.48x
 5 nm          1.066             1.903          98.8 MV/m            10.12x
10 nm          1.025             1.882          96.1 MV/m            10.41x
20 nm          0.958             1.867          90.2 MV/m            11.08x
```

The breakdown comparison uses the greater-than-10-MV/cm BOX result from
Schwarzenbach et al., DOI `10.1109/S3S.2015.7333496`.

The 20 nm radius lowers mid-sidewall p99 field by 7.5% relative to the sharp
mesa and has the best sidewall convergence: its 5 nm to 2.5 nm p99 change is
2.6%. It is therefore the preferred tested geometry for surface control.

The BOX p99 appears to decrease by only 3.0% from 0 to 20 nm. That trend is
not resolved: the coarse-to-fine BOX p99 change remains 16% to 27%, larger
than the complete radius effect. The code records
`box_radius_effect_resolved = false` rather than presenting a false optimum.

All fine-grid BOX peaks retain at least 10.1x margin to the sourced breakdown
field. The non-monotonic sampled peak at 5 nm is a rasterized-edge result, not
evidence that small rounding physically makes breakdown worse.

## Leakage implication

T35's qualified surface result was `0.30 nA` at 1 V for a trench surface sheet
resistance of `1e10 ohm/sq`. Scaling that value by the fine-grid sidewall p99
enhancement gives approximately `0.31 nA` for the sharp geometry and
`0.29 nA` for 20 nm rounding. This is a field-weighted diagnostic, not a new
transport law; both remain below the provisional `1 nA` criterion.

The correct hardware qualification remains an electrical measurement at the
full 1 V bias on the etched geometry. A low-field sheet-resistance measurement
cannot be safely extrapolated through an unknown field-dependent surface
conduction mechanism.

## Impact

The T35 isolation result survives a three-dimensional field solve. Corner
rounding is useful but not transformative at this pitch: it reduces the
sidewall field modestly, while BOX field concentration remains comfortably
below breakdown for every tested radius.

The next simulation should not spend more global-grid computation chasing an
ideal edge singularity. It should use local mesh refinement or an imported FEM
mesh around the Si/BOX/air triple line, including a finite sidewall damage
layer. In parallel, the fabrication coupon should specify a 20 nm plan radius
and directly measure leakage after etch, clean, vacuum exposure, and anneal.
