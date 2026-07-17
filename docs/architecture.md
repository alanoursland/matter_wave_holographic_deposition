# Code architecture

The active implementation is organized under `src/iqs`. Historical module
names remain as compatibility wrappers, but reusable code must not import
from them.

## Dependency direction

```text
constants / numerics
        |
sources / actuators / lattices / deposition
        |
holography
        |
pipelines
        |
experiments / plotting / scripts
```

The two integrated configurations are:

- `iqs.pipelines.patterned_substrate.PatternedSubstratePipeline`, the
  historical v9 patterned-substrate workflow.
- `iqs.pipelines.holographic_caging.HolographicCagingPipeline`, the
  historical v10 source-to-holography workflow.

The original class names remain aliases for compatibility.

## Public modules

- `iqs.sources`: source parameters and coherent matter-wave source models.
- `iqs.actuators`: ideal, electrostatic, FEM, and field-solver actuator models.
- `iqs.holography`: inverse solver, target generation, metrics, and validation.
- `iqs.lattices`: diamond-network models and density mapping.
- `iqs.deposition`: stochastic surface and electrical extraction models.
- `iqs.experiments`: reusable sweeps that return data.
- `iqs.plotting`: figure-producing helpers.

Top-level `coherent_matterwave_beam.py`, `inverse_holography.py`, `sim_v9.py`,
and `sim_v10.py` exist only so older reports, notebooks, and tests keep their
original imports. New studies should import from `iqs`.

## Direct execution

Repository-local demo entry points live in `scripts/`. They bootstrap the
`src` directory so an installed package is not required.
