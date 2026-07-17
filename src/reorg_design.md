# Reorganization Design: Library-Oriented Integrated Quantum Substrate Codebase

> Implementation status (2026-07-16): the active source, holography, and v9/v10
> pipeline implementations now live under `iqs`; top-level historical modules
> are compatibility wrappers. Dedicated experiment, plotting, and demo modules
> expose the package-oriented API. See `docs/architecture.md` for the current
> dependency map. This document is retained as the migration rationale.

## Purpose

This document proposes a staged reorganization of the current simulation code into a reusable library with thin experiment scripts. The goal is not to rewrite the physics model first. The goal is to separate stable computational components from versioned demos so that future physics improvements can be made without repeatedly untangling pipeline code, plotting code, validation code, and experiment harnesses.

The current codebase has several strong reusable pieces:

- `coherent_matterwave_beam.py`: charged coherent matter-wave source model.
- `inverse_holography.py`: SQUID phase-screen model, inverse holography solvers, target generation, validation, plotting.
- `diamond_caging.py`: 2D diamond-network tight-binding model for flux-dependent Aharonov-Bohm caging.
- `sim_v9.py`: older integrated patterned-substrate pipeline with Floquet dressing, binding filter, caging, sweeps, plots, and demos.
- `sim_v10.py`: newer end-to-end source-to-holography-to-caging pipeline using CMWB and inverse holography.

The primary architectural problem is that `sim_v9.py` and `sim_v10.py` are doing double duty as both versioned experiments and importable libraries. The reorganization should make versioned scripts depend on the library, not the other way around.

## Target Package Layout

```text
iqs/
  __init__.py

  constants.py
  types.py
  config.py

  numerics/
    __init__.py
    device.py
    propagation.py
    metrics.py
    normalization.py

  sources/
    __init__.py
    species.py
    cavity.py
    coherent_matterwave.py

  devices/
    __init__.py
    squid_array.py
    inductance.py

  holography/
    __init__.py
    solver.py
    targets.py
    validation.py

  floquet/
    __init__.py
    spatial_dressing.py
    binding_filter.py
    validation.py

  lattices/
    __init__.py
    diamond.py
    lieb.py
    density_mapping.py

  pipelines/
    __init__.py
    holographic_caging.py
    patterned_substrate.py

  experiments/
    __init__.py
    pressure_sweep.py
    sideband_selectivity.py
    disorder_robustness.py
    decoherence_sweep.py
    ablation.py

  plotting/
    __init__.py
    pipeline_plots.py
    holography_plots.py
    caging_plots.py
    sweep_plots.py

scripts/
  run_v9_demo.py
  run_v10_demo.py
  run_inverse_holography_demo.py
  run_cmwb_demo.py
```

The dependency direction should be:

```text
constants / types / numerics
        ↓
sources      devices      lattices      floquet
        ↓        ↓            ↓            ↓
        holography            caging adapters
              ↓                 ↓
              pipelines / experiments
                       ↓
                    plotting / scripts
```

Library modules should not import from `scripts/`, versioned demos, or plotting modules. Plotting and scripts may import anything below them.

## Design Principles

1. Preserve behavior first. Initial milestones should move code with minimal logic changes.
2. Remove version files from the import graph. `sim_v9.py` and `sim_v10.py` should become scripts or compatibility wrappers, not dependencies.
3. Separate physics models from pipeline adapters. For example, a diamond lattice model should not be responsible for deciding how a 2D deposition image becomes an initial lattice state.
4. Separate computation from plotting. Plotting should consume result objects; simulation should not create figures.
5. Replace ad hoc dictionaries with typed result objects where practical.
6. Make configuration explicit and serializable.
7. Keep validation close to the component being validated, but keep demos out of library modules.
8. Keep historical models, such as the Lieb lattice, but mark them as experimental or archival rather than default.

## Milestone 0: Freeze Current Behavior and Add Golden Runs

### Goal

Before moving code, capture enough baseline behavior to detect accidental changes during refactoring.

### Work

Create a small set of reproducible smoke tests or golden-run scripts using fixed seeds and reduced sizes where needed.

Recommended baselines:

- CMWB synchronization result for one electron case and one He+ case.
- Diamond spectrum validation at `phi = pi` for `4x4`.
- Diamond caging run on a simple synthetic density.
- Inverse holography roundtrip validation.
- v10 dots target run with caging disabled, if runtime is acceptable.

### Acceptance Criteria

- Baseline scripts run without needing earlier sim versions.
- Key scalar outputs are recorded in a simple JSON or text artifact.
- Refactor milestones can compare against these outputs with tolerances.

### Notes

This milestone does not require a full test suite. It only needs enough checks to catch obvious breakage.

## Milestone 1: Extract Shared Constants, Metrics, and Propagation

### Goal

Create the lowest-level reusable utilities and remove duplicated numerical code.

### New Modules

```text
iqs/constants.py
iqs/numerics/device.py
iqs/numerics/metrics.py
iqs/numerics/propagation.py
iqs/numerics/normalization.py
```

### Move

From `sim_v9.py` and `inverse_holography.py`:

- `hbar`, `k_B`, `m_He`, and other shared constants.
- `michelson_contrast`.
- `ssim_score`.
- `min_feature_size`.
- Angular-spectrum propagation logic.
- Torch device selection.

### Proposed API

```python
from iqs.numerics.propagation import AngularSpectrumPropagator

prop = AngularSpectrumPropagator(N=N, L=L, k0=k0, z=z, device=device)
psi_out = prop.forward(psi_in)
psi_back = prop.backward(psi_out)
```

### Acceptance Criteria

- `sim_v9.py` and `inverse_holography.py` use the same propagation implementation.
- Existing propagation outputs match baseline within numerical tolerance.
- Metrics functions are imported from one place.

### Risks

The two existing propagators may differ slightly in assumptions or transfer-function construction. Preserve the current behavior by matching the implementation used in `inverse_holography.py` first, then adapt v9 to it carefully.

## Milestone 2: Extract the Diamond Network Lattice Module

### Goal

Make the diamond caging model a standalone library component.

### New Modules

```text
iqs/lattices/diamond.py
iqs/lattices/density_mapping.py
```

### Move

From `diamond_caging.py`:

- `DiamondNetworkSimulator` to `iqs/lattices/diamond.py`.
- Hamiltonian construction.
- Spectrum validation.
- Exact time evolution.

Then split density mapping into a separate helper:

- Downsample density.
- Smooth density.
- Find local peaks.
- Load peaks onto A sites.
- Build localization mask around peaks.

### Proposed API

```python
from iqs.lattices.diamond import DiamondNetwork
from iqs.lattices.density_mapping import DensityPeakMapper

lattice = DiamondNetwork(Lx=16, Ly=16)
mapper = DensityPeakMapper(N_grid=256, n_peaks=16)
state = mapper.map_to_a_sites(density_2d, lattice)
result = lattice.evolve(state, phi=np.pi, T=40.0, dt=0.05)
```

### Acceptance Criteria

- `DiamondNetwork` can be used without any deposition image.
- The old image-to-caging behavior is preserved through `DensityPeakMapper`.
- Returned results still contain `times`, `fidelity`, `spread`, `localization`, `snapshots`, `phi`, `n_peaks`, and `flat_ok` for compatibility.
- `sim_v9.py` and `sim_v10.py` no longer need to import diamond code through versioned files.

### Risks

The current caging result visualizes A-site snapshots only. Keep this behavior initially, but document it clearly. Later, add full five-sublattice snapshot options.

## Milestone 3: Extract Source Models

### Goal

Turn the coherent matter-wave beam code into a reusable source package.

### New Modules

```text
iqs/sources/species.py
iqs/sources/cavity.py
iqs/sources/coherent_matterwave.py
```

### Move

From `coherent_matterwave_beam.py`:

- `SPECIES` and `add_species` to `species.py`.
- `CavityGeometry` to `cavity.py`.
- `CoherentMatterwaveBeam` to `coherent_matterwave.py`.
- Plotting functions to plotting or scripts.
- Standalone demos to `scripts/run_cmwb_demo.py`.

### Proposed API

```python
from iqs.sources import CavityGeometry, CoherentMatterwaveBeam

cavity = CavityGeometry(pressure_Pa=100)
source = CoherentMatterwaveBeam(
    species="He+",
    E_kinetic_eV=k_B * T_beam / e_C,
    B_field=0.01,
    cavity=cavity,
    dE_frac=0.01,
    beam_current_A=1e-6,
)
sync = source.synchronize(seed=42)
beam = source.build_beam(sync, N_grid=256, L=400e-9)
```

### Acceptance Criteria

- Source models can be imported without matplotlib.
- Plotting and demo execution are not triggered on import.
- v10 imports source classes from `iqs.sources`, not from the old file.
- Synchronization outputs match baseline for fixed seeds and parameters.

### Risks

`build_beam` imports `gaussian_filter` locally. This is acceptable, but dependencies should be made explicit in package setup.

## Milestone 4: Extract SQUID Device and Holography Solver

### Goal

Split the inverse holography file into device model, solver, target generation, and validation.

### New Modules

```text
iqs/devices/squid_array.py
iqs/devices/inductance.py
iqs/holography/solver.py
iqs/holography/targets.py
iqs/holography/validation.py
```

### Move

From `inverse_holography.py`:

- `SQUIDArray` to `devices/squid_array.py`.
- Mutual inductance and current conversion either remain on `SQUIDArray` or move to `devices/inductance.py`.
- `InverseHolographySolver` to `holography/solver.py`.
- Target generators to `holography/targets.py`.
- `bandlimit_target` and `smooth_target` to `holography/targets.py`.
- `validate_roundtrip` to `holography/validation.py`.
- Plotting to `plotting/holography_plots.py`.
- Main demo to `scripts/run_inverse_holography_demo.py`.

### Proposed API

```python
from iqs.devices import SQUIDArray
from iqs.holography import InverseHolographySolver, smooth_target

squid = SQUIDArray(N_loops=32, N_grid=256, L_grid=400e-9)
solver = InverseHolographySolver(squid_array=squid, grid=grid, beam=beam_params)
target = smooth_target(raw_target, N_loops=squid.N_loops, corner_radius=0.03, sigma=2)
result = solver.solve_gradient_descent(target, n_iter=500, lr=0.05)
```

### Acceptance Criteria

- Holography solver no longer imports `sim_v9.py`.
- Solver uses shared constants, metrics, and propagator modules.
- Current mapping remains available as a post-processing utility.
- Roundtrip validation works from the new module.

### Risks

The solver currently assumes a clean Gaussian input beam internally, but v10 mutates `psi_in_t`, `A_in_t`, `psi_in_np`, and `A_in_np` after construction. Replace this with an explicit method:

```python
solver.set_input_beam(psi_profile)
```

This preserves behavior while removing private attribute mutation from pipeline code.

## Milestone 5: Extract Floquet and Binding Stages from v9

### Goal

Make the Floquet and binding stages reusable without requiring `IntegratedQuantumSubstrate`.

### New Modules

```text
iqs/floquet/spatial_dressing.py
iqs/floquet/binding_filter.py
iqs/floquet/validation.py
```

### Move

From `sim_v9.py`:

- `validate_floquet` to `floquet/validation.py`.
- `stage2_floquet_dress_spatial` to `floquet/spatial_dressing.py`.
- `stage3_binding_filter` to `floquet/binding_filter.py`.

### Proposed API

```python
from iqs.floquet import SpatialFloquetDresser, BindingFilter

dresser = SpatialFloquetDresser(n_sidebands=2, device=device)
psi_dressed, floquet_info = dresser.dress(psi, phase_map, V_max=0.9)
psi_ads, binding_info = BindingFilter(n_vals=dresser.n_vals).apply(
    psi_dressed,
    n_resonant=0,
)
```

### Acceptance Criteria

- Floquet dressing can be used by both v9-style and v10-style pipelines.
- The caveat remains documented: this is a local dressed-state decomposition unless the substrate is explicitly time-periodic.
- Existing v9 outputs remain approximately stable.

### Risks

The current implementation depends on `self.N`, `self.fl_dim`, `self.n_vals`, and grid shape from `IntegratedQuantumSubstrate`. These should become explicit constructor arguments or inferred from input arrays.

## Milestone 6: Introduce Typed State and Result Objects

### Goal

Reduce dictionary sprawl and make pipeline boundaries explicit.

### New Module

```text
iqs/types.py
```

### Proposed Dataclasses

```python
@dataclass
class Grid2D:
    N: int
    L: float
    dx: float
    x: np.ndarray
    y: np.ndarray

@dataclass
class BeamState:
    psi: np.ndarray
    wavelength: float
    k0: float
    coherence_r: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass
class HolographyResult:
    phase_screen: np.ndarray
    loop_phases: np.ndarray | None
    target: np.ndarray
    achieved: np.ndarray
    metrics: dict[str, float]
    convergence: list[float]
    method: str

@dataclass
class FloquetResult:
    density: np.ndarray
    avg_pops: np.ndarray | None
    entropy: float | None
    adsorption_fraction: float | None
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass
class CagingResult:
    times: np.ndarray
    fidelity: np.ndarray
    spread: np.ndarray
    localization: np.ndarray
    snapshots: list[tuple[float, np.ndarray]]
    phi: float
    n_peaks: int
    flat_ok: bool

@dataclass
class PipelineResult:
    source: BeamState
    holography: HolographyResult
    density_actual: np.ndarray
    density_filtered: np.ndarray
    floquet: FloquetResult | None
    cage_pi: CagingResult | None
    cage_0: CagingResult | None
    metadata: dict[str, Any] = field(default_factory=dict)
```

### Acceptance Criteria

- New pipeline code returns dataclasses.
- Compatibility wrappers can still expose dict-like results for old plotting code during migration.
- Plotting functions are updated to consume dataclasses.

### Risks

Converting everything at once will be noisy. Start with internal dataclasses in new pipeline modules, then add compatibility adapters.

## Milestone 7: Build New Pipeline Modules

### Goal

Replace versioned integrated classes with reusable pipeline orchestration.

### New Modules

```text
iqs/pipelines/holographic_caging.py
iqs/pipelines/patterned_substrate.py
```

### Holographic Caging Pipeline

This is the v10-style pipeline:

```text
CMWB source
  → clean beam profile for inverse design
  → inverse SQUID holography
  → noisy source through optimized screen
  → optional Floquet filter
  → diamond caging at Phi=pi and Phi=0
```

### Patterned Substrate Pipeline

This is the v9-style pipeline:

```text
Kuramoto/noisy beam
  → hand-designed AB phase pattern
  → propagation
  → spatial Floquet dressing
  → binding filter
  → diamond caging
```

### Acceptance Criteria

- `HolographicCagingPipeline` does not import `sim_v9.py` or `sim_v10.py`.
- `PatternedSubstratePipeline` does not import `sim_v9.py`.
- Both pipelines share source, propagation, Floquet, metrics, caging, and plotting infrastructure.
- v10 demo can be recreated through `scripts/run_v10_demo.py`.
- v9 demo can be recreated through `scripts/run_v9_demo.py`.

### Risks

The optional Floquet stage in v10 currently reconstructs phase from `phase_screen` and density rather than using the propagated complex wavefunction. Preserve that behavior initially, but label it as a modeling choice. Later, add an option to pass the actual propagated complex wavefunction.

## Milestone 8: Move Experiments Out of Version Files

### Goal

Turn sweeps and ablations into reusable experiment functions.

### New Modules

```text
iqs/experiments/pressure_sweep.py
iqs/experiments/sideband_selectivity.py
iqs/experiments/disorder_robustness.py
iqs/experiments/decoherence_sweep.py
iqs/experiments/ablation.py
```

### Move

From `sim_v9.py` and `sim_v10.py`:

- Pressure sweep.
- Sideband selectivity sweep.
- Multi-species deposition.
- Disorder robustness.
- Decoherence sweep.
- Ablation study.
- Multi-target holography benchmark.

### Acceptance Criteria

- Experiments accept pipeline/config objects rather than hardcoding all parameters internally.
- Experiments return structured results without plotting as a side effect.
- Plotting is a separate call.

### Example

```python
sweep = run_pressure_sweep(
    base_config=config,
    pressures=[1e-3, 1.0, 100, 1000, 10000, 101325],
    target=target_dots,
)
plot_pressure_sweep(sweep)
```

## Milestone 9: Move Plotting to Dedicated Modules

### Goal

Keep visualization useful while removing it from core computation.

### New Modules

```text
iqs/plotting/pipeline_plots.py
iqs/plotting/holography_plots.py
iqs/plotting/caging_plots.py
iqs/plotting/sweep_plots.py
```

### Move

- `plot_v10_pipeline`.
- `plot_multi_target`.
- `plot_pressure_sweep`.
- `plot_pipeline`.
- `plot_selectivity`.
- `plot_multi_species`.
- `plot_disorder`.
- `plot_decoherence`.
- `plot_inverse_holography`.

### Acceptance Criteria

- Core modules can be imported without importing matplotlib.
- Plotting functions accept result objects or experiment result tables.
- Scripts call plotting explicitly.

### Risks

Some current plots rely on dictionary keys from old result structures. Add adapters or update plots after dataclass migration.

## Milestone 10: Convert Versioned Files into Scripts or Compatibility Wrappers

### Goal

End the use of versioned simulation files as library modules.

### Work

Convert:

- `sim_v9.py` → `scripts/run_v9_demo.py` plus optional `legacy/sim_v9.py` wrapper.
- `sim_v10.py` → `scripts/run_v10_demo.py` plus optional `legacy/sim_v10.py` wrapper.
- `inverse_holography.py` main block → `scripts/run_inverse_holography_demo.py`.
- `coherent_matterwave_beam.py` main block → `scripts/run_cmwb_demo.py`.

### Compatibility Strategy

For a transition period, old files may remain as wrappers:

```python
# legacy/sim_v10.py
from iqs.pipelines.holographic_caging import HolographicCagingPipeline
from iqs.holography.targets import target_grid_of_dots

# provide aliases or thin wrappers for old names if needed
IntegratedPipelineV10 = HolographicCagingPipeline
```

### Acceptance Criteria

- No library module imports `sim_v9.py` or `sim_v10.py`.
- Old demos still run through scripts.
- Public examples use the new package API.

## Milestone 11: Improve Documentation and Examples

### Goal

Make the project understandable as a library and as a research scaffold.

### New Documents

```text
docs/architecture.md
docs/physics_model_boundaries.md
docs/pipeline_v9_patterned_substrate.md
docs/pipeline_v10_holographic_caging.md
docs/caging_models.md
docs/holography_model.md
docs/source_model.md
```

### Important Documentation Boundaries

Clearly distinguish:

- Simulated tight-binding AB caging versus physical stabilization of deposited matter.
- Local spatial Floquet decomposition versus rigorous time-periodic Floquet theory.
- CMWB/Kuramoto source coherence model versus full beamline simulation.
- Phase-only holography feasibility versus physical SQUID fabrication feasibility.
- Designed clean-source performance versus noisy-source actual performance.

### Acceptance Criteria

- A new reader can identify which modules are core library, experiments, plotting, and scripts.
- Each major model states its assumptions and outputs.
- Known modeling jumps are explicit rather than buried in code comments.

## Milestone 12: Physics-Facing Follow-Up Refactors

### Goal

After the code is library-like, improve the scientific modeling without destabilizing the architecture.

### Candidate Improvements

1. Beam state model
   - Include angular spread, transverse coherence length, velocity spread, and emittance-like parameters.
   - Let `r` affect more than random phase noise.

2. SQUID/device constraints
   - Add current limits, phase quantization, noise, mutual inductance conditioning, and fabrication-scale constraints.
   - Penalize unreachable current maps during optimization.

3. Holography solver
   - Support constrained optimization over loop phases and currents.
   - Add multi-objective losses: target fidelity, efficiency, smoothness, current feasibility, robustness to source noise.

4. Density-to-lattice mapping
   - Support A-site loading, bridge-site loading, compact localized state loading, and direct lattice-state initialization.
   - Compare how caging depends on initialization.

5. Floquet/binding stage
   - Separate static spatial filter mode from explicit time-periodic Floquet mode.
   - Add a mode flag so plots and summaries do not overstate rigor.

6. Caging observables
   - Add full sublattice snapshots.
   - Add participation ratio, inverse participation ratio, leakage outside cages, and disorder ensemble statistics.

7. Validation
   - Make spectrum checks size-aware and flux-aware.
   - Add analytical checks for small lattices.

## Proposed Migration Order Summary

```text
0. Freeze behavior with golden runs.
1. Extract constants, metrics, device, propagation.
2. Extract diamond lattice and density mapping.
3. Extract source models.
4. Extract SQUID device and holography solver.
5. Extract Floquet and binding stages.
6. Introduce typed state/result objects.
7. Build new v9/v10-style pipeline modules.
8. Move sweeps and ablations into experiments.
9. Move plotting into plotting modules.
10. Convert versioned files into scripts/wrappers.
11. Add architecture and model-boundary docs.
12. Begin physics-facing improvements.
```

## Immediate First Pull Request

The first PR should be intentionally boring:

1. Create package skeleton.
2. Move constants and metrics.
3. Move device selection helper.
4. Move angular-spectrum propagator.
5. Update `inverse_holography.py` and `sim_v9.py` to import those shared pieces.
6. Run baseline smoke checks.

Do not change the physics in the first PR. The goal is to prove the code can be moved safely.

## Immediate Second Pull Request

The second PR should extract caging:

1. Move `DiamondNetworkSimulator` into `iqs/lattices/diamond.py`.
2. Add `DensityPeakMapper` in `iqs/lattices/density_mapping.py`.
3. Keep a compatibility wrapper in `diamond_caging.py` if desired.
4. Update v9 and v10 to import from the new location.
5. Verify that caging outputs still provide the same keys.

## Immediate Third Pull Request

The third PR should extract source and holography:

1. Move CMWB source components into `iqs/sources`.
2. Move SQUID array into `iqs/devices`.
3. Move inverse solver and targets into `iqs/holography`.
4. Add `solver.set_input_beam(...)` and update v10 to use it instead of mutating solver internals.
5. Keep old top-level files as import-compatible wrappers until scripts are ready.

## Final Desired State

After the reorganization, the main v10-style workflow should look like this:

```python
from iqs.config import HolographicCagingConfig
from iqs.pipelines import HolographicCagingPipeline
from iqs.holography.targets import target_grid_of_dots
from iqs.plotting.pipeline_plots import plot_holographic_caging

config = HolographicCagingConfig.default_v10()
pipeline = HolographicCagingPipeline.from_config(config)

target = target_grid_of_dots(config.grid.N, config.grid.L, n_dots=3)
result = pipeline.run(target, seed=42)

plot_holographic_caging(result, "results/v10_pipeline_dots.png")
```

The version number should live in configs and scripts, not in the core architecture.
