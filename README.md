# Inertial Spin Model: Flocking with Positional Forces

A C++17 implementation of a two-dimensional, constant-speed Inertial Spin Model (ISM) flock, with stochastic dynamics, metric neighbor alignment, and position-dependent force fields. The repository provides a headless parameter-sweep executable and an SFML interactive visualizer.

## Features

- Constant-speed agents with inertia, friction, and correlated Gaussian noise.
- Velocity alignment with neighbors selected by metric distance.
- Uniform cell-list index for faster metric-radius neighbor searches.
- Confinement and other force-field implementations: generalized power-law, harmonic, quartic, central radial/swirl, and file-table forces.
- A constrained stochastic velocity-Verlet-style update.
- CSV metrics and periodic state snapshots for batch runs.
- Interactive sliders, FPS display, reset control, and a full-box view in the visualizer.

## Requirements

- CMake 3.20 or newer.
- A C++17 compiler.
- GSL development libraries.
- SFML graphics, window, and system components to build the visualizer.

GSL is required even for a headless build. SFML is optional when configuring with `ENABLE_VISUALIZATION=OFF`.

## Build

From the repository root:

```sh
cmake -S . -B build -DENABLE_VISUALIZATION=ON
cmake --build build --config Release
```

For a headless batch-only build:

```sh
cmake -S . -B build -DENABLE_VISUALIZATION=OFF
cmake --build build --config Release
```

The targets are `ism_visual` (when visualization is enabled) and `ism_batch`. Depending on the CMake generator, executables may be in `build/` or a configuration subdirectory such as `build/Release/`.

## Interactive visualizer

Run from the repository root, optionally passing a configuration file:

```sh
./build/ism_visual config.txt
```

On Windows, the executable is typically `.\build\ism_visual.exe config.txt` (or under `build\Release\` for a multi-configuration generator). If no file is supplied, the visualizer looks for `config.txt` in the current working directory.

### Visualizer unit convention

The visualizer uses normalized units with

\[
v_0=J=\eta=1,
\]

where \(v_0\) is the speed, \(J\) the alignment coupling, and \(\eta\) the friction. Its time and length units are \(t_0=\eta/J=1\) and \(L_0=v_0\eta/J=1\). Thus its `TimeStep`, `BoxSize`, and well parameters are dimensionless values in these units.

The controls represent:

- **Inertia \(\Omega\):** \(\Omega=J\chi/\eta^2\). The visualizer passes this value as the inertia; with the fixed scales above, its numerical value equals \(\chi\).
- **Temperature \(T/J\):** the dimensionless noise temperature.
- **Well radius \(r_0^*\):** \(Jr_0/(v_0\eta)\).
- **Well strength \(s^*\):** \(v_0s/\eta\) for the generalized well used by the visualizer.
- **Well exponent \(n\):** the power-law exponent. The potential scales as \(r^{2n}\): \(n=1\) is harmonic and \(n=2\) is quartic.

Friction, interaction strength, and maximum speed are fixed at 1 in the visualizer; `Friction`, `Interaction`, and `MaxSpeed` values in the config file do not change them. The visualizer accepts `VisualOmega` (or `Omega`) and `VisualTemp` for its dimensionless inertia and temperature. For compatibility, it also falls back to `VisualInertia`/`InertiaInit` and `TempInit`. The sample config sets `VisualOmega 1.0` separately so the batch sweep's `InertiaInit` does not unintentionally start the visualizer at an extremely small inertia.

The sidebar shows the fixed scales and normalized parameters. Drag the sliders to change parameters in real time; **Reset flock** respawns the agents while retaining the current slider values. Toggle **Radius: ON** and click a bird to show its light-grey unit interaction-radius circle; its current neighbors are outlined. The simulation advances 100 physics steps per rendered frame and the window is capped at 120 FPS. The square simulation box spans \([-B,B]\) in each direction, where \(B\) is `BoxSize`; the display fits the whole box while preserving its aspect ratio. In open-boundary mode, the view also expands as needed to keep particles that leave the box visible.

Set `PeriodicBoundary 1` in the config to use periodic boundaries in either executable; `0` or an omitted key selects open boundaries. In open mode, positions are not wrapped or reflected and neighbor distances are Euclidean. In periodic mode, positions are wrapped into \([-B,B)\) after every step and neighbor distances use minimum-image separation.

Set `RandomInitialConditions 1` to initialize positions uniformly by area inside a disk of radius `min(WellRadius, BoxSize)`, with independent, uniformly random headings at the model's fixed speed. Set it to `0` (or omit it) to retain the swirling initial ring. This option applies to both the visualizer and batch executable.

The visualizer uses a fixed metric interaction cutoff of 1.0. It does not provide a topological-neighbor mode. Existing position-dependent force fields are not periodic; periodic mode therefore rejects active positional fields. For the visualizer, set `WellStrength 0` to use periodic boundaries; the well controls are hidden in that mode.

## Batch parameter sweeps

The batch executable requires an output directory and a parameter file:

```sh
./build/ism_batch output config.txt
```

On Windows:

```powershell
.\build\ism_batch.exe output config.txt
```

The runner makes a grid sweep over `InertiaInit` to `InertiaFinal` and `TempInit` to `TempFinal`, using the corresponding `InertiaSteps` and `TempSteps`. Each pair produces a separate `run_NNN` subdirectory under the requested output directory. Inertias and temperatures are linearly spaced. The runner currently assigns seeds starting at 1337 and increments them for each run, so the config file's `Seed` value does not control the batch sweep.

Batch simulations use the parameter values in the config file as supplied; unlike the visualizer, they do not impose \(v_0=J=\eta=1\). Use a consistent unit system for `TimeStep`, `Inertia`, `Friction`, `Interaction`, `Temperature`, and force parameters.

### Configuration keys

The sample [`config.txt`](./config.txt) includes the main settings:

| Key | Meaning |
|---|---|
| `InertiaInit`, `InertiaFinal`, `InertiaSteps` | Inertia sweep range and number of values. |
| `TempInit`, `TempFinal`, `TempSteps` | Temperature sweep range and number of values. |
| `Seed` | Seed in the config; used by the visualizer, but overridden by the batch runner. |
| `NumBirds` | Number of agents. |
| `RandomInitialConditions` | `1` uses area-uniform positions inside the well and independent random headings; `0` or omitted uses the swirling initial ring. |
| `TotalSteps` | Integration steps per batch run. |
| `LogSteps` | Write one metrics row every this many steps. |
| `SaveSteps` | Write a full state snapshot every this many steps. |
| `BoxSize` | Half-width \(B\) of the square simulation box; its coordinate range is \([-B,B]\). |
| `PeriodicBoundary` | `1` enables position wrapping and periodic neighbor distances; `0` or omitted uses open boundaries. |
| `TimeStep` | Integration step in the selected unit system; in the visualizer this is dimensionless \(\Delta t^*=J\Delta t/\eta\). |
| `Inertia`, `Friction`, `Interaction`, `Temperature`, `MaxSpeed` | Batch-run dynamical parameters. The sweep runner supplies `Inertia` and `Temperature` from the sweep values. |
| `WellStrength`, `WellRadius`, `WellExponent` | Generalized power-law confinement parameters. A positive strength activates the field. |
| `HarmonicK`, `QuarticK`, `RadialStr`, `SwirlStr` | Optional force-field parameters read by the batch simulation. |

The batch entry point currently sets `QuarticK` to 0.5 for open runs (0 for periodic runs) and sets `HarmonicK`, `RadialStr`, and `SwirlStr` to zero before each sweep. The generalized well is controlled by `WellStrength`, `WellRadius`, and `WellExponent`. Because these positional fields are not periodic, periodic runs fail with a clear error if any positional field is active; set `WellStrength` to 0 for periodic runs.

Neighbor selection is metric: agents within distance 1.0 are included, excluding the agent itself. A uniform cell list limits distance checks to nearby cells and is rebuilt when forces are refreshed; periodic mode wraps the cell indices across the box edges. The cutoff is currently fixed in both the batch runner and visualizer; topological neighbor selection is not implemented. Open boundaries do not impose reflecting walls; use positional confinement if the flock should remain near the origin.

## Output files

Each batch run directory contains:

- `metrics.csv` — rows sampled at `LogSteps`, with total spin, mean spin, mean velocity components, tangential speed, mean radius, center of mass, and total energy.
- `env_step_N.csv` — model parameters, including the periodic-boundary flag, saved with each state snapshot.
- `birds_step_N.csv` — per-agent position, velocity, acceleration, constraint multipliers, force, and noise values.
- `rng_step_N.bin` — GSL random-number-generator state saved alongside each snapshot. A complete restart/load workflow is not currently implemented.

Snapshots are written at step 0 and then every `SaveSteps`. Metric rows include the simulation step index; convert steps to time using `step * TimeStep`.

## Implementation notes

The time-step update is a constrained stochastic velocity-Verlet-style scheme. It automatically subdivides requested time steps when the predicted deterministic velocity change or an eight-standard-deviation estimate of the velocity noise is large. This allows much smaller inertia values and stronger fields without clipping forces; exceptionally stiff, non-finite, or statistically extreme updates still report an error, and the visualizer pauses so parameters can be adjusted.