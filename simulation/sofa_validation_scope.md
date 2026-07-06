# Scope: SOFA Validation of the Catheter Impedance MPC

Scoping note for validating the controller on a **SOFA** FEM continuum-catheter plant —
the highest-fidelity, most domain-credible target of the three engines considered
(RK4 → MuJoCo → SOFA). Companion to `mujoco_benchmark_scope.md`. Decision document.

## 1. Why SOFA (and why it beats MuJoCo here)

| | RK4 (paper) | MuJoCo (done) | **SOFA** |
|---|---|---|---|
| Catheter model | scalar double integrator | 8-link rigid PRB *approximation* | **FEM Cosserat beam (real continuum)** |
| Tendon actuation | lumped into `u` | hand-built moment arm | **Lagrangian cable constraint** |
| Tissue | analytic Kelvin–Voigt | analytic Kelvin–Voigt | FEM contact **or** analytic KV |
| Domain credibility | — | low (general robotics) | **high — the interventional/continuum standard** |

SOFA (Inria DEFROST) is *the* framework for soft/continuum medical-robot simulation. The
catheter is modeled as the actual deformable rod the paper's §II reduces *from* (it cites
Cosserat). For a catheter paper, a SOFA validation is what continuum-robotics reviewers
recognize — more than MuJoCo or Isaac.

**What it does NOT buy:** it still does not close the real gap (no **hardware**). SOFA
contact is FEM but still a model, not endocardium. Treat SOFA as the strongest *in-silico*
credibility step, not a hardware substitute.

## 2. Plugin stack

- **SOFA core** + **SofaPython3** (scene scripting + closed-loop control from Python).
- **`Cosserat`** (SofaDefrost) — FEM beam/rod = the catheter body. *Or* **`BeamAdapter`**,
  the classic SOFA catheter/guidewire plugin (interventional-radiology lineage).
- **`SoftRobots`** (SofaDefrost) — `CableConstraint`/`CableActuator` = the tendon, applied
  as a Lagrangian constraint (principled `J_κ`, not a hand-built moment arm).
- Contact: SOFA collision pipeline against a wall, **or** (recommended, for fidelity to the
  paper) apply the analytic Kelvin–Voigt force `F=k_t·pen+b_t·pen_rate` as an external
  force on the tip node — exactly as in the RK4 and MuJoCo benchmarks, giving exact
  `k_t=5000` and a like-for-like comparison.

## 3. Scene + control architecture

- **Catheter:** a Cosserat beam of length `L` discretized into `N` frames; bending stiffness
  set from the segment's flexural rigidity `EI` (calibrate to the RK4 `k_eff`).
- **Tendon:** one `CableConstraint` routed along the beam, offset `d_tend` from the
  neutral axis; a position/force actuator drives cable displacement/tension.
- **Tissue:** Kelvin–Voigt wall on a moving node (cardiac), analytic force on the tip.
- **Controller:** a `Sofa.Core.Controller` subclass with an `onAnimateBeginEvent` callback
  at the SOFA timestep; read tip-node position/velocity, run the **same** Impedance-MPC
  (impedance + offset-free + hard bound, ported from `catheter_mpc_mujoco.py`), map
  `F_mpc → cable command` via the measured `J_κ(κ)`, write the actuator input.
- Reuse wholesale from the MuJoCo port: the Λ-rescaled gains, the elasticity feedforward,
  the tendon-space hard bound `T ≤ (F_safe+F̂_fric)/J_κ(κ)`, the cardiac wall, the
  friction/hysteresis test.

## 4. What to measure (re-validate the MuJoCo results on FEM)

| Quantity | MuJoCo result to reproduce | SOFA expectation |
|---|---|---|
| `Λ(κ)`, `J_κ(κ)` | 3.5e-3, 0.087 | read from FEM mass/Jacobian; values differ, structure holds |
| Table II force safety | no-FC 0.60 N / FC 0.47 N | hard bound still the only safe controller |
| Table V cardiac | bound holds (0.47 N) | holds; FEM beam may add lateral modes |
| Offset-free free-space | 96 % (constant) / 42 % (hysteresis) | FEM cable friction → re-examine the hysteresis number |
| Constant-curvature kinematics | approx (over-curls) | FEM gives true non-constant curvature — a *better* test |

The highest-value SOFA-specific result: the **`J_κ` and offset-free-hysteresis findings on
a true FEM beam** — if they survive on real continuum mechanics, the paper-revision items
(force-constraint remap, hysteresis-aware observer) are much better supported.

## 5. Effort, risk, platform

- **Platform:** unlike Isaac, SOFA *can* run on macOS/Apple-Silicon, but install is the main
  risk. Options, fastest first: (a) **conda-forge** `sofa-*` packages if they cover
  `osx-arm64` + the Python ABI; (b) SOFA release binaries (Python-version-matched); (c)
  build from source (CMake + Qt + plugins) — multi-hour, error-prone on arm64. The
  SofaDefrost plugins (SoftRobots/Cosserat/BeamAdapter) are **not** core and may need a
  source build even if core SOFA is conda-installable.
- **Python ABI:** SofaPython3 is built against a specific Python (commonly 3.9–3.11); the
  system Python here is **3.13**, likely too new — a matched conda env is probably required.
- **Learning curve:** SOFA's scene-graph/component model and the `Controller` event loop are
  a steeper step than MuJoCo's flat API. Budget for it.
- **Estimate:** ~1 day if conda binaries cover everything; **multi-day** if any plugin needs
  a source build. A Linux box sidesteps most of the install risk.

## 6. Recommendation
**Not required for the current submission** — RK4 stands, MuJoCo cross-validation is already
optional. SOFA is the right investment **if** the goal becomes a high-credibility continuum
validation (journal/T-RO version, continuum-robotics venue, or a reviewer demanding FEM).
Then build it as its own milestone — ideally on Linux, or via a matched conda env here — and
port the controller wholesale from `catheter_mpc_mujoco.py`. The install-feasibility check on
this machine is recorded below.

## 7. Install-feasibility check on this machine (2026-06-13)

**Verdict: a quick spike to a runnable Python catheter scene is NOT feasible on this Mac
(Apple M1 Max, macOS arm64, system Python 3.13). A full source build is required — do it
on Linux or in a matched conda env, as its own milestone.**

Evidence (all tried here):
- **No pip wheel** — `pip install sofapython3` → *No matching distribution*.
- **conda-forge has core SOFA C++ only.** `conda create -c conda-forge sofa-devel sofa-gl
  python=3.11` succeeds and installs the SOFA C++ libraries (`libSofa.Component.*.dylib`,
  v25.12, osx-arm64), but provides **no `runSofa` binary and no Python binding** —
  `import Sofa` → *ModuleNotFoundError*. These are dev libs for *building plugins against
  SOFA*, not a runnable Python SOFA.
- **SofaPython3 binding: not packaged** on conda-forge (no `*sofa*python*` match) → source
  build, ABI-matched to the Python (3.13 here is likely too new; use 3.10/3.11).
- **Catheter plugins not packaged:** `SoftRobots`, `Cosserat`, `BeamAdapter` all return no
  conda-forge match → each needs a source build against the SOFA build. This is the wall.

**Path to actually do it (when ready):** a Linux box (sidesteps the arm64 build risk) or a
matched conda env on macOS, then build SOFA + SofaPython3 + Cosserat (or BeamAdapter) +
SoftRobots from source per the SofaDefrost build docs. Estimate **multi-day**, dominated by
the plugin builds, not the controller port (which transfers wholesale from
`catheter_mpc_mujoco.py`). Not warranted for the current submission (see §6).
