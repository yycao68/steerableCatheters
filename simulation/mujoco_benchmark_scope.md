# Scope: MuJoCo Catheter Benchmark

Scoping note for replacing / augmenting the RK4 benchmark (`catheter_benchmark.py`)
with a MuJoCo physics benchmark. Decision document, not an implementation.

## 1. What this buys — and what it does not

**Buys:**
- A contact solver and a *distributed-compliance* plant instead of a lumped scalar
  double integrator — so the single-DOF projection, the `Λ(κ)` claim, and the
  "disturbance state absorbs the residual" claim are tested against a plant that
  genuinely violates the modelling assumptions, not one built to satisfy them.
- The phenomena the paper currently **defers**: tendon–sheath friction (real
  hysteresis, not a constant `F_fric=0.45 N` bias), tendon backlash/slack, non-constant
  curvature under load, and off-tip contact.
- A directly *measurable* operational-space inertia `Λ(κ)` (via the MuJoCo mass matrix
  + tip Jacobian), turning the §VI-F(1) pole-invariance claim from an analytic check
  into a plant-grounded one.

**Does NOT buy:**
- It does **not** close the real gap, which is **no hardware**. MuJoCo contact is still
  a tuned spring-damper, not endocardium. A reviewer who wants wet-lab/bench evidence
  is not answered by a sim-engine swap.
- It does not make the *scalar reduced-order* claim more credible on its own — RK4 on
  the post-feedforward double integrator is already the honest fidelity for that claim.
  MuJoCo's value is in stressing the **uncancelled** terms, not the cancelled ones.

> Read this as: MuJoCo is worth doing to *harden the robustness story and pre-empt the
> "RK4 is toy physics" reviewer*, not as a substitute for hardware.

## 2. Model architecture (MJCF)

### 2a. Catheter body — recommended: `cable` elasticity plugin
MuJoCo ≥3.x ships `plugin="mujoco.elasticity.cable"` (libelasticity): an inextensible
elastic rod discretized as a chain of bodies with bending/twisting stiffness — a
Cosserat-like model, which is exactly the `\cite{rucker2011}` class the paper reduces
*from*. This is more faithful and less hand-tuning than pseudo-rigid-body (PRB).

- Discretize the single segment (length `L`) into ~20 elements.
- Set bending stiffness so aggregate flexural rigidity `EI` matches the nominal
  `k_eff = 50` restoring behaviour of the scalar model (calibrate by commanding a known
  tip load and matching tip deflection).
- **Fallback (PRB):** serial chain of N rigid links + torsional joint springs,
  `k_joint = N·EI/L`, joint damping for `b_eff`. Use only if the cable plugin's tendon
  coupling proves awkward.

### 2b. Tendon actuation
- One `<tendon><spatial>` routed through `<site>`s along the rod, anchored at the tip
  (positive pull → curvature `κ>0`, matching the paper's sign convention).
- A `<motor tendon=...>` actuator supplies the tendon force `u` (the controller output).
- **Tendon `frictionloss`** → models sheath friction → *real* hysteresis. **Tendon
  slack / `springlength`** → backlash. These replace the constant `F_fric` bias and are
  the headline new disturbance the Kalman state must absorb.

### 2c. Tissue wall + cardiac motion
- Wall = a `box`/`plane` geom on a `slide` joint (contact-normal axis).
- **Contact stiffness/damping map directly:** `solref="-5000 -40"` gives `k_t=5000 N/m`,
  `b_t=40 N·s/m` (negative solref = direct stiffness/damping), matching the paper's
  Kelvin–Voigt wall exactly. `solimp` tuned for minimal penetration.
- Cardiac motion: drive the wall slide joint with `y_wall(t)=A·sin(2πft)` (mocap or a
  position actuator), reusing the `0.3 mm/1 Hz` and `0.5 mm/1.2 Hz` cases of Table V.

### 2d. Sensing & timestep
- Tip position from a tip `<site>`; **no force sensor** (consistent with the paper). A
  `touch`/force sensor may be added for *ground-truth validation only*, not fed to the
  controller.
- Timestep: contact stiffness 5 kN/m needs `dt ≲ 1e-4 s` (or `implicitfast`/`implicit`
  integrator) for stable contact. Control loop stays at 500 Hz (sub-sample the physics).
- **Known gotcha (from project memory):** add `<compiler angle="radian"/>` or joint
  ranges come out ~53× too small.

## 3. Controller integration

Reuse the existing `ImpedanceMPC`/Kalman/QP code almost verbatim — only the plant I/O
changes:
- Extract scalar tip-normal `y, ẏ` from MuJoCo state (project tip site motion onto the
  contact normal).
- Get `J_κ` and `Λ(κ)` from `mj_jac` (tip Jacobian) + `mj_fullM` (mass matrix) each step
  — now *measured*, not analytic.
- Partial-physics feedforward now cancels the **known cable/joint stiffness+damping**;
  the larger residual (distributed dynamics, friction, contact) flows into `d̂`. This is
  the intended stress test of the "absorb the residual" design rule.
- Map `F_mpc → tendon force` via `J_κ`, apply through the tendon actuator.

## 4. Validation matrix — what gets re-checked and expected outcome

| Claim / table | Expected under MuJoCo | Risk |
|---|---|---|
| Table II hard-constraint story (only +FC ≤0.5 N; no-FC ~7 N) | Qualitatively **survives** — QP still caps controller-induced force | low |
| Table II peak/hold *numbers* | **Shift** ± with distributed compliance & contact tuning; hold error likely larger | expected, must re-report |
| Table V cardiac safety-mode | Bound holds for moderate wall; breach threshold **moves** with real contact damping | low-med |
| §VI-F(1) `Λ(κ)` invariance | **Stronger** — now measured from mass matrix; 1.4× variation re-derived | could change number |
| §VI-F(4) 7 mm/s velocity limit | **Shifts** with actual `b_t` contact transient | expected |
| §VI-F(3) 27 µs timing | **Unchanged** (controller-side, engine-independent) | none |
| Offset-free 98% free-space cut | **HIGHEST-VALUE NEW TEST** — does `d̂` absorb tendon hysteresis/backlash? | **high** — may degrade |

The bottom row is the point of the whole exercise: it's the strongest new evidence *and*
the biggest risk to an existing claim. If the Kalman state cleanly absorbs real
sheath-friction hysteresis, that's a much stronger paper. If it doesn't, you learn that
before a reviewer (or hardware) does.

## 5. Where claims could break (be honest up front)
- **Contact stiffness ceiling:** 5 kN/m is stiff for MuJoCo at large `dt`; if it needs an
  implicit integrator or tiny `dt`, the "stiff-contact tension" result is still valid but
  costs sim time.
- **Single-DOF projection leakage:** MuJoCo's multi-link state has lateral/coupling
  modes the 1-DOF controller can't see. Bounded → a *good* robustness result; if large →
  exposes the over-reduction the paper flags as future work (multi-segment).
- **Feedforward mismatch:** larger uncancelled residual could push `d̂` past the Kalman
  bandwidth → the §V-C "design rule" gets a real numeric test.

## 6. Effort & sequencing (~4–6 focused days)
1. **Model build + calibrate** (cable rod + tendon + wall, match `EI`/contact to scalar
   params): ~1.5–2 d. *Highest-risk step — contact + tendon tuning.*
2. **Controller integration** (state I/O, `mj_jac`/`mj_fullM`, tendon mapping): ~1 d.
3. **Re-run + re-validate Tables II, V, §VI-F; regenerate figures:** ~1–1.5 d.
4. **Reconcile changed numbers in `body.tex`, add a "MuJoCo verification" subsection or
   replace the RK4 results:** ~0.5–1 d.

Reuse: borrow MJCF/loop patterns from the verified `knee_rehab` (myoLeg) and
`dexterous_hand` (LEAP) sims rather than starting cold.

## 6b. Spike results (run 2026-06-13 — `catheter_mujoco.py`)

Step-1 spike built: an 8-link pseudo-rigid-body (PRB) tendon-driven catheter + a
Kelvin–Voigt wall on a z-slide (cardiac drive), MuJoCo 3.8.1, no cable plugin needed.
Generated model saved to `catheter_mujoco_model.xml`. Verdict: **GO with caveats.**

**Works (de-risked):**
- Tendon-driven bending, tip-normal `y/ẏ` readout, `Λ(κ)` from the mass matrix, and
  cardiac wall drive (cmd 0.30 mm → meas 0.29 mm) all run cleanly.
- `Λ(κ)` varies ~2.1× over a wide curvature sweep (paper's 1.4× is over a narrower
  `κ∈[2,25]` range — plausibly consistent).

**Contact stiffness — CALIBRATED to 5 kN/m (resolved 2026-06-13):**
- MuJoCo's regularized soft contact (`solref`/`solimp`) **cannot** faithfully represent a
  stiff 5 kN/m wall — direct-mode `solref="-5000 -40"` measured only ~100 N/m force/pen
  even with impedance→1, because the contact is a regularized constraint, not a literal
  spring.
- **Faithful fix (now in `catheter_mujoco.py`):** apply the tissue reaction
  **analytically** as `F = k_t·pen + b_t·pen_rate` on the tip (the paper's Kelvin–Voigt
  model, identical to the RK4 benchmark). MuJoCo carries the catheter; the wall geom is
  visual-only. Verified `F/pen = 4998 N/m` (target 5000) at every press point, exact by
  construction. The tip rests at µm-scale penetration (2–36 µm) exerting 0.01–0.18 N.
- **Numerical cost:** the analytic force is explicit on a near-massless tip, so stiff
  k_t=5000 needs `dt=2e-5` **plus joint `armature=2e-5`** (reflected/structural inertia)
  for stability; without these the tendon-driven tip joint NaNs. Worth flagging as a
  real cost of the MuJoCo route.

**Force-constraint J_κ remap — DEMONSTRATED (resolved 2026-06-13):**
- The controller commands tendon tension `T`; the tip force is `F_tip = J_κ(κ)·T` with
  `J_κ ≈ d_tend/L`. The scalar bound `u ≤ F_safe` implicitly assumes `J_κ = 1`.
- On a higher-transmission variant (arm 2.5 mm, `J_κ ≈ 0.10`) where the bound binds
  (`F_tip` reaches 0.67 N at the 8 N tendon limit), the three strategies give:
  | force constraint | tip force | verdict |
  |---|---|---|
  | (a) none | **0.67 N** | UNSAFE (33 % over) |
  | (b) scalar `T ≤ F_safe` | **0.012 N** | 2 % of budget — unusably over-conservative |
  | (c) `J_κ` remap `T ≤ (F_safe−b)/J_κ` | **0.49 N** | ~F_safe ✓ |
  The scalar bound is wrong in **both** directions; only the remap lands at `F_safe`. The
  remap inverts the local *affine* tendon→tip map (the map is mildly superlinear and not
  through the origin); the small residual vs 0.5 N is `J_κ`'s curvature dependence, which
  the controller removes by evaluating `J_κ(κ)` online.
- **Paper implication:** the force-safety constraint must be written
  `T ≤ (F_safe + F_fric)/J_κ(κ)`, configuration-dependent — a genuine revision to the
  force-constraint section (currently `u_ff + J_κᵀF_mpc ≤ F_safe + F̂_fric`, which is
  dimensionally right *if* `J_κ` is the true transmission, but the scalar benchmark runs
  with `J_κ = 1`). The qualitative "only the hard constraint is safe" story survives; the
  numbers and the constraint expression do not port directly.

**MPC ported with Λ-rescaled gains — DONE (`catheter_mpc_mujoco.py`, 2026-06-13):**
Closed-loop Impedance MPC running on the MuJoCo PRB plant at 500 Hz. Measured
`Λ=3.5e-3, J_κ=0.087, k_eff=8.4 N/m`. Three changes the port required:
- **Λ-rescaled gains:** `Kd→Kd·Λ=3.15`, `Dd→Dd·Λ=0.21` (keeping Kd=900 would put poles at
  ω≈500 rad/s, unstable at 500 Hz). MPC works in **tip-force units** so the bound
  `|F_mpc|≤F_safe` is dimensionally correct; the J_κ remap lives in the actuator map
  `T = −(ff+F_mpc)/J_κ`.
- **Elasticity feedforward is essential (new finding):** the Λ-rescaled gains alone are
  *softer* than the catheter's own `k_eff=8.4 N/m`, so without `ff=k_eff·y_ref` the tip
  reaches only ~27 % of the reference (approach RMS 1.77 mm). Adding the Layer-1 FF
  restores the double-integrator residual → **approach RMS 0.018 mm**. The FF matters
  *more* than the scalar model implied.
- **Result reproduces Table II on the physics plant:**
  | controller | approach RMS | peak contact F | hold err |
  |---|---|---|---|
  | MPC + force constraint | 0.018 mm | **0.467 N** (safe) | 1.41 mm |
  | MPC, no constraint | 0.018 mm | **0.595 N** (VIOLATES) | 1.41 mm |
  The offset-free integral winds up against the 1.5 mm-penetrating target and the
  unconstrained drive breaches 0.5 N; only the J_κ-remapped hard bound keeps it safe.
  Hold error 1.41 mm is contact-limited (the target needs ~7.5 N tip force, the catheter
  is tendon-limited to ~0.7 N) — matches the scalar paper's ~1.3–1.5 mm.
- **Caveat:** with gentle integral gain the constraint is *non-binding* (force stays well
  under 0.5 N) — the danger appears only under aggressive offset-free windup. So the
  "only the hard constraint is safe" claim is now **conditional on the offset-free
  aggressiveness and J_κ geometry**, a more nuanced statement than the scalar Table II.
- **Cardiac re-validation (Table V) on the physics plant — DONE, with an important
  caveat found in deep-check:** force-constrained MPC with a beating wall (tissue damping
  uses the *relative* tip–wall velocity). Bound **holds at all three conditions** — static
  0.467 N, 0.3 mm/1 Hz 0.468 N, 0.5 mm/1.2 Hz+noise 0.466 N. **BUT** instrumenting the
  0.5 mm/1.2 Hz case shows the wall genuinely moves (0.468 mm achieved) yet the contact
  force is nearly unchanged from static (range 0.11–0.47 N either way). Reason: the
  catheter tip (k_eff≈8.4 N/m) is ~600× softer than the tissue (5 kN/m), so the compliant
  tip **rides with** the wall — a 0.47 mm excursion that would add ~2.35 N to a rigid tip
  changes the force by only ~0.35 N here. **The cardiac motion is BUFFERED by catheter
  compliance, not stress-testing the constraint** — the opposite regime from the stiffer
  RK4 scalar plant (where the transient could breach). So this is a *refinement* of the
  RK4 Table V, not a match: the real compliant catheter is inherently more forgiving of
  cardiac motion than the reduced model suggests. (Earlier "slightly stronger than RK4 /
  mirrors RK4" framing was wrong and has been corrected here, in the script, and in §VI-G.)

**Tendon-friction hysteresis vs offset-free — DONE (the original fidelity motivation):**
Real tendon-sheath friction is Coulomb (sign-dependent → hysteretic), not the scalar
model's constant 0.45 N bias. Added MuJoCo tendon `frictionloss`; controller calibrated on
the nominal (frictionless) plant, run on the friction plant. Free-space sinusoid tracking
(0.5–3.5 mm @ 0.5 Hz, many reversals):
| plant / controller | track RMS |
|---|---|
| no friction, impedance only | 0.273 mm |
| no friction, **+ offset-free** | **0.011 mm** (96 % cut — reproduces the paper's "98 % free-space" claim) |
| friction 0.3 N, impedance only | 0.806 mm |
| friction 0.3 N, **+ offset-free** | **0.468 mm** (only 42 % cut) |
**Finding:** offset-free's near-perfect rejection holds for a *constant* bias but degrades
to ~42 % under realistic *hysteretic* friction — the integral `d̂` lags each reversal and
cannot cancel the sign-dependent part. The scalar benchmark's constant-`F_fric` hid this.
Motivates a friction/hysteresis-aware observer (a concrete new future-work item with
teeth, beyond the paper's generic §V-C remark).

**Implication for the paper:** the qualitative hard-constraint story **survives on the
physics plant**, but direct number-for-number reuse does not — the `J_κ` unit issue, the
elasticity-FF dependence, and the offset-free hysteresis degradation mean the MuJoCo port
is a *substantive re-derivation*, not a re-host. Real value (strengthens the paper) **and**
real risk (refines/conditions several current claims). Budget accordingly.

**Port status: all originally-scoped items complete** — plant built+calibrated (contact
5 kN/m), J_κ remap, MPC ported with Λ-rescaled gains, Tables II & V reproduced, tendon
hysteresis tested. Files: `catheter_mujoco.py` (plant/spike/calibration),
`catheter_mpc_mujoco.py` (MPC port + Tables II/V + friction). RK4 benchmark untouched.

## 7. Recommendation
Do it **only if** the target venue is one that will reject on "RK4 is not a physics
engine," **or** if you want the tendon-hysteresis offset-free result as a genuine new
contribution. For an RA-L / robotics-conference submission with cardiac+hardware already
scoped as future work, the RK4 benchmark is defensible and MuJoCo is a *strengthening*
exercise, not a *blocker*. If you green-light it, start with step 1 (model + calibration)
as a spike — if the contact/tendon tuning matches the scalar params cleanly, the rest is
mechanical; if it fights, re-evaluate before sinking the full 6 days.
