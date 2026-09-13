# Peer-Review and Reproducibility Audit

**Manuscript:** [paper/catheter_ieee.tex](paper/catheter_ieee.tex) with [paper/body.tex](paper/body.tex)  
**Review date:** 2026-09-11  
**Recommendation:** Rebuild and source-verify before submission

## Revision Status -- 2026-09-11

The main claim-alignment and wording issues identified below have been addressed in
the maintained TeX source:

- the title now describes force-limited control rather than claiming a validated
  predictive safety controller;
- the derived horizon QP is separated from the executable MuJoCo benchmark,
  which uses static impedance feedback, integral residual compensation,
  pointwise corrective-force clipping, and tendon saturation;
- the corrective-force proxy is distinguished from measured Kelvin--Voigt tissue
  reaction;
- the kinematic normal Jacobian and calibrated tendon-to-tip-force transmission
  now use separate symbols;
- unsupported recursive-feasibility and constrained-ISS theorems were replaced by
  the verified unconstrained DARE equilibrium and fixed-gain LPV analysis;
- solver timing is machine-qualified, and the approach sweep now reports the
  observed transition from 0.485 N at 12 mm/s to 0.591 N at 15 mm/s;
- the categorical literature capability table was removed;
- the benchmark and generated figures now use implementation-accurate controller
  names;
- avoidably self-critical phrases were recast as validation scope and next-stage
  work without removing safety-relevant qualifications;
- MuJoCo 3.13 packed-mass compatibility and an exact Python dependency manifest
  were added.

Submission still requires rebuilding both PDFs, checking the page limit, and
verifying clinical quantities and literature-comparison values against primary
sources. If predictive constrained control is to remain a central experimental
claim, the online horizon QP, Kalman estimator, tendon/curvature rows, and
contact-output model must be implemented and evaluated on the distributed plant.

## Executive Summary

The project contains a useful, reproducible simulation result. On an eight-link
MuJoCo catheter with analytically applied Kelvin--Voigt tissue reaction, integral
residual compensation reduces free-space approach RMS from 0.284 to 0.028 mm.
For a penetrating reference, adding a pointwise 0.5 N corrective-force limit
reduces measured peak tissue force from 0.595 to 0.467 N without changing the
approach RMS. The saved output and a fresh run agree.

The original manuscript attributed these results to a finite-horizon QP and an
augmented Kalman filter. The benchmark instead computed
`Kd*e + Dd*de + dhat`, clipped that scalar command, and mapped it to tendon
tension. It did not execute a horizon QP, Kalman update, online gain schedule,
curvature row, or horizon tendon row. A separate analytical script does solve a
reduced-model OSQP problem, but that is not the controller used for the headline
MuJoCo table.

The original force language was also too strong. Clipping corrective force does
not bound the complete Kelvin--Voigt tissue reaction under arbitrary approach
speed, contact-normal error, or inter-sample impact. The supplied tests themselves
show this: measured force is 0.485 N at 12 mm/s, 0.591 N at 15 mm/s, and 0.840 N
at 20 mm/s. With 30-degree normal error, the nominal cap gives a 0.539 N realized
estimate; cosine tightening reduces it to 0.458 N.

## Major Findings

### 1. The headline MuJoCo controller was not the stated horizon QP

The benchmark controller in
[simulation/catheter_benchmark_mujoco.py#L50](simulation/catheter_benchmark_mujoco.py#L50)
uses impedance gains plus an integral residual state. Its constrained branch
applies `np.clip` to the scalar corrective force before tendon mapping. The
headline loop contains no horizon construction, OSQP call, Kalman filter, or
curvature constraint.

The actual finite-horizon matrices and OSQP solve exist only in
[simulation/catheter_verify.py#L34](simulation/catheter_verify.py#L34) and
[simulation/catheter_verify.py#L85](simulation/catheter_verify.py#L85). That
script verifies a reduced two-state model with a corrective-force box; it does
not drive the MuJoCo plant in the headline comparison.

**Correction applied:** [paper/body.tex#L260](paper/body.tex#L260) now identifies
the benchmark as the static impedance specialization and lists the horizon QP,
Kalman filter, gain schedule, and curvature rows as unevaluated implementation
stages. Table and figure labels were updated accordingly.

### 2. The corrective-force clip was not a tissue-force certificate

The implementation clips `F_mpc`, then maps feedforward plus corrective force to
tendon tension. Tissue reaction is independently generated from penetration and
relative velocity in
[simulation/catheter_mujoco.py#L157](simulation/catheter_mujoco.py#L157).
Consequently, measured tissue force also depends on elastic feedforward,
contact damping, approach speed, normal alignment, actuator saturation, and the
sampled trajectory.

The fresh approach sweep confirms the distinction: the 0.5 N corrective-force
clip coexists with measured peaks of 0.591 N at 15 mm/s and 0.840 N at 20 mm/s.
The normal-misalignment sweep also exceeds 0.5 N without cosine tightening.

**Correction applied:** [paper/body.tex#L191](paper/body.tex#L191) defines a
corrective-force proxy and [paper/body.tex#L206](paper/body.tex#L206) explicitly
separates it from complete tissue reaction. The abstract, results, cardiac study,
and conclusion now describe measured outcomes under tested conditions.

### 3. The constrained stability and ISS claims were not established

The earlier nominal-stability theorem assumed recursive feasibility without
constructing a terminal set or invariant feasible region. The ISS theorem did
not provide observer-error dynamics, a Lyapunov argument for the constrained
closed loop, or an executable implementation of the stated horizon controller.
Citing standard MPC and offset-free results does not establish those missing
conditions for this implementation.

The closed-form Jury analysis does support the narrower fixed-gain LPV statement:
for the stated reduced model, the first instability occurs at
$\rho\approx3.39$, while the analyzed workspace has $\rho\in[0.70,1.0]$.

**Correction applied:** [paper/body.tex#L218](paper/body.tex#L218) now presents
only the verified unconstrained equilibrium and LPV margin. The additional
requirements for recursive feasibility and constrained ISS are stated at
[paper/body.tex#L224](paper/body.tex#L224).

### 4. Two different Jacobians were conflated

The manuscript originally used $J_n=dy/d\kappa$, with units of length squared,
both for operational-space inertia and for mapping tendon tension to tip force.
The benchmark value 0.087 is instead a dimensionless local slope
$\partial F_n/\partial T$, calibrated from press tests. Using one symbol for both
maps obscured the physical units and the source of the force proxy.

**Correction applied:** the tendon transmission is now denoted
$\eta_T=\partial F_n/\partial T$ in
[paper/body.tex#L191](paper/body.tex#L191), while $J_n$ remains the kinematic
normal Jacobian.

### 5. Timing and approach-speed statements were stale

Repeated runs on an Intel Xeon w7-2475X, Windows 11, Python 3.12.12, and OSQP
1.1.3 measured 1.5--5.8 microseconds for the pre-inverted matrix--vector step and
32.4--37.4 microseconds for warm-started OSQP update/solve. These are
solver-level timings, not complete controller latency.

The approach sweep was previously summarized as safe across tested speeds, but
the script tests unsafe 15 and 20 mm/s cases as well.

**Correction applied:** exact machine context and solver scope are reported at
[paper/body.tex#L370](paper/body.tex#L370); all speed/force pairs are reported at
[paper/body.tex#L372](paper/body.tex#L372).

### 6. The literature capability table overreached its evidence

The removed table assigned categorical rates, sensing requirements, hard-force
capability, and stability labels to broad method classes, including uncited
"Learning" and "Cosserat" categories. The repository records that these entries
had not been independently checked against all primary sources. Such a table is
not defensible as a substitute for a systematic review.

**Correction applied:** the table was replaced by a scoped objective-level
comparison. Hardware force-regulation RMSE and this paper's simulation peak are
explicitly described as different metrics and evidence levels.

### 7. A clean environment had no declared dependencies

Before this review the project contained no `requirements.txt`, lock file, or
`pyproject.toml`. The global Python lacked MuJoCo and OSQP; the VS Code-selected
environment contained the needed packages. The verification script also failed
against MuJoCo 3.13 because `MjData.qM` was renamed and the `mj_fullM` signature
changed.

**Correction applied:** [requirements.txt](requirements.txt) pins the verified
Python packages, and
[simulation/catheter_mujoco.py#L138](simulation/catheter_mujoco.py#L138)
supports both legacy and MuJoCo 3.13 mass-matrix APIs.

## Claims Verified by Fresh Execution

1. Main benchmark: 0.284/0.028/0.028/1.822 mm approach RMS for classical
   impedance, impedance plus integral, force-limited impedance plus integral,
   and joint-space PD.
2. Peak tissue force: 0.169, 0.595, 0.467, and 0.000 N for the same controllers.
3. Cardiac cases: 0.467 N static, 0.468 N at 0.3 mm/1 Hz, and 0.466 N at
   0.5 mm/1.2 Hz plus 0.2 mm position noise.
4. Force-regulation characterization: 141.9 mN RMSE/0.347 N peak for a static
   wall and 166.6 mN/0.434 N with cardiac motion and noise.
5. Reduced-model QP: 1.9857 N unconstrained first command versus 0.5000 N with
   the OSQP box; inactive-case infinity-norm difference $1.77\times10^{-3}$.
6. Reduced-model pole spread: $1.78\times10^{-6}$ with per-configuration DARE
   redesign versus $1.50\times10^{-4}$ with a fixed gain, an 84-fold reduction.
7. LPV boundary: $\rho^\star\approx3.39$.
8. Approach-speed sweep: threshold respected through the tested 12 mm/s case and
   exceeded at 15 and 20 mm/s.

## Self-Critical Language Pass

The maintained manuscript no longer uses phrases such as "three limitations
matter clinically," "only the proposed method," "does not support a blanket
claim," "worse than the reported hardware controllers," or "genuine cost the
assumption hid." These were replaced by neutral statements of validated scope,
measured operating boundaries, and next implementation stages. Necessary
qualifications about simulation-only evidence, force-proxy scope, approach
velocity, alignment, and hardware validation remain.

## Remaining Submission Work

1. Rebuild `catheter_ieee.pdf` and `catheter_arxiv.pdf`; the stored PDFs predate
   this revision. No LaTeX engine is installed in the current environment.
2. Check the resulting page count against the selected venue. The prior build was
   documented as nine pages, while the repository notes an eight-page RA-L limit.
3. Verify bibliographic metadata, force thresholds, tissue parameters, sensing
   rates, and external comparison values against primary sources.
4. Decide whether the paper's principal contribution is the derived predictive
   formulation or the implemented force-limited static specialization. An
   experimental predictive-control claim requires implementing and rerunning the
   horizon controller on the distributed plant.
5. Add integrated-platform or hardware evidence before using clinical-safety or
   certified-force language.

## Validation Record

- Selected interpreter:
  `C:\Users\HILDEV01\Documents\Yongyan\.venv\Scripts\python.exe`
  (Python 3.12.12).
- Packages: NumPy 2.5.3, SciPy 1.18.1, OSQP 1.1.3, MuJoCo 3.13.0,
  Matplotlib 3.10.5.
- `catheter_verify.py`: passed after the MuJoCo 3.13 compatibility correction.
- `catheter_benchmark_mujoco.py`: passed; all saved benchmark numbers reproduced.
- Regenerated paper figures are byte-identical to the validated simulation
  outputs.
- TeX source checks: balanced braces/environments, no duplicate labels, no
  missing references, and no undefined citation keys.
- TeX/PDF build: unavailable because `latexmk` and `pdflatex` are not installed.
