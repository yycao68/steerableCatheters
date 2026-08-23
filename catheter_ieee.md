# Interaction Dynamics Modeling and Predictive Control for Safe Steerable Catheter-Tissue Interaction

**Yongyan Cao and Jinshan Tang**

*Voryx Robotics LLC, San Jose, CA 95136, USA — yongyancao@gmail.com*


---

*Abstract*—Safe steerable catheter control is fundamentally a problem of **interaction dynamics**: the tip must follow a planned motion, remain compliant against moving tissue, reject friction and hysteresis, and respect a clinically meaningful never-exceed contact-force bound. We formulate catheter-tissue interaction dynamics in the scalar tip-normal coordinate of a single-segment single-tendon catheter. A partial-physics feedforward cancels only the reliable nominal bending dynamics, exposing a configuration-invariant linear interaction-dynamics model whose input gain varies through the scalar catheter inertia. A predictive optimizer then regulates this interaction state subject to hard contact-force, tendon-force, and curvature constraints. An augmented Kalman filter compresses tissue contact, friction, hysteresis, and modeling error into a sensor-free disturbance state, giving nominal offset-free regulation in free space while leaving force safety to the explicit constraint. The unconstrained and disturbance-free limit recovers classical catheter impedance as a special realization of the same interaction dynamics, rather than as the main design object. MuJoCo distributed-compliance simulation shows that the disturbance estimate reduces free-space approach error by 90% and that only the force-constrained predictive controller simultaneously tracks accurately and respects the 0.5 N safety bound (0.47 N peak); the unconstrained offset-free controller tracks similarly but violates the bound (0.60 N). These results show that offset-free motion regulation and contact-force safety are coupled interaction-dynamics objectives, and that the explicit predictive constraint is what resolves their tension under stiff tissue contact.

*Index Terms*—Steerable catheter, continuum robot, interaction dynamics, predictive control, Kalman filter, contact-force safety, offset-free regulation, cardiac ablation.

---

## I. Introduction

### A. Clinical Motivation

Steerable catheters are the primary tool for cardiac electrophysiology (EP) procedures including radiofrequency ablation, where the tip must be positioned precisely at target tissue while maintaining controlled, stable contact. The central control problem is therefore not merely tip tracking and not merely force regulation; it is the regulation of **catheter-tissue interaction dynamics**. The interaction state must encode how the tip moves relative to tissue, how persistent friction and contact forces bias that motion, and how safety limits reshape what motion is physically allowable.

Existing methods regulate these interaction dynamics through different mechanisms. Classical impedance control [1] shapes the tip port as a virtual mechanical impedance $Z(s) = M_d s^2 + D_d s + K_d$, providing passive compliance without an explicit contact model. It is an important interaction-dynamics realization, but three limitations matter clinically: **(i)** the contact-force safety bound $\|F_\text{tip}\| \le F_\text{safe}$ is not enforced as a prediction-horizon constraint; **(ii)** persistent contact force produces a steady-state tip error $e_\infty = K_d^{-1} F_\text{contact}$; and **(iii)** the controller cannot anticipate trajectory curvature, impending contact, or force saturation.

### B. Challenges Specific to Catheters

Catheters present challenges absent in rigid-body robots. **Model uncertainty:** the mechanics follow a variant of the Cosserat rod equations [3], [4], with distributed elasticity, sheath friction, tendon backlash, and hysteresis that vary widely across specimens, making full-model inversion impractical. **Cardiac motion:** the heart moves ~10–15 mm per cycle at ~1 Hz, so contact forces are quasi-periodic with known frequency but unknown phase and amplitude. **Limited sensing:** most clinical catheters carry no force sensor; tip position is available via electromagnetic (EM) tracking, fluoroscopy, or intracardiac echocardiography at moderate latency (10–50 ms). **Safety criticality:** a perforation force threshold of ~0.3–0.5 N is clinically relevant [12], making hard force-constraint enforcement—not soft penalization—essential.

### C. Related Work

**Catheter and continuum-robot control** has progressed from PID and Jacobian-based kinematic regulation [5], [6] to model-less feedback [6] and predictive interaction-dynamics optimization. Cosserat-based predictive control can model rich catheter dynamics but requires online nonlinear-program solution, limiting update rates to ~10–20 Hz. Constant-curvature kinematics [2] underlie most reduced-order interaction models. Learning-based approaches improve robustness to model uncertainty but usually do not provide formal stability or hard interaction-constraint guarantees.

**Sensor-free contact-force control** is an active and directly relevant line. Jolaei *et al.* [14] regulate the contact force of a tendon-driven ablation catheter *without a force sensor* by combining position control with a displacement-based viscoelastic contact model, reporting force-regulation RMS error of 0.03–0.05 N. Kesner and Howe [15] combine ultrasound guidance with force control on a moving cardiac target, achieving ~0.08 N RMS force tracking. These works regulate one component of the interaction dynamics—force—to a setpoint. The present framework instead treats force as a hard never-exceed interaction-dynamics constraint while regulating the tip-motion state, a complementary objective (Section VI-C). Note that [15] utilized ultrasound guidance to actively track the moving target, which explains their low force RMS.

**Virtual-impedance interaction shaping** has been demonstrated for catheter tip-force regulation with empirically tuned parameters. In the terminology of this paper, impedance is one way to realize interaction dynamics, but the configuration-dependent effective inertia and compliance are typically not addressed, giving inconsistent behavior across the workspace.

**Disturbance rejection in flexible robots** can also be interpreted as interaction-dynamics regulation. Active disturbance rejection control (ADRC) and its extended state observer [8] lump unknown dynamics into an estimated disturbance state. The connection to offset-free predictive control [9], [10] provides a formal framework for combining disturbance estimation with constraint-aware interaction optimization [7].

**Base framework.** This work builds on the interaction-dynamics framework introduced for redundant-manipulator physical human–robot interaction [13]: nonlinear robot dynamics are transformed into a configuration-invariant linear interaction model, then regulated by predictive optimization with disturbance augmentation and safety constraints. We carry that hierarchy—Interaction Dynamics → Configuration-Invariant Dynamics → Predictive Optimization—to the steerable-catheter setting, where the catheter-specific challenges are single-DOF actuation, severe model uncertainty, and hard contact-force safety.

### D. Contributions

1. **Correct single-DOF formulation.** A single-segment single-tendon catheter has one controllable degree of freedom (the curvature). We develop the controller in the scalar tip-normal coordinate, avoiding the common over-parameterization that treats the 2-D tip as independently actuable.

2. **Configuration-invariant interaction dynamics.** A partial-physics feedforward cancels only the *known* nominal bending stiffness and damping, reformulating the uncertain nonlinear catheter mechanics into a configuration-invariant linear interaction-dynamics model. The double integrator is the resulting model, not the contribution itself; the contribution is the interaction-dynamics reduction that makes prediction and constraints simple.

3. **Predictive interaction-dynamics optimization.** We establish (Theorem 1) that the unconstrained, disturbance-free realization recovers classical catheter impedance, while the constrained predictive realization adds offset-free disturbance rejection and explicit interaction-constraint enforcement.

4. **Hard interaction-dynamics constraints.** The QP enforces $\|F_\text{tip}\| \le F_\text{safe}$ as a hard interaction constraint over the horizon, realized through the transmitted tendon force together with tendon and curvature limits.

5. **Honest characterization of predictive value.** We show, by simulation, that the dominant benefit of prediction here is interaction-constraint enforcement, that offset-free rejection helps mainly in free space, and that offset-free motion regulation *without* a force constraint is unsafe under stiff contact.

---

## II. Interaction Dynamics Formulation

### A. Catheter Mechanics and Degrees of Freedom

Consider a planar single-segment tendon-actuated catheter; tip position $p \in \mathbb{R}^2$ is set by tendon displacement $u \in \mathbb{R}$ (positive pull gives curvature $\kappa > 0$). The dominant bending dynamics take the rigid-body-analogous form

$$M(\kappa)\ddot{\kappa} + C(\kappa,\dot\kappa)\dot{\kappa} + K(\kappa) = u + d_\text{cat}, \tag{1}$$

with effective bending inertia $M$, damping $C$, nonlinear (hysteretic) stiffness $K$, and a lumped term $d_\text{cat}$ collecting tissue contact, friction, backlash, and hysteresis. This equation is not used as a full Cosserat model; it is the starting point for constructing a low-dimensional interaction-dynamics state. The correspondence $\kappa\leftrightarrow q$, $K\leftrightarrow G$ bridges to the configuration-invariant interaction-dynamics framework established for redundant manipulators in physical human–robot interaction [13]; the present paper specializes that hierarchy to the single-DOF catheter setting.

**Degrees of freedom.** With a single tendon, the only controllable coordinate is $\kappa$; the tip pose $p(\kappa)$ traces a one-parameter curve, so the two tip coordinates are *not* independently actuable—only motion along the Jacobian direction $J_\kappa$ is. We therefore control a scalar coordinate (the curvature, or equivalently the tip displacement $y$ along the contact normal) and recover the 2-D pose through the kinematics of Section II-B.

*Remark 1 (Disturbance compression).* The full mechanics are an infinite-dimensional Cosserat PDE; (1) is a deliberate reduction in which we control only the tip state $[\kappa,\dot\kappa]$ and absorb higher-order distributed effects into $d_\text{cat}$. Not every phenomenon must be modeled—only observed and slowly varying relative to the control bandwidth.

### B. Tip Kinematics

For a constant-curvature arc of length $L$,

$$p_x = \frac{\sin(\kappa L)}{\kappa},\qquad p_z = \frac{1-\cos(\kappa L)}{\kappa}, \tag{2}$$

with continuous limits $p_x\to L$, $p_z\to0$ as $\kappa\to0$. Let $J_\kappa(\kappa)=dp/d\kappa\in\mathbb{R}^{2\times1}$ be the translational Jacobian and let $n$ be the controlled tip-normal direction. The scalar normal coordinate is $y=n^\top p$, with

$$J_n(\kappa)=\frac{dy}{d\kappa}=n^\top J_\kappa(\kappa),\qquad
\Lambda_n(\kappa)=\left(J_n(\kappa)M^{-1}(\kappa)J_n(\kappa)\right)^{-1}>0. \tag{2a}$$

This scalar operational-space inertia is defined only away from configurations where $J_n=0$ and normalizes the effective tip impedance in the controllable direction.

### C. Contact Model and the Force-Safety Interaction Constraint

During contact, a tip force $F_\text{tip}$ acts at the tip; in simulation the tissue is Kelvin–Voigt, $F_{\text{tissue}} = k_t\delta + b_t\dot\delta$, with penetration $\delta=\max(0,p_\text{surf}-p)$, tissue stiffness $k_t\approx 2$–$20$ kN/m, and damping $b_t$. The contact maps to curvature through virtual work, $\tau_\text{ext}=J_\kappa^\top F_\text{tip}$; along the controlled normal this reduces to $\tau_\text{ext}=J_n F_n$, where $F_n=n^\top F_\text{tip}$. The clinically relevant perforation threshold $F_\text{safe}\approx 0.5$ N [12] motivates a hard predicted normal-force bound, $|F_n(t)|\le F_\text{safe}$, with additional approach-velocity limits needed for environment-induced impact transients.

### D. Interaction-Dynamics Objective

Given a reference $(p_d,\dot p_d,\ddot p_d)$, define the scalar interaction state $x=[e,\dot e]^\top$, where $e=y_d-y$ is the error in the controllable normal coordinate. This state characterizes the regulated catheter-tissue interaction dynamics: motion error, velocity error, persistent unknown loading, and the constraints that limit allowable contact. The objective is to choose $u(t)$ so $e(t)\to0$ in the nominal free-space limit while respecting tendon limits $u_\text{lo}\le u\le u_\text{hi}$, predicted normal-force safety $|F_n|\le F_\text{safe}$, and curvature limits $\kappa_\text{lo}\le\kappa\le\kappa_\text{max}$ under bounded unknown contact, friction, and hysteresis.

---

## III. Predictive Interaction Dynamics Control

### A. Two-Layer Architecture

$$u = \underbrace{u_\text{ff}}_{\text{Layer 1: interaction-dynamics normalization}} + \underbrace{J_n F_\text{mpc}}_{\text{Layer 2: predictive correction}}, \tag{3}$$

with feedforward $u_\text{ff} = \hat C(\kappa,\dot\kappa)\dot\kappa + \hat K(\kappa) + J_n\Lambda_n(\kappa)\ddot y_d$ canceling the *nominal* dynamics. In the single-DOF setting $F_\text{mpc}\in\mathbb{R}$ is the scalar corrective force along the controlled tip-normal; the transmitted generalized input is $J_nF_\text{mpc}$, dimensionally consistent with $u$. For a rigid robot Layer 1 is exact; for a catheter $\hat C,\hat K$ are uncertain, so the cancellation is partial and the residual is absorbed by the disturbance estimate—more robust than inverting an unreliable model.

### B. Configuration-Invariant Interaction Dynamics

Since the system is single-DOF (Section II-A), the error $e=y_d-y$ and the corrective input $F_\text{mpc}\in\mathbb{R}$ are scalar. With $x_e=[e,\dot e]^\top\in\mathbb{R}^2$,

$$\dot x_e = \underbrace{\begin{bmatrix}0&1\\0&0\end{bmatrix}}_{A_c\,(\text{const})} x_e + \underbrace{\begin{bmatrix}0\\-\Lambda_n^{-1}(\kappa)\end{bmatrix}}_{B_c(\kappa)} F_\text{mpc} + \underbrace{\begin{bmatrix}0\\1\end{bmatrix}}_{E_c\,(\text{const})} d(t), \tag{4}$$

where $d$ lumps modeling error, contact, and the kinematic coupling $\dot J_n\dot\kappa$, normalized to acceleration units. This is the key configuration-invariant interaction-dynamics model: $A_c$ is **constant**, while only the scalar gain $B_c(\kappa)$ varies through $\Lambda_n(\kappa)$. The double integrator is the resulting normalized interaction model, not the main claim. The same constant-$A_d$, parameter-varying-$B_d$ structure exploited in the base framework [13] and, earlier, in min–max MPC under input saturation [11], enables offline precomputation of the prediction matrices. Because $A_c$ is nilpotent, the exact ZOH discretization is closed-form:

$$A_d=\begin{bmatrix}1&\Delta t\\0&1\end{bmatrix},\quad B_d(\kappa_k)=\begin{bmatrix}-\tfrac12\Lambda_n^{-1}(\kappa_k)\Delta t^2\\-\Lambda_n^{-1}(\kappa_k)\Delta t\end{bmatrix}. \tag{5}$$

### C. Disturbance Augmentation

Augment with a scalar integrating disturbance $\hat d\in\mathbb{R}$:

$$\begin{bmatrix}x_e(k+1)\\\hat d(k+1)\end{bmatrix}=\begin{bmatrix}A_d&G_d\\0&1\end{bmatrix}\begin{bmatrix}x_e(k)\\\hat d(k)\end{bmatrix}+\begin{bmatrix}B_d(\kappa_k)\\0\end{bmatrix}F_\text{mpc}(k), \tag{6}$$

with $G_d=[\tfrac12\Delta t^2,\Delta t]^\top$ the ZOH of $E_c$—**distinct from** $B_d$, which carries the $-\Lambda_n^{-1}$ input gain ($d$ is in acceleration units, so $G_d$ is not $\Lambda_n^{-1}$-scaled). A steady-state Kalman filter estimates $\hat d$ from tip-error measurements; the pair is observable with $C_\text{aug}=[I_2,0]$ since $G_d\neq 0$. No force sensor is required.

### D. Predictive Interaction-Dynamics Optimization

For a frozen or predicted inertia sequence, let $U_d(\hat d)=[\Lambda_{n,0}\hat d,\ldots,\Lambda_{n,N-1}\hat d]^\top$ be the steady force sequence that cancels the estimated acceleration disturbance. The condensed input-centered QP is

$$\min_{U}\ \tfrac12 U^\top H\,U
+ \Big(\Gamma^\top\bar Q\big(\Phi x_e(k)+\Delta(\hat d)\big)-\bar R U_d(\hat d)\Big)^\top U, \tag{7}$$

with $U=[F_\text{mpc}(0);\ldots;F_\text{mpc}(N-1)]\in\mathbb{R}^N$, $H=\Gamma^\top\bar Q\Gamma+\bar R$, $Q=\text{diag}(K_d,D_d)$, and $\Phi\in\mathbb{R}^{2N\times2}$, $\Gamma$ precomputed once (constant $A_d$). When no interaction constraint is active the solution is the closed form

$$U^\star=-H^{-1}\!\left[\Gamma^\top\bar Q\big(\Phi x_e+\Delta(\hat d)\big)-\bar R U_d(\hat d)\right], \tag{8}$$

a matrix-vector multiply. With constant $\Lambda_n$ over the horizon, $\Delta(\hat d)+\Gamma U_d(\hat d)=0$, so $V=U-U_d$ reduces the optimizer to the nominal interaction-dynamics regulator in $V$. When a constraint binds, OSQP solves the condensed QP reusing the cached factorization. The optimization is subject to tendon limits, curvature limits, and the contact-force interaction bound.

**Force-safety interaction constraint.** The corrective input $F_{\text{mpc}}$ is the tip-normal force, and the two-layer command is realized as the tendon tension $T_k=-(u_{\text{ff},k}+F_{\text{mpc},k})/J_n(\kappa_k)$ with the tip-space elastic feedforward $u_{\text{ff}}=k_{\text{eff}}\,y_d$ ($J_n$ bounded away from zero). The predicted normal contact force—the applied tip force beyond the catheter's own elastic restoring—is then

$$\hat F_{n,k}=J_n(\kappa_k)\,|T_k|-k_{\text{eff}}\,|y_k|=k_{\text{eff}}\,e_k+F_{\text{mpc},k}, \tag{9}$$

with any calibrated free-space friction bias absorbed into the disturbance estimate $\hat d$ (Disturbance Compression Principle). Because the compliant elastic term $k_{\text{eff}}e_k$ is small over the workspace ($|k_{\text{eff}}e_k|\ll F_{\text{safe}}$), the hard predicted-force bound $|\hat F_{n,k}|\le F_{\text{safe}}$ is enforced to leading order by the box constraint on the corrective force,

$$-F_{\text{safe}}\le F_{\text{mpc},k}\le F_{\text{safe}}, \tag{10}$$

applied directly in the QP. This bounds the controller-induced/predicted quasi-static normal force; fast environment-induced damping spikes require approach-velocity limits.

### E. LQR-Impedance Equivalence

**Theorem 1 (LQR-impedance equivalence).** *In the unconstrained, disturbance-free limit ($d\equiv0$, no active constraints) the receding-horizon law is the static linear state feedback $F_\text{mpc}=K_\text{eff}e+D_\text{eff}\dot e$—the LQR-realized instance of the classical catheter impedance law $u_\text{imp}=\hat C\dot\kappa+\hat K+J_n(K_\text{eff}e+D_\text{eff}\dot e)$, realizing the effective tip-normal impedance $Z_\text{eff}(s,\kappa)\approx\Lambda_n(\kappa)s^2+D_\text{eff}s+K_\text{eff}$—where the realized gains $(K_\text{eff},D_\text{eff})$ are the unconstrained LQR feedback for the weights $(Q,R)$. The unconstrained predictive interaction controller is therefore an LQR-tuned classical impedance, not a free rendering of a prescribed $(K_d,D_d)$ (following the pHRI base, Theorem 1).*

*Proof sketch.* Without constraints and $d\equiv0$ the QP minimizer is the stationary point $U^\star=-H^{-1}\Gamma^\top\bar Q\Phi x_e$, whose first block is a static gain $F_\text{mpc}=K_\text{eff}e+D_\text{eff}\dot e$; with the DARE terminal cost this is the infinite-horizon LQR feedback for $(Q,R)$, independent of $N$. After Layer-1 cancellation $\ddot e=-\Lambda_n^{-1}(\kappa)F_\text{mpc}+d$; multiplying by $\Lambda_n(\kappa)$ gives $\Lambda_n\ddot e+D_\text{eff}\dot e+K_\text{eff}e=\Lambda_n d$, the stated impedance, with $u_\text{imp}$ recovered after adding $u_\text{ff}$. $\square$

**Remark (Prescribed vs. realized gains).** The cost weights $Q=\text{diag}(K_d,D_d)$, $R$ are design *penalties*; the realized impedance $(K_\text{eff},D_\text{eff})$ is their LQR (Riccati) image, so in general $(K_\text{eff},D_\text{eff})\neq(K_d,D_d)$—in particular the cheap-control limit $R\to0$ yields a plant-determined gain, not $(K_d,D_d)$. To render a specific $(K_d,D_d)$ exactly, prescribe the impedance gain directly ($F_\text{mpc}=K_d e+D_d\dot e$), making the equivalence exact at the cost of the predictive look-ahead (cf. pHRI base, Remark 2).

When constraints are active or $d\neq0$, the predictive interaction controller departs from the static impedance law—precisely its advantage.

---

## IV. Stability Analysis

**Theorem 2 (Nominal stability).** *With terminal cost $Q_f$ the DARE solution at a nominal configuration $\kappa_0$, if the QP is feasible at $k=0$ and the LPV variation $\|B_d(\kappa_k)-B_d(\kappa_0)\|$ is sufficiently small—quantified for this scalar plant by $\rho=\Lambda_{n,\mathrm{ref}}/\Lambda_{n,\mathrm{true}}<\rho^\star\approx3.4$, comfortably satisfied by the $1.4\times$ workspace variation—the closed loop is asymptotically stable and $x_e(k)\to0$ (no disturbance), provided recursive feasibility holds.*

**Theorem 3 (ISS / offset-free).** *If $\|d(t)\|\le\bar d$, the augmented system is observable, and the input-centered QP is feasible, $(x_e,\hat d)$ is input-to-state stable with an error bound that scales with the unmodeled time variation of $d$. For constant matched disturbances in the nominal, free-space limit, $\hat d\to d_\infty$ implies $U_d=\mathbf{1}_N\Lambda_n d_\infty$ and $\Delta(d_\infty)+\Gamma U_d=0$, hence the unconstrained law becomes $F_\text{mpc}=\Lambda_n d_\infty+K_\text{eff}e+D_\text{eff}\dot e$ and the steady error is exactly zero. For slowly varying disturbances the residual error scales with $\|\dot d\|$. Under stiff contact the achievable steady-state error is additionally bounded by the force constraint and is therefore nonzero (Section VI-B confirms a contact-limited residual).*

*Proof status.* Theorem 1 is supported by the sketch above. Theorems 2–3 follow standard arguments—terminal-cost/recursive-feasibility stability for constrained MPC [7] and the augmented-observer offset-free result [9]—specialized to the constant-$A_d$ scalar plant with the input-centering identity $\Delta(\hat d)+\Gamma U_d(\hat d)=0$.

*Remark 2 (Quantified LPV margin).* The "sufficiently small variation" condition of Theorem 2 admits a concrete bound for this scalar plant. Since $B_d(\kappa)\propto\Lambda_n^{-1}(\kappa)$, a gain $K$ designed at $\Lambda_{n,\text{ref}}$ applied where the true inertia is $\Lambda_{n,\text{true}}$ yields the closed loop $A_d-\rho\,B_d(\Lambda_{n,\text{ref}})K$ with $\rho=\Lambda_{n,\text{ref}}/\Lambda_{n,\text{true}}$. Computing the spectral radius (`catheter_verify.py`), the loop remains Schur for $\rho<3.4$, i.e. for $\Lambda_{n,\text{true}}>0.29\,\Lambda_{n,\text{ref}}$; under-estimates ($\Lambda_{n,\text{true}}>\Lambda_{n,\text{ref}}$) are always stable. Over the workspace $\Lambda_n(\kappa)$ varies only $1.4\times$ ($\rho\in[0.70,1.0]$)—comfortably inside the stable region—so a single fixed gain is robustly stable, and the $\Lambda_n(\kappa)$-normalized gain is configuration-invariant to machine precision (Section VI-F). The nominal-limit caveat on Theorem 3 is essential: as the simulation shows (Sections VI-B, VI-E), offset-free tracking does **not** hold under stiff contact or fast cardiac motion.

---

## V. Discussion

### A. The Role of Predictive Interaction Optimization

| Prediction source | Quality | Benefit |
|---|---|---|
| Force prediction $\hat F(k+i)$ | poor (tissue unknown) | limited |
| Trajectory feedforward $\ddot p_d(k+i)$ | exact (offline plan) | **high** |
| Constraint anticipation | exact (known bounds) | **high** |
| Kalman disturbance estimate at $k$ | moderate | **high** |

The primary value of prediction here is **interaction-constraint enforcement** and **trajectory feedforward**, not contact-force forecasting; a zero-order-hold $\hat d(k+i)=\hat d(k)$ is adopted, and the offset-free guarantee is preserved at steady state regardless of prediction quality.

### B. The $\dot d=0$ Assumption and Cardiac Motion

The integrating model is a worst-case approximation guaranteeing offset-free tracking for any slowly varying disturbance. For the ~1 Hz quasi-periodic cardiac force, two refinements are possible: (A) increase the disturbance process-noise covariance for faster adaptation (trading noise sensitivity); (B) with ECG-gated phase, a periodic disturbance model $\hat d(k+1)=\hat d_\text{DC}(k)+A_\text{card}\sin(\omega_\text{heart}(k{+}1)\Delta t+\phi)$ turns cardiac force into a predictable component—left to future work. During first contact the estimate lags by several samples; the hard force constraint provides the safety backstop in this transient.

### C. Partial vs. Full Physics Cancellation

Catheter stiffness $K(\kappa)$ is subject to hysteresis, temperature dependence, and specimen variability; inverting with an inaccurate $\hat K$ can *increase* the effective disturbance. Canceling only reliably estimated components is more robust—the residual simply adds to $d(t)$. Design rule: include a feedforward term only if its omission would push the steady disturbance beyond the Kalman tracking bandwidth.

### D. Sensing Requirements

Tip position for Kalman updates is available from EM tracking (10–40 Hz, ~1 mm), fluoroscopy (intermittent), ICE (operator-dependent), or FBG shape sensing (continuous, ~0.5 mm). A force sensor is *not* required; if FBG force estimation is available, $\hat F_\text{tip}$ can be injected directly to reduce contact-transient lag.

### E. Robustness to Contact-Normal Misalignment

The scalar coordinate $y$ is taken along an assumed contact normal. Proximal deflections along the femoral access path can rotate the true tissue normal by an angle $\theta$ relative to this assumed direction. Two effects follow, both benign for small-to-moderate $\theta$. First, the *effective* contact stiffness seen along the assumed axis scales as $k_t\cos^2\theta$: a $15^\circ$ misalignment changes it by $\cos^2 15^\circ\approx 0.93$, a $7\%$ error that is exactly the kind of slowly-varying model mismatch the disturbance state $\hat d$ absorbs (Remark 1). Second, a tangential component $\sim F\sin\theta$ appears off-axis; for the single-DOF controller this is an unobserved direction, but it does not enter the regulated coordinate and, being bounded, leaves the offset-free and ISS properties intact. For the force-safety bound the projection is *conservative in the wrong direction*: the estimated normal force underestimates the true magnitude by $1/\cos\theta$, so a misalignment-aware design should either tighten $F_\text{safe}$ by $\cos\theta_\text{max}$ or estimate $\theta$ from the shape sensor. Quantifying tracking under large, time-varying $\theta$ (and the multi-DOF extension that makes the normal observable) is future work.

---

## VI. Interaction-Dynamics Simulation Results

### A. Setup

We evaluate the controller on a *distributed-compliance* plant in MuJoCo: an eight-link pseudo-rigid-body tendon-driven catheter (total length 5 cm), actuated through a routed spatial tendon and pressing on a Kelvin-Voigt tissue wall ($k_t=5$ kN/m, $b_t=40$ N·s/m). This plant violates the reduced-order assumptions: eight elastic DOFs, configuration-dependent tendon transmission, and distributed bending that the scalar model lumps away. The scalar operational-space inertia and tendon-to-tip-normal transmission are read from the model: $\Lambda_n\approx3.5\times10^{-3}$, $J_n\approx0.087$, catheter tip stiffness $k_\text{eff}\approx8.4$ N/m. Control rate is 500 Hz; tendon limit is $\pm8$ N; $F_\text{safe}=0.5$ N. The reference approaches the surface (1 s), presses to a target 1.5 mm *past* it, holds 1 s, and retracts, deliberately stressing the safety mechanism. Reproducible: `simulation/catheter_benchmark_mujoco.py`.

### B. Controller Comparison

**TABLE I: Interaction-Dynamics Comparison Under Free Motion and Contact**

| Controller | Approach RMS (mm) | Max contact force (N) | $F_\text{safe}$ violated? | Hold pos. error (mm) |
|---|:---:|:---:|:---:|:---:|
| Classical impedance | 0.28 | 0.17 | No | 1.48 |
| Predictive ID (no FC) | 0.03 | 0.60 | **Yes** | 1.41 |
| **Predictive ID (with FC)** | **0.03** | **0.47** | **No** | 1.41 |
| Joint-space PD | 1.82 | 0.00 | No | 3.76 |

Key results:
- **Offset-free rejection of the model-reduction residual (free space).** The uncancelled higher-order bending dynamics and configuration-dependent $J_n$ act as a slowly varying residual; classical impedance leaves a 0.28 mm approach error, which the offset-free integrator drives to 0.03 mm (**90% reduction**).
- **Force safety is the decisive interaction constraint.** Only the force-constrained predictive interaction controller achieves both accurate tracking (0.03 mm) and the 0.5 N bound (peak 0.47 N). Classical impedance and joint-space PD stay under the bound only by being too soft to track the penetrating target.
- **Offset-free motion regulation *without* a force limit overshoots the bound.** The unconstrained predictive controller, driving hard to eliminate position error against the penetrating target, pushes the contact force to 0.60 N (0.595 N precisely)—19% over $F_\text{safe}$. Offset-free motion regulation and force safety remain in tension, resolved only by the hard interaction constraint.
- **Hold error is contact-limited.** The tracking controllers retain ~1.41 mm hold error because the commanded depth needs more force than the tendon can deliver while respecting the force bound. The simulation therefore does **not** support a blanket "zero steady-state error" claim under stiff contact—offset-free holds in free space and in the nominal limit.

### C. Comparison with Published Catheter-Control Approaches

Table II positions the proposed predictive interaction-dynamics framework against representative classes of catheter/continuum-robot control from the literature. Direct numerical comparison is limited—reported conditions (platform, sensing, contact, metric) differ substantially—so the table compares *capabilities* and order-of-magnitude characteristics rather than asserting head-to-head accuracy. The literature entries (e.g., the "~10–20 Hz" Cosserat-NMPC rate) are **representative characterizations and have not been independently verified** against the primary sources; only the row for this work is from our own simulation (Section VI-B).

**TABLE II: Capability Comparison with Representative Catheter-Control Methods**

| Method (class) | Model | Update rate | Hard force bound | Offset-free | Force sensor | Formal stability |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| Classical impedance [1] | port impedance | kHz | **No** | No | optional | passivity |
| Jacobian / model-less [5], [6] | kinematic | ~10–100 Hz | No | No | local | partial |
| Cosserat NMPC [4] | full Cosserat | ~10–20 Hz | soft only | varies | yes | NLP-dependent |
| Learning-based | data-driven | high | No | No | varies | none |
| ADRC/ESO [8] | lumped-disturbance | kHz | No | yes | yes (ESO) | yes |
| **Predictive interaction dynamics (this work)** | reduced + Kalman | **500 Hz** | **Yes (hard)** | **yes (nominal)** | **no** | ISS + DARE |

**Observations.** (1) Among the surveyed classes, only the proposed framework provides a *hard* contact-force interaction constraint with formal feasibility, the clinically essential capability (Section VI-B). (2) It retains a high update rate (500 Hz) because the configuration-invariant $A_d$ structure makes the unconstrained step a matrix–vector multiply, unlike Cosserat predictive control whose online NLP limits rates to ~10–20 Hz. (3) It requires no force sensor, unlike ADRC formulations that key off a measured force channel. (4) Its formal guarantees (ISS, DARE-based stability) exceed those of learning-based and purely kinematic controllers. The trade-off is reduced model fidelity relative to full Cosserat predictive control—mitigated by the disturbance-compression principle (Remark 1).

**Quantitative benchmark.** The closest published baselines are the sensor-free tendon-driven ablation-catheter force controllers of Jolaei *et al.* [14] and Kesner and Howe [15]. Table III places the present work against their reported numbers. A direct head-to-head is *not* meaningful because the objectives differ: [14], [15] **regulate** contact force to a setpoint, whereas this work tracks a position trajectory while **bounding** the contact force as a hard constraint. The comparison should therefore be read by capability, not by the force-error column alone.

**TABLE III: Quantitative Comparison with Reported Sensor-Free / Force-Controlled Catheter Results**

| Work | Objective | Sensing | Force result | Position result | Hard force bound |
|---|---|:---:|:---:|:---:|:---:|
| Jolaei *et al.* [14] | force regulation | sensor-free (contact model) | 0.03–0.05 N RMSE | — | no (setpoint) |
| Kesner & Howe [15] | force regulation | force sensor + ultrasound | ~0.08 N RMSE | — | no (setpoint) |
| **This work** (safety mode) | position tracking + force safety | sensor-free (Kalman) | bound held, peak **0.47 N** ($\le 0.5$) | 0.03 mm approach RMS | **yes (hard)** |
| **This work** (force-reg mode) | force regulation to 0.2 N + bound | sensor-free (Kalman) | ~140 mN RMSE (idealized); ~170 mN (cardiac) | — | **yes (hard)** |

**Honest reading.** Two points must be stated plainly. (1) The same formulation can also *regulate* force, and the hard safety bound coexists with a regulation setpoint in one controller; however, the position-tracking objective and the force-regulation objective compete, and the simple integral force law regulates to only ~140 mN RMSE even in idealized simulation—the force-regulation mode requires dedicated tuning or an explicit position/force task hierarchy. We therefore make **no** claim of force-regulation parity with [14], [15] (whose 30–80 mN are *hardware* results). (2) Under 1 Hz cardiac motion, force-regulation error is ~170 mN, motivating a cardiac-phase-aware periodic model. The distinct contribution is the unification of offset-free position tracking with a *hard never-exceed* force bound and feasibility guarantees in a single constrained optimization—a capability setpoint force regulation does not provide.

### D. Disturbance Estimation

The Kalman filter converges within a few sample periods; during the contact-onset transient the estimate lags, but the tendon-force cap prevents bound violation before convergence—the constraint, not the estimate, provides the safety guarantee. The disturbance estimate improves *tracking*; the hard constraint provides *safety*.

### E. Force-Regulation Mode and the Cardiac-Motion Limit

To compare like-for-like with the setpoint force controllers of [14], [15], we add a force-regulation mode: during free-space approach the controller tracks position as before, and during contact it regulates a sensor-free contact-force estimate (the commanded tip force beyond the modeled elasticity, $\hat F_c = J_n|T| - k_\text{eff}|y|$) to a clinical target of $0.2\,\text{N}$ via an integral law on the tendon command, while the hard predicted-force bound (10) remains active. Two conditions are run over the contact-hold window:

- **Idealized (slowly varying disturbance):** force RMSE **~140 mN** about the 0.2 N target. With the disturbance propagation corrected, the offset-free position drive and the force-regulation integral law compete (the controller simultaneously tries to reach the penetrating position target and hold 0.2 N), so the simple integral law no longer regulates tightly; recovering low force RMSE requires re-tuning the force-regulation gain or an explicit position/force task hierarchy (future work). The hard never-exceed bound is unaffected.

- **Realistic (1 Hz cardiac wall motion, 0.3 mm amplitude, plus 0.2 mm position noise):** force RMSE is **~170 mN**. The zero-order-hold $\dot d = 0$ model cannot anticipate the periodic motion (Section V-B), so the estimate lags. This is *worse* than the reported hardware controllers and identifies the cardiac-phase-aware periodic disturbance model (Section V-B, Option B) as **essential** future work for force-regulation use—not merely an enhancement.

In both force-regulation conditions the *peak* measured force stays under $F_\text{safe}$ (~0.35 N idealized, ~0.43 N with cardiac motion): the compliant catheter closes on the wall gently, so the regulation error shows up as a poorly-*regulated* (not unsafe) force, while the hard predicted-force bound caps the controller-induced command throughout.

The safety-mode results (Table I) are unaffected by this limitation: the hard bound concerns controller-induced force and is enforced regardless of disturbance model, whereas tight force *regulation* under cardiac motion additionally requires predicting the motion.

### F. Verification of Structural Claims

Three claims not exercised by the constant-$\Lambda$ position-tracking benchmark above are verified separately (`simulation/catheter_verify.py`):

1. **Configuration-adaptive compliance.** Computing $\Lambda_n(\kappa)$ from the scalar normal Jacobian (2a) over $\kappa\in[2,25]\,\text{m}^{-1}$ (a 1.4$\times$ inertia variation), the $\Lambda_n(\kappa)$-normalized controller holds the closed-loop pole spread to $\sim2\times10^{-6}$—roughly 80$\times$ tighter than a fixed-$\Lambda_n$ gain. This confirms that the $\Lambda_n(\kappa)$ normalization—not manual gain scheduling—is what renders the tip compliance configuration-independent. The magnitude of the fixed-gain drift is modest for a single segment's 1.4$\times$ range but grows for the larger inertia variation of multi-segment catheters (Section VII).

2. **Hard constraint via the QP.** For a 3 mm error the unconstrained step demands $|F_\text{mpc}| = 1.99\,\text{N}$; solving the condensed QP with OSQP under a 0.5 N limit returns $|F_\text{mpc}| = 0.50\,\text{N}$ (constraint respected), and reproduces the closed-form solution to $10^{-3}$ when the constraint is inactive. This verifies that the actuator/safety limits are genuinely enforced by the QP, with the closed-form step (Section III-D) recovered when no constraint binds.

3. **Real-time feasibility.** The unconstrained closed-form step runs in $\approx$0.9 $\mu$s and the warm-started OSQP solve in $\approx$27 $\mu$s per step on a desktop CPU—both far inside the 2 ms (500 Hz) budget, supporting the real-time claim with margin for embedded hardware. (Embedded-target timing remains future work.)

4. **Approach-velocity safety limit.** A sweep of constant-velocity impacts (FC controller) shows the contact force stays within $F_\text{safe}$ for approach speeds up to $\approx$7 mm/s (peak 0.49 N at 7 mm/s) and breaches it at 8 mm/s (0.54 N). This is consistent with the damping-only bound $v < F_\text{safe}/b_t = 12.5$ mm/s, the simulated limit being lower because penetration adds to the $b_t\dot y$ term. The qualification is important for the safety claim: **the hard constraint bounds *controller-induced* force, but a fast impact injects an *environment-induced* $b_t\dot y$ transient that no input limit can prevent within one sample.** Safe operation therefore requires the approach to be velocity-limited (or the contact onset to be detected and the reference ramped)—a planning constraint, not a controller failure.

---

## VII. Conclusion

Rather than viewing steerable catheter control as an impedance-control problem, this work formulates it as an interaction-dynamics regulation problem. The catheter-tissue interaction is first reduced to the correct scalar tip-normal state, then normalized into configuration-invariant linear dynamics, and finally regulated by predictive optimization with disturbance augmentation and hard interaction constraints. Classical impedance appears as the unconstrained, disturbance-free special case; the clinically relevant controller departs from it when offset-free motion regulation, tendon limits, curvature limits, and contact-force safety compete. Simulation shows that disturbance augmentation reduces free-space approach error, but the decisive result is interaction safety: only the force-constrained predictive controller keeps contact force within the 0.5 N bound, whereas an unconstrained offset-free controller drives through tissue. The resulting framework provides a unified foundation for compliant catheter motion, constraint handling, and safety-critical tissue interaction, and is extensible to multi-segment catheters, surgical robotics, rehabilitation, whole-body control, and dexterous manipulation.

**Future work:** (i) hardware validation on a tendon-actuated catheter with EM tracking; (ii) FBG-based force estimation as a direct disturbance channel; (iii) cardiac-phase-aware periodic disturbance models; (iv) extension to multi-segment catheters with $m$ independent tendons (an $m$-DOF version with configuration-varying $\Lambda_n(\kappa)$); (v) energy-budget augmentation toward certified passivity during aggressive contact; and (vi) a hysteresis-aware disturbance model for sign-dependent tendon-sheath friction.

---

## References

[1] N. Hogan, "Impedance control: An approach to manipulation, Parts I–III," *ASME J. Dyn. Syst. Meas. Control*, vol. 107, no. 1, pp. 1–24, 1985.

[2] R. J. Webster III and B. A. Jones, "Design and kinematic modeling of constant curvature continuum robots: A review," *Int. J. Robot. Res.*, vol. 29, no. 13, pp. 1661–1683, 2010.

[3] S. S. Antman, *Nonlinear Problems of Elasticity*, 2nd ed. New York: Springer, 1995.

[4] D. C. Rucker and R. J. Webster III, "Statics and dynamics of continuum robots with general tendon routing and external loading," *IEEE Trans. Robot.*, vol. 27, no. 6, pp. 1033–1044, 2011.

[5] D. B. Camarillo, C. F. Milne, C. R. Carlson, M. R. Zinn, and J. K. Salisbury, "Mechanics modeling of tendon-driven continuum manipulators," *IEEE Trans. Robot.*, vol. 24, no. 6, pp. 1262–1273, 2008.

[6] M. C. Yip and D. B. Camarillo, "Model-less feedback control of continuum manipulators in constrained environments," *IEEE Trans. Robot.*, vol. 30, no. 4, pp. 880–889, 2014.

[7] J. B. Rawlings, D. Q. Mayne, and M. Diehl, *Model Predictive Control: Theory, Computation, and Design*, 2nd ed. Madison, WI: Nob Hill, 2017.

[8] J. Han, "From PID to active disturbance rejection control," *IEEE Trans. Ind. Electron.*, vol. 56, no. 3, pp. 900–906, 2009.

[9] G. Pannocchia and J. B. Rawlings, "Disturbance models for offset-free model-predictive control," *AIChE J.*, vol. 49, no. 2, pp. 426–437, 2003.

[10] U. Maeder, F. Borrelli, and M. Morari, "Linear offset-free model predictive control," *Automatica*, vol. 45, no. 10, pp. 2214–2222, 2009.

[11] Y.-Y. Cao and Z. Lin, "Min–max MPC algorithm for LPV systems subject to input saturation," *IEE Proc. Control Theory Appl.*, vol. 152, no. 3, pp. 266–272, 2005.

[12] H. Yokoyama, H. Nakagawa, J. Kautzner, et al., "Novel contact force sensor incorporated in irrigated radiofrequency ablation catheter predicts lesion size and incidence of steam pop and thrombus," *Circ. Arrhythm. Electrophysiol.*, vol. 1, no. 5, pp. 354–362, 2008.

[13] Y. Cao and J. Tang, "Impedance MPC for physical human–robot interaction: Predictive disturbance rejection with joint-limit safety," *IEEE Trans. Robot.*, under review, 2026. (Base framework; companion paper.)

[14] M. Jolaei, A. Hooshiar, A. Sayadi, J. Dargahi, and M. Packirisamy, "Sensor-free force control of tendon-driven ablation catheters through position control and contact modeling," in *Proc. 42nd Annu. Int. Conf. IEEE Eng. Med. Biol. Soc. (EMBC)*, 2020, pp. 5248–5251, doi: 10.1109/EMBC44109.2020.9176019.

[15] S. B. Kesner and R. D. Howe, "Robotic catheter cardiac ablation combining ultrasound guidance and force control," *Int. J. Robot. Res.*, vol. 33, no. 4, pp. 631–644, 2014, doi: 10.1177/0278364913511350.

---

*Note on verification status. Reference [14] was verified against its primary source; [1]–[12] are standard works whose bibliographic details should be checked before submission, and [13] is the under-review companion. Clinical figures (perforation threshold, cardiac-motion amplitude, tissue stiffness, sensing rates) and the Table II literature characterizations are representative and not independently verified here. Simulation results are produced by `simulation/catheter_benchmark.py` (Sections VI-B/E) and `simulation/catheter_verify.py` (Section VI-F), both single-DOF, NumPy/SciPy/OSQP.*
