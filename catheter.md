# Impedance MPC for Steerable Catheter Tip Control Under Unknown Contact Forces

**Draft — For Discussion**

---

## Abstract

Steerable catheters used in cardiac electrophysiology and endovascular interventions must simultaneously track a planned trajectory and maintain safe, compliant contact with tissue — two objectives that classical impedance control cannot fulfill without steady-state error under persistent contact forces, and without any formal mechanism for enforcing safety constraints. We present an Impedance MPC framework for catheter tip control that generalizes classical impedance control to the catheter setting: a feedforward partial-physics cancellation layer reduces the uncertain, nonlinear catheter dynamics to a linear residual plant, and a receding-horizon QP optimizes over this residual subject to hard contact-force safety constraints and tip curvature limits. An augmented Kalman filter estimates the lumped disturbance — comprising tissue contact forces, friction, hysteresis, and modeling error — and renders the tracking offset-free in the nominal limit. We establish a formal equivalence between the proposed MPC and classical catheter impedance control in the unconstrained, disturbance-free case, and prove input-to-state stability under bounded disturbances. The effective tip impedance is automatically normalized by the configuration-dependent operational-space inertia $\Lambda(\kappa)$, providing curvature-adaptive compliance without manual gain scheduling. Because a single-tendon single-segment catheter is one controllable degree of freedom, we develop the controller in the scalar tip-normal coordinate. Simulation on a single-DOF catheter model shows that the Kalman augmentation rejects a constant friction/hysteresis bias, cutting free-space approach error by 40% relative to classical impedance, and — the decisive result — that **only the force-constrained MPC keeps the contact force within the 0.5 N safety bound** (peak 0.36 N vs. 0.79 N for classical impedance). We further show that offset-free disturbance rejection and contact-force safety are in direct tension: an unconstrained offset-free controller, driving to eliminate position error against stiff tissue, pushes the contact force to over 2 N. The explicit hard force constraint is what resolves this — a capability no impedance controller possesses.

---

## 1. Introduction

### 1.1 Clinical Motivation

Steerable catheters are the primary tool for cardiac electrophysiology (EP) procedures including radiofrequency ablation, where the catheter tip must be positioned precisely at target tissue sites while maintaining controlled, stable contact. Two competing objectives define the control problem:

1. **Trajectory tracking**: the tip must reach and hold a planned position with sub-millimeter accuracy to ensure ablation lesions are placed correctly.
2. **Safe compliance**: the tip must not exert excessive force against the endocardium, as over-contact causes perforation (a potentially fatal complication), while under-contact leads to inadequate lesion formation.

Classical impedance control [HOGAN1985] shapes the port behavior of the tip as a mechanical impedance $Z(s) = M_d s^2 + D_d s + K_d$, providing passive compliance without requiring an explicit contact model. It has been applied to catheter tip force regulation [REF] because of its simplicity and inherent passivity. However, it suffers from three fundamental limitations in the catheter setting:

- **(i) No constraint enforcement**: the contact force safety bound $\|F_\text{tip}\| \leq F_\text{safe}$ cannot be enforced; post-hoc force clipping destroys passivity.
- **(ii) Steady-state position error**: any persistent contact force produces a nonzero tip position error $e_\infty = K_d^{-1} F_\text{contact}$, preventing accurate lesion placement.
- **(iii) No look-ahead**: the controller cannot anticipate trajectory curvature or predict approaching tissue contact.

### 1.2 Challenges Specific to Catheters

Beyond the limitations shared with any impedance controller, catheters present additional challenges absent in rigid-body robots:

**Model uncertainty.** The catheter mechanics are governed by a variant of the Cosserat rod equations [ANTMAN1995], which involve distributed elasticity, friction at the sheath interface, tendon backlash, and hysteresis. These phenomena are highly variable across catheter specimens and operating conditions, making full-model inversion (as used in computed-torque control for rigid robots) impractical.

**Cardiac motion.** The heart moves ~10–15 mm per cycle at ~1 Hz. Contact forces therefore oscillate periodically, creating a disturbance that is neither constant nor unpredictable — it is quasi-periodic with known frequency but unknown phase and amplitude.

**Limited sensing.** Most clinical catheters do not carry force sensors. Position sensing is typically available via electromagnetic (EM) tracking, fluoroscopy, or intracardiac echocardiography (ICE), providing tip position at moderate latency (10–50 ms). Fiber Bragg Grating (FBG) sensors can provide shape and strain information, enabling force estimation, but are not yet clinically standard.

**Safety criticality.** Unlike industrial robots, errors in the catheter setting have immediate patient safety implications. A perforation force threshold of approximately 0.3–0.5 N is clinically relevant [REF], making hard force constraint enforcement — not just soft penalization — essential.

### 1.3 Related Work

**Catheter control.** Model-based catheter control has progressed from PID [REF] to Jacobian-based controllers [REF] and more recently to MPC formulations [REF]. Cosserat-based model predictive control [REF] achieves high tip accuracy but requires online solution of a nonlinear program (NLP), limiting update rates to 10–20 Hz. Learning-based approaches [REF] improve on model uncertainty but lack formal stability and constraint guarantees.

**Impedance control for catheters.** Classical impedance control for catheter tip force regulation has been demonstrated in [REF], where the impedance parameters are tuned empirically. The configuration-dependent effective impedance is not addressed, leading to inconsistent compliance across the workspace.

**Disturbance rejection in flexible robots.** Active disturbance rejection control (ADRC) [HAN2009] and its extended state observer (ESO) share the core idea of lumping unknown dynamics into an estimated disturbance state. The connection to MPC-based offset-free tracking [PANNOCCHIA2003, MAEDER2010] provides a formal framework for combining disturbance estimation with constraint-aware optimization.

### 1.4 Contributions

This paper makes the following contributions:

1. **Partial-physics feedforward for catheters.** We formulate a feedforward cancellation layer that removes only the *known* components of catheter dynamics (a nominal bending stiffness and damping), leaving the residual plant as a double integrator in tip-position error space. Unlike full computed-torque cancellation for rigid robots, this partial cancellation is robust to the large modeling uncertainty inherent in catheter mechanics.

2. **Impedance MPC unification.** We establish that Impedance MPC with an impedance-shaped cost matrix recovers the classical catheter impedance law in the unconstrained, disturbance-free case, while providing constraint enforcement and offset-free disturbance rejection beyond this case.

3. **Hard contact-force safety constraints.** The QP enforces $\|F_\text{tip}\| \leq F_\text{safe}$ as a hard constraint over the prediction horizon — a capability unavailable to any impedance controller.

4. **Disturbance compression.** Tissue contact forces, friction, hysteresis, and modeling error are unified into a single augmented disturbance state estimated by a Kalman filter, eliminating the need for explicit modeling of each phenomenon.

5. **Analysis of prediction quality.** We explicitly characterize when MPC look-ahead provides benefit over a pure reactive controller in the catheter setting, and show that the primary value of MPC here lies in constraint enforcement and trajectory feedforward rather than force prediction accuracy.

---

## 2. System Model and Problem Formulation

### 2.1 Catheter Mechanics

We consider a single-segment tendon-actuated steerable catheter operating in the plane. The catheter tip position $p \in \mathbb{R}^2$ is controlled via tendon displacement $u \in \mathbb{R}$ (positive pull produces positive curvature $\kappa > 0$).

**Degrees of freedom.** A single-segment, single-tendon catheter has *one* controllable degree of freedom: the curvature $\kappa$. The tip pose $p(\kappa)$ traces a one-parameter curve as $\kappa$ varies, so the two tip coordinates cannot be commanded independently — only motion along the instantaneous Jacobian direction $J_\kappa$ is actuated. We therefore develop the controller in a scalar coordinate: either the curvature $\kappa$ or, equivalently for tracking against tissue, the tip displacement $y$ along the contact normal. The 2-D tip pose is recovered from the scalar state through the kinematics of Section 2.2. (Multi-segment catheters with $m$ independent tendons recover an $m$-DOF version; see Section 7.)

The dominant bending dynamics of a tendon-actuated catheter can be written in a form analogous to the Euler-Bernoulli beam with tip loading:

$$M(\kappa)\ddot{\kappa} + C(\kappa, \dot{\kappa})\dot{\kappa} + K(\kappa) = u + d_\text{cat}$$

where:
- $\kappa \in \mathbb{R}$ is the tip curvature (control output via forward kinematics)
- $M(\kappa) > 0$ is the effective bending inertia
- $C(\kappa, \dot{\kappa})$ captures velocity-dependent damping
- $K(\kappa)$ is the bending stiffness (nonlinear due to hysteresis)
- $u$ is the tendon-induced moment (control input)
- $d_\text{cat}$ lumps tissue contact forces, friction, backlash, and hysteresis residuals

This structure is isomorphic to the rigid-body manipulator equation

$$M(q)\ddot{q} + C(q,\dot{q})\dot{q} + G(q) = \tau + d$$

with the correspondence $\kappa \leftrightarrow q$, $K \leftrightarrow G$, identifying the theoretical bridge to the Impedance MPC framework developed for rigid-body systems.

**Remark (Modeling fidelity).** The full catheter mechanics are governed by the Cosserat rod equations, an infinite-dimensional PDE system. The lumped-parameter model above represents a deliberate model reduction: we control only the tip state $[\kappa, \dot\kappa]$ and absorb all higher-order distributed effects into the disturbance $d_\text{cat}$. This is a *Disturbance Compression* principle — not every physical phenomenon needs to be modeled explicitly; it only needs to be observable and slowly varying relative to the control bandwidth.

### 2.2 Tip Kinematics

For a constant-curvature segment of arc length $L$, the tip position in the base frame is:

$$p_x = \frac{\sin(\kappa L)}{\kappa}, \qquad p_z = \frac{1 - \cos(\kappa L)}{\kappa}$$

The translational Jacobian $J_\kappa \in \mathbb{R}^{2 \times 1}$ maps curvature rate to tip velocity:

$$\dot{p} = J_\kappa(\kappa)\,\dot{\kappa}, \qquad \ddot{p} = J_\kappa(\kappa)\,\ddot{\kappa} + \dot{J}_\kappa(\kappa, \dot{\kappa})\,\dot{\kappa}$$

The operational-space inertia scalar is:

$$\Lambda(\kappa) = \left(J_\kappa\,M(\kappa)^{-1}\,J_\kappa^\top\right)^{-1} \in \mathbb{R}_{>0}$$

This scalar normalizes the effective tip impedance to be configuration-adaptive: the same tendon force produces different tip accelerations depending on the current curvature, and $\Lambda(\kappa)$ accounts for this automatically.

### 2.3 Contact Model

During tissue contact, a tip contact force $F_\text{tip} \in \mathbb{R}^2$ acts at the catheter tip. In simulation we model the tissue as a Kelvin–Voigt spring-damper:

$$F_{\text{tissue},n}(t) = k_t\,\delta_n(t) + b_t\,\dot{\delta}_n(t)$$

where $\delta_n = \max(0, p_{n,\text{surface}} - p_n)$ is the normal penetration depth, $k_t$ is tissue stiffness (clinically $k_t \approx 2\text{–}20\,\text{kN/m}$ for cardiac tissue), and $b_t$ is tissue damping. These parameters are unknown in practice and are treated as part of the disturbance $d_\text{cat}$.

The contact force maps into curvature space via the Jacobian transpose:

$$\tau_\text{ext} = J_\kappa^\top\,F_\text{tip}$$

**Safety constraint.** The clinically relevant perforation threshold is approximately $F_\text{safe} = 0.5\,\text{N}$ [REF]. This constraint must be enforced as a hard bound:

$$\|F_\text{tip}(t)\| \leq F_\text{safe} \quad \forall t$$

### 2.4 Control Objective

Given a reference trajectory $(p_d(t), \dot{p}_d(t), \ddot{p}_d(t))$ planned offline, design control input $u(t)$ such that the tracking error $e(t) = p_d(t) - p(t)$ satisfies:

$$\lim_{t \to \infty} e(t) = 0$$

subject to:
- Tendon tension limits: $u_\text{lo} \leq u \leq u_\text{hi}$
- Contact force safety: $\|F_\text{tip}\| \leq F_\text{safe}$
- Curvature limits: $\kappa_\text{lo} \leq \kappa \leq \kappa_\text{max}$ (catheter physical limits)

in the presence of bounded, unknown contact forces, friction, and hysteresis.

---

## 3. Impedance MPC Controller Design

### 3.1 Two-Layer Architecture

The controller decomposes the tendon input as:

$$u = \underbrace{u_\text{ff}}_{\text{Layer 1: partial physics cancellation}} + \underbrace{J_\kappa^\top F_\text{mpc}}_{\text{Layer 2: MPC correction}}$$

**Layer 1 (feedforward)** cancels the *known* nominal dynamics at the current state:

$$u_\text{ff} = \hat{C}(\kappa, \dot{\kappa})\dot{\kappa} + \hat{K}(\kappa) + J_\kappa^\top \Lambda(\kappa)\,\ddot{p}_d$$

where $\hat{C}, \hat{K}$ denote nominal (possibly imperfect) model estimates.

**Layer 2 (receding-horizon QP)** acts on the linear residual plant, enforcing constraints and rejecting disturbances.

**Key distinction from rigid-body case.** For a well-modeled rigid robot, Layer 1 achieves *exact* cancellation, reducing the plant to a double integrator. For a catheter, $\hat{C}$ and $\hat{K}$ contain significant uncertainty. The partial cancellation reduces — but does not eliminate — the nonlinear terms. The residual modeling error is absorbed into the disturbance state $\hat{d}$, estimated by the Kalman filter. This is more robust than attempting exact cancellation with an unreliable model.

### 3.2 Error Dynamics and LPV Model

Substituting the two-layer input into the catheter equation of motion and forming $\ddot{e} = \ddot{p}_d - \ddot{p}$:

$$\ddot{e} = -\Lambda^{-1}(\kappa)\,F_\text{mpc} + \underbrace{J_\kappa M^{-1}\left[(\hat{C} - C)\dot\kappa + (\hat{K} - K)\right] + J_\kappa M^{-1}\tau_\text{ext} - \dot{J}_\kappa \dot\kappa}_{d(t)}$$

The disturbance $d(t)$ lumps:
- **Modeling error**: $(\hat{C} - C)\dot\kappa + (\hat{K} - K)$ — the mismatch between nominal and true catheter stiffness/damping
- **Contact force**: $J_\kappa M^{-1} \tau_\text{ext}$ — tissue, friction, hysteresis
- **Velocity coupling**: $\dot{J}_\kappa \dot\kappa$ — kinematic coupling term

Because the system is single-DOF (Section 2.1), the error $e = y_d - y$ along the controlled tip-normal coordinate is **scalar**, and $F_\text{mpc}\in\mathbb{R}$ is the scalar corrective force. Defining the error state $x_e = [e,\, \dot{e}]^\top \in \mathbb{R}^2$:

$$\dot{x}_e = \underbrace{\begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix}}_{A_c\,(\text{constant})} x_e + \underbrace{\begin{bmatrix} 0 \\ -\Lambda^{-1}(\kappa) \end{bmatrix}}_{B_c(\kappa)} F_\text{mpc} + \underbrace{\begin{bmatrix} 0 \\ 1 \end{bmatrix}}_{E_c\,(\text{constant})} d(t)$$

$A_c$ is **constant** — independent of catheter configuration; only the scalar input gain $B_c(\kappa)$ is parameter-varying through $\Lambda(\kappa)$. This is the same constant-$A_d$ LPV structure exploited in the rigid-body cases, enabling offline precomputation of the prediction matrices.

ZOH discretization at MPC sample time $\Delta t$ (exact, since $A_c$ is nilpotent):

$$x_e(k+1) = \underbrace{\begin{bmatrix} 1 & \Delta t \\ 0 & 1 \end{bmatrix}}_{A_d\,(\text{constant})} x_e(k) + \underbrace{\begin{bmatrix} -\tfrac{1}{2}\Lambda^{-1}(\kappa_k)\Delta t^2 \\ -\Lambda^{-1}(\kappa_k)\,\Delta t \end{bmatrix}}_{B_d(\kappa_k)} F_\text{mpc}(k)$$

### 3.3 Disturbance Augmentation

To drive the steady-state error to zero (in the nominal limit) under persistent contact forces, we augment the state with a scalar integrating disturbance $\hat{d} \in \mathbb{R}$:

$$\begin{bmatrix} x_e(k+1) \\ \hat{d}(k+1) \end{bmatrix} = \begin{bmatrix} A_d & G_d \\ 0 & 1 \end{bmatrix} \begin{bmatrix} x_e(k) \\ \hat{d}(k) \end{bmatrix} + \begin{bmatrix} B_d(\kappa_k) \\ 0 \end{bmatrix} F_\text{mpc}(k)$$

where the disturbance couples through $G_d = [\tfrac{1}{2}\Delta t^2,\, \Delta t]^\top$, the ZOH discretization of $E_c$ — **distinct from** the control matrix $B_d$, which carries the $-\Lambda^{-1}(\kappa_k)$ input gain. ($d$ is defined in acceleration units, so $G_d$ is *not* scaled by $\Lambda^{-1}$.) A steady-state Kalman filter estimates $\hat{d}(k)$ from the measured tip error $x_e(k)$, available from the EM tracker or vision system.

**Observability.** The augmented pair $(A_\text{aug}, C_\text{aug})$ with $C_\text{aug} = [I_2,\, 0_{2\times1}]$ (error-state measurement) is observable since $A_d$ has no repeated eigenvalues on the unit circle and $G_d \neq 0$. The disturbance state is therefore identifiable from tip-position measurements alone — no force sensor is required.

### 3.4 Receding-Horizon QP

At each MPC step, the controller solves:

$$\min_{U}\;\frac{1}{2}U^\top H\,U + x_e(k)^\top\,\Phi^\top \bar{Q}\,\Gamma\,U$$

where $U = [F_\text{mpc}(0); \ldots; F_\text{mpc}(N-1)] \in \mathbb{R}^{N}$ (scalar inputs), $H = \Gamma^\top \bar{Q}\,\Gamma + \bar{R}$, $Q = \text{diag}(K_d, D_d)$ encodes the desired impedance parameters, and $\Phi \in \mathbb{R}^{2N \times 2}$ is precomputed once at startup (since $A_d$ is constant). When no constraint is active the solution is the closed form $U^\star = -H^{-1}\Gamma^\top\bar Q(\Phi x_e + \Delta(\hat d))$, a matrix–vector multiply; when a constraint binds, OSQP solves the condensed QP reusing the cached factorization.

Subject to:
$$u_\text{lo} \leq u_\text{ff}(k) + J_\kappa^\top F_\text{mpc}(k) \leq u_\text{hi}, \quad \forall k \in [0, N-1]$$
$$\|F_\text{tip}(k)\|_\infty \leq F_\text{safe}, \quad \forall k \in [0, N-1]$$
$$\kappa_\text{lo} \leq \kappa(k) \leq \kappa_\text{max}, \quad \forall k \in [0, N-1]$$

The contact-force constraint $\|F_\text{tip}\| \leq F_\text{safe}$ is enforced through the transmitted tendon force. At contact equilibrium the feedforward cancels the nominal stiffness and damping, so the steady relation is $F_\text{tip} = u - \hat{F}_\text{fric}$, where $\hat{F}_\text{fric} = \Lambda\,\hat{d}_\text{free}$ is the (calibrated, free-space) friction/hysteresis bias separated from contact by the Kalman estimate. The hard bound therefore becomes a linear cap on the total command, $u_\text{ff}(k) + J_\kappa^\top F_\text{mpc}(k) \leq F_\text{safe} + \hat{F}_\text{fric}$, applied over the horizon. In the simulation of Section 6 this cap holds the contact force at $0.36\,\text{N}$ while every unconstrained controller exceeds the $0.5\,\text{N}$ bound.

### 3.5 Impedance MPC Equivalence

**Theorem 1 (Impedance MPC Equivalence for Catheters).** Let $Q = \text{blkdiag}(K_d, D_d)$ and $R \to 0$. In the absence of constraints and disturbances ($d \equiv 0$), the receding-horizon solution as $N \to \infty$ converges to the classical catheter impedance law:

$$u_\text{imp} = \hat{C}(\kappa, \dot\kappa)\dot\kappa + \hat{K}(\kappa) + J_\kappa^\top\left(K_d\,e + D_d\,\dot{e}\right)$$

with effective tip impedance:

$$Z_\text{eff}(s, \kappa) \approx \Lambda(\kappa)\left(s^2 + D_\text{eff}(\kappa)\,s + K_\text{eff}(\kappa)\right)$$

where $\Lambda(\kappa)$ provides automatic configuration-adaptive mass normalization.

*Proof sketch.* Without constraints and with $R \to 0$, the QP collapses to unconstrained minimization of $\sum_k x_e(k)^\top Q\,x_e(k)$. The optimal infinite-horizon solution is the LQR gain $K_\text{LQR} = [K_d, D_d]$ for the double-integrator $B_d$ plant, recovering $F_\text{mpc} = K_d e + D_d \dot{e}$ at the nominal configuration $\kappa_0$.  $\square$

When constraints are active or $d \neq 0$, the MPC departs from the unconstrained impedance law — this departure is precisely its advantage over classical impedance control.

---

## 4. Stability Analysis

### 4.1 Nominal Stability

**Theorem 2 (Closed-Loop Stability).** Consider the system controlled by the receding-horizon law with terminal cost $Q_f$ chosen as the DARE solution at a nominal configuration $\kappa_0$. If the QP is feasible at $k=0$ and the LPV variation $\|B_d(\kappa_k) - B_d(\kappa_0)\|$ is sufficiently small, then the closed-loop system is asymptotically stable and $x_e(k) \to 0$ as $k \to \infty$ in the absence of disturbances.

### 4.2 Offset-Free Tracking Under Contact Forces

**Theorem 3 (ISS and Offset-Free Tracking).** If the lumped disturbance satisfies $\|d(t)\| \leq \bar{d}$ for all $t$, and the augmented system is observable, then the augmented state $(x_e, \hat{d})$ is input-to-state stable and the steady-state position error satisfies $\|e_\infty\| \to 0$ as the Kalman filter converges.

*The key condition* — $d(t)$ bounded and slowly varying — requires discussion in the catheter context (Section 5.2).

---

## 5. Discussion

### 5.1 The Role of MPC Prediction in the Catheter Setting

A natural concern is whether MPC look-ahead provides benefit when contact force predictions over the horizon are unreliable. We distinguish three sources of prediction information:

| Prediction source | Quality | Benefit |
|---|---|---|
| Force prediction $\hat{F}(k+i)$ | Poor (tissue unknown) | Limited |
| Trajectory feedforward $\ddot{p}_d(k+i)$ | Exact (offline plan) | **High** |
| Constraint anticipation | Exact (known bounds) | **High** |
| Kalman disturbance estimate at $k$ | Moderate | **High** |

The primary value of MPC in this setting is **constraint enforcement** (preventing force limit violation before it occurs) and **trajectory feedforward** (anticipating reference curvature), not contact force prediction. A zero-order hold $\hat{d}(k+i) = \hat{d}(k)$ is adopted over the horizon; the offset-free guarantee is preserved at steady state regardless of prediction quality.

This is a fundamental difference from the excavator setting, where the Kelvin-Voigt soil model provides a physically meaningful prediction structure. In the catheter setting, the force prediction is essentially constant over the horizon, and the MPC's advantage over a reactive impedance controller comes primarily from the other two channels.

### 5.2 The $\dot{d} = 0$ Assumption and Cardiac Motion

The integrating disturbance model $\hat{d}(k+1) = \hat{d}(k)$ is a worst-case approximation that guarantees offset-free tracking for *any* slowly varying disturbance. For cardiac catheter procedures, the contact force has a quasi-periodic component at ~1 Hz (heart rate). Two refinements are possible:

**Option A (engineering pragmatics): Kalman tuning.** Increase the process noise covariance $Q_w$ of the disturbance state to make the Kalman filter more responsive to rapid changes. This trades estimation smoothness for faster adaptation to cardiac-induced force transients, at the cost of increased noise sensitivity.

**Option B (periodic disturbance model).** If the cardiac cycle phase $\phi(t)$ is available (e.g., from ECG gating), the disturbance model can be extended to:

$$\hat{d}(k+1) = \hat{d}_\text{DC}(k) + A_\text{card}\sin(\omega_\text{heart}(k+1)\Delta t + \phi)$$

This turns the cardiac force into a predictable component, improving MPC look-ahead quality substantially. This is an extension left for future work.

**Contact force transients.** When the catheter makes first contact with tissue, the force rises rapidly — this is the scenario where $\dot{d} = 0$ is most violated. The Kalman filter will lag behind by several sample periods. The hard constraint $\|F_\text{tip}\| \leq F_\text{safe}$ in the QP provides a safety backstop during this transient: even if the disturbance estimate has not yet converged, the MPC will reduce $F_\text{mpc}$ to respect the force bound.

### 5.3 Partial vs. Full Physics Cancellation

Unlike the excavator case where $M(q), C(q,\dot q), G(q)$ can be computed with reasonable accuracy from CAD models, catheter stiffness $K(\kappa)$ is subject to hysteresis, temperature dependence, and specimen variability. Attempting full inversion with an inaccurate $\hat{K}$ can *increase* the effective disturbance rather than reduce it.

The partial cancellation approach adopted here — canceling only the components that can be estimated reliably — is more robust. The residual modeling error simply adds to the disturbance $d(t)$ that the Kalman filter estimates. As long as the total disturbance remains bounded and observable, the offset-free guarantee holds regardless of how large the modeling error is.

A practical design rule: include a term in $u_\text{ff}$ only if its omission would cause the steady-state disturbance to exceed the Kalman filter's tracking bandwidth.

### 5.4 Sensing Requirements

The framework requires tip position measurements $p(t)$ for Kalman filter updates. These are available from:

- **EM tracking**: 10–40 Hz, ~1 mm accuracy — sufficient for Kalman correction
- **Fluoroscopy**: intermittent, radiation dose constraints limit continuous use
- **ICE / ultrasound**: operator-dependent positioning
- **FBG shape sensing**: continuous, ~0.5 mm accuracy — preferred if available

A force sensor is *not required* by the framework. If FBG-based force estimation is available, the estimated $\hat{F}_\text{tip}$ can be injected directly into the augmented state as in the excavator pressure-sensor path, reducing Kalman estimation lag during contact transients.

---

## 6. Simulation Results

### 6.1 Simulation Setup

We implement the single-DOF tip-normal model of Section 2 in Python (NumPy/SciPy); the residual plant after partial-physics feedforward is the scalar double integrator $\Lambda\,\ddot y = u + F_\text{ext} - b_\text{eff}\dot y - k_\text{eff} y$, integrated with RK4. Parameters: effective operational-space tip inertia $\Lambda = 1.0$, nominal damping $b_\text{eff} = 5$ and stiffness $k_\text{eff} = 50$ (both cancelled by the feedforward), and a constant tendon-hysteresis bias $F_\text{fric} = 0.45\,\text{N}$ representing the unmodelled friction/hysteresis the disturbance state must absorb. Tissue is a Kelvin–Voigt wall at $6\,\text{mm}$ with $k_t = 5\,\text{kN/m}$, $b_t = 40\,\text{N\,s/m}$. The MPC runs at 500 Hz ($\Delta t = 2\,\text{ms}$) with horizon $N = 20$; tendon force is limited to $\pm 8\,\text{N}$ and the safety bound is $F_\text{safe} = 0.5\,\text{N}$.

The reference trajectory approaches the tissue surface over 1 s, presses to a firm-contact target $0.9\,\text{mm}$ *past* the surface (a depth that would require $\approx$4.5 N if tracked rigidly), holds for 2 s (ablation dwell), and retracts. The penetrating target deliberately stresses the force-safety mechanism: reaching it exactly is unsafe, so a correct controller must trade position for force.

The benchmark is reproducible: `simulation/catheter_benchmark.py`.

### 6.2 Controller Comparison

Four controllers are compared, all sharing the partial-physics feedforward where applicable: **classical impedance** ($K_d = 900\,\text{N/m}$, $D_d = 60\,\text{N\,s/m}$); **Impedance MPC (no force constraint)** with Kalman augmentation; **Impedance MPC (with force constraint)** enforcing $F_\text{safe} = 0.5\,\text{N}$ via a tendon-force cap derived from the Kalman friction/contact estimate; and **joint-space PD** (no feedforward).

| Controller | Approach RMS (mm) | Max contact force (N) | $F_\text{safe}$ violated? | Hold pos. error (mm) |
|---|:---:|:---:|:---:|:---:|
| Classical impedance | 0.50 | 0.79 | **Yes** | 1.35 |
| Impedance MPC (no FC) | 0.30 | 2.04 | **Yes** | 1.30 |
| Impedance MPC (with FC) | 0.30 | **0.36** | No | 1.46 |
| Joint-space PD | 0.69 | 0.55 | **Yes** | 1.40 |

*Approach RMS* is the free-space tracking error before contact; *max contact force* and *hold position error* are during the contact dwell. Key observations:

- **Offset-free rejection of the friction bias (free space).** Classical impedance leaves a free-space error of $0.50\,\text{mm}$, matching the analytic prediction $e_\infty = F_\text{fric}/K_d = 0.45/900 = 0.50\,\text{mm}$. The Kalman augmentation rejects this constant bias, cutting the approach error to $0.30\,\text{mm}$ — a 40% reduction. This is where offset-free tracking helps.

- **The force-safety constraint is the decisive contribution.** Only the force-constrained MPC respects the $0.5\,\text{N}$ bound (peak $0.36\,\text{N}$). Classical impedance reaches $0.79\,\text{N}$ (58% over), and joint-space PD $0.55\,\text{N}$.

- **Offset-free tracking *without* a force limit is dangerous in contact.** The unconstrained Impedance MPC, driving toward the penetrating target, pushes to $2.04\,\text{N}$ — over 4× the safety bound — precisely because it is trying to eliminate the position error against stiff tissue. This is the central cautionary result: offset-free disturbance rejection and contact-force safety are in direct tension, and only the explicit hard constraint resolves it.

- **Contact-phase position error is contact-limited, not controller-limited.** All controllers retain $\approx$1.3–1.5 mm hold error because the commanded depth cannot be reached without exceeding $F_\text{safe}$; the force-constrained MPC accepts marginally more error ($1.46\,\text{mm}$) as the explicit, intentional price of safety. The simulation thus does **not** support a blanket "zero steady-state error" claim under stiff contact — offset-free holds in free space and in the nominal infinite-horizon limit, but stiff-contact tracking is fundamentally force-limited.

Figure (`simulation/catheter_results.png`) shows the tip trajectories and contact-force histories: the unconstrained MPC's force transient peaks at 2 N while the constrained MPC alone remains below the 0.5 N line throughout.

### 6.3 Disturbance Estimation Quality

The Kalman filter converges to the friction/contact disturbance estimate within a few sample periods. During the contact-onset transient the estimate lags, but the tendon-force cap prevents force-bound violation even before convergence — the constraint, not the estimate, provides the safety guarantee. This mirrors the role separation found in the rigid-body cases: the disturbance estimate improves *tracking*, while the hard constraint provides *safety*.

---

## 7. Conclusion

We have presented an Impedance MPC framework for steerable catheter tip control that extends classical impedance control with constraint enforcement, offset-free disturbance rejection, and configuration-adaptive tip compliance. The key theoretical result establishes that Impedance MPC recovers the classical impedance law in the unconstrained case while adding capabilities that are clinically essential: hard contact force safety constraints and zero steady-state position error under persistent tissue contact.

The framework inherits the two-layer feedforward-correction architecture from the rigid-body case, adapted to catheter mechanics through a partial physics cancellation that is robust to the large modeling uncertainty inherent in catheter systems. The disturbance compression principle — lumping tissue contact forces, friction, hysteresis, and modeling error into a single Kalman-estimated state — eliminates the need for explicit modeling of these phenomena.

We have honestly characterized the limitations of MPC in this setting: contact force prediction quality over the horizon is fundamentally limited by unknown tissue geometry and cardiac motion. However, the primary value of MPC in the catheter application lies in constraint enforcement and trajectory feedforward rather than force prediction — and these benefits are preserved regardless of prediction quality.

**Future work** includes: (i) hardware validation on a physical tendon-actuated catheter with EM tracking; (ii) integration of FBG-based force estimation as a direct disturbance channel; (iii) cardiac-phase-aware periodic disturbance models for improved prediction during ablation; (iv) extension to multi-segment catheters with configuration-varying $\Lambda(\kappa)$; and (v) energy-tank augmentation for certified passivity during aggressive contact maneuvers.

---

## References

- [HOGAN1985] N. Hogan, "Impedance control: An approach to manipulation," ASME J. Dyn. Syst. Meas. Control, 1985.
- [ANTMAN1995] S. S. Antman, *Nonlinear Problems of Elasticity*, Springer, 1995.
- [HAN2009] J. Han, "From PID to active disturbance rejection control," IEEE Trans. Ind. Electron., 2009.
- [PANNOCCHIA2003] G. Pannocchia and J. B. Rawlings, "Disturbance models for offset-free model-predictive control," AIChE J., 2003.
- [MAEDER2010] U. Maeder, F. Borrelli, and M. Morari, "Linear offset-free model predictive control," Automatica, 2010.
- [RAWLINGS2017] J. B. Rawlings, D. Q. Mayne, and M. Diehl, *Model Predictive Control: Theory, Computation, and Design*, Nob Hill, 2017.
- [CAO2005] Y. Cao and A. Bhatt, "Min-max MPC for LPV systems," *Proc. ACC*, 2005.

---

*— Draft for internal discussion. Simulation results in Section 6 are produced by `simulation/catheter_benchmark.py` (single-DOF model, NumPy/SciPy). —*
