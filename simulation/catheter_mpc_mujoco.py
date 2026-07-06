#!/usr/bin/env python3
"""
Impedance MPC ported to the MuJoCo PRB catheter — step-2 of the MuJoCo port
==========================================================================
Closes the loop: the Impedance-MPC of `catheter_benchmark.py` (scalar RK4) is run
against the distributed-compliance MuJoCo PRB catheter of `catheter_mujoco.py`,
with the two changes the spike found necessary:

  (1) Lambda-rescaled gains. The MuJoCo tip inertia is Lambda ~ 3.5e-3, ~280x lighter
      than the scalar benchmark's normalized 1.0, so keeping Kd=900 would place the
      closed loop at omega ~ 500 rad/s (unstable at the 500 Hz control rate). We rescale
      Kd -> Kd*Lambda, Dd -> Dd*Lambda to preserve the validated omega ~ 30 rad/s.

  (2) J_k-remapped actuator mapping. The MPC works in tip-normal FORCE units (F_mpc is a
      real tip force, so the hard bound |F_mpc| <= F_safe is dimensionally correct); the
      tendon command is T = -F_mpc / J_k(kappa), so the J_k remap lives in the actuator
      map, not the constraint. Elasticity is absorbed by the offset-free integral state.

We run an approach->press->hold->retract trajectory against a wall and compare the
force-constrained vs unconstrained controller, reproducing the Table II safety story on
the physics plant. Uses the d_tend=2.5 mm transmission variant where the bound binds.

Run: python3 catheter_mpc_mujoco.py
"""
import numpy as np
import mujoco
import catheter_mujoco as C       # plant: build_xml, tip_state, tissue_force, apply_tissue, etc.

F_SAFE = C.F_SAFE                 # 0.5 N
DT_CTRL = 0.002                   # 500 Hz control
Y_WALL = 0.004                    # wall surface depth (tip-normal, downward +)

# ── desired tip-normal trajectory: approach -> press -> hold -> retract ────────
def y_ref(t):
    y_press = Y_WALL + 1.5e-3                      # 1.5 mm PAST surface (stresses safety)
    if t < 1.0:                                    # approach to surface
        a = 0.5*(1-np.cos(np.pi*t)); return a*Y_WALL, 0.5*np.pi*np.sin(np.pi*t)*Y_WALL
    if t < 1.5:                                    # press into penetrating target
        s = (t-1.0)/0.5; a = 0.5*(1-np.cos(np.pi*s))
        return Y_WALL + a*1.5e-3, 0.5*np.pi*np.sin(np.pi*s)/0.5*1.5e-3
    if t < 2.5:                                    # hold
        return y_press, 0.0
    if t < 3.5:                                    # retract
        s = (t-2.5)/1.0; a = 0.5*(1-np.cos(np.pi*s))
        return y_press*(1-a), -0.5*np.pi*np.sin(np.pi*s)*y_press
    return 0.0, 0.0


def measure_plant(model, data):
    """Lambda (tip inertia), local J_k = dF_tip/dT, and tip-normal elastic stiffness
    k_eff = F_tip/deflection (for the Layer-1 feedforward that cancels the catheter's
    own restoring force -- which the scalar benchmark lumped into k_eff and cancelled)."""
    mujoco.mj_resetData(model, data); data.joint("wall_z").qpos[0] = -0.06
    C.settle(model, data, -2.0, -0.06)
    lam = C.operational_inertia(model, data)
    def press(T):                                  # in-contact -> J_k = dF_tip/dT
        mujoco.mj_resetData(model, data); return C.settle(model, data, -abs(T), 0.0)
    J_k = (press(6.0) - press(4.0)) / 2.0
    def deflect(T):                                # free-space deflection at tension T
        mujoco.mj_resetData(model, data); data.joint("wall_z").qpos[0] = -0.06
        C.settle(model, data, -abs(T), -0.06); return C.tip_state(model, data)[0]
    # k_eff_tip from a small free-space deflection near the 4-5 mm operating range
    y_lo = deflect(0.5)
    k_eff = (J_k*0.5) / y_lo if y_lo > 1e-6 else 0.0
    return lam, J_k, k_eff


def run_mpc(model, data, lam, J_k, k_eff, constrain, Ki=300.0, T_end=3.5,
           cardiac_amp=0.0, f_heart=1.0, pos_noise=0.0, seed=0, ref=y_ref):
    """Closed-loop Impedance MPC on the MuJoCo plant. Returns time, e(mm), F_contact(N), y(mm).
    cardiac_amp/f_heart drive the wall (z-slide); pos_noise adds tip-measurement noise;
    ref(t)->(y,yd) is the tip-normal reference (default approach-press-hold)."""
    Kd, Dd = 900.0*lam, 60.0*lam      # Lambda-rescaled impedance (preserves omega~30 rad/s)
    # Ki = offset-free integral; winds up against the (unreachable) penetrating target,
    # which is what drives the contact force toward the limit and stresses safety.
    n_sub = int(round(DT_CTRL/model.opt.timestep))
    nctrl = int(T_end/DT_CTRL)
    mujoco.mj_resetData(model, data)
    rng = np.random.default_rng(seed)
    dhat = 0.0
    ts, es, Fs, ys = [], [], [], []
    for k in range(nctrl):
        t = k*DT_CTRL
        y, yd = C.tip_state(model, data)
        y_meas = y + (rng.normal(0, pos_noise) if pos_noise > 0 else 0.0)
        yr, yrd = ref(t)
        e, de = yr - y_meas, yrd - yd
        dhat += Ki*e*DT_CTRL                          # integral disturbance estimate
        F_mpc = Kd*e + Dd*de + dhat                   # impedance + offset-free, tip-force units
        if constrain:                                 # hard bound on the CONTACT force
            F_mpc = float(np.clip(F_mpc, -F_SAFE, F_SAFE))
        ff = k_eff*yr                                 # Layer 1: cancel catheter elasticity
        T = float(np.clip(-(ff+F_mpc)/J_k, -C.TENDON_MAX, 0.0))  # J_k remap -> tendon tension
        data.ctrl[0] = T
        Fc = 0.0
        for _ in range(n_sub):                        # step physics with analytic tissue
            data.ctrl[1] = cardiac_amp*np.sin(2*np.pi*f_heart*data.time)  # wall qpos (cardiac)
            f, _ = C.tissue_force(model, data); C.apply_tissue(model, data, f)
            mujoco.mj_step(model, data)
            Fc = max(Fc, f)
        ts.append(t); es.append(e*1e3); Fs.append(Fc); ys.append(C.tip_state(model, data)[0]*1e3)
    return map(np.array, (ts, es, Fs, ys))


def ref_sin(t):
    """Free-space sinusoidal tip-normal reference (0.5-3.5 mm, below the 4 mm wall) ->
    many reversals, where Coulomb tendon friction bites hardest."""
    w = 2*np.pi*0.5                       # 0.5 Hz
    y0, A = 2.0e-3, 1.5e-3
    return y0 + A*np.sin(w*t), A*w*np.cos(w*t)


def friction_test(lam, J_k, k_eff):
    """Does offset-free d^ absorb a HYSTERETIC (Coulomb) tendon friction, or only the
    scalar model's constant bias? Controller is calibrated on the nominal (frictionless)
    plant; the real plant adds tendon frictionloss -> the disturbance d^ must absorb it."""
    print("\n" + "="*68)
    print("Tendon-friction (hysteresis) test — offset-free d^ vs Coulomb friction")
    print("   free-space sinusoid 0.5-3.5 mm @ 0.5 Hz; controller uses nominal model")
    print(f"   {'plant / controller':<42}{'track RMS (mm)':>16}")
    for fr, fr_lbl in [(0.0, "no friction"), (0.3, "tendon friction 0.3 N")]:
        model = mujoco.MjModel.from_xml_string(
            C.build_xml(d_tend=0.0025, tendon_max=8.0, tendon_friction=fr))
        data = mujoco.MjData(model)
        for Ki, ki_lbl in [(0.0, "impedance only"), (300.0, "+ offset-free")]:
            t, e, _, _ = run_mpc(model, data, lam, J_k, k_eff, constrain=False,
                                 Ki=Ki, T_end=4.0, ref=ref_sin)
            ss = t >= 1.0                                 # after first transient
            rms = float(np.sqrt(np.mean(e[ss]**2)))
            print(f"   {fr_lbl+', '+ki_lbl:<42}{rms:>16.3f}")
    print("   -> constant bias is trivially absorbed; the question is the SIGN-DEPENDENT")
    print("      part. Offset-free cuts the error but cannot fully cancel it (d^ lags each")
    print("      reversal) -> residual hysteresis error = the real cost the scalar model's")
    print("      constant-F_fric hid. Motivates a friction/hysteresis-aware observer.")
    print("="*68)


def main():
    model = mujoco.MjModel.from_xml_string(C.build_xml(d_tend=0.0025, tendon_max=8.0))
    data = mujoco.MjData(model)
    lam, J_k, k_eff = measure_plant(model, data)

    print("="*68)
    print("Impedance MPC ported to MuJoCo PRB catheter (d_tend=2.5mm variant)")
    print(f"  measured Lambda={lam:.3e}  J_k={J_k:.3f}  k_eff={k_eff:.1f} N/m  "
          f"F_safe={F_SAFE} N  ctrl=500 Hz")
    print(f"  Lambda-rescaled gains: Kd={900*lam:.2f}  Dd={60*lam:.3f} (from 900/60); "
          f"FF cancels k_eff")
    print("="*68)

    res = {}
    for name, con in [("MPC + force constraint", True), ("MPC, no constraint", False)]:
        t, e, F, y = run_mpc(model, data, lam, J_k, k_eff, con)
        appr = (t >= 0.3) & (t < 0.95)
        hold = (t >= 1.8) & (t < 2.5)
        rms_appr = float(np.sqrt(np.mean(e[appr]**2)))
        maxF = float(np.max(F))
        hold_e = float(np.mean(np.abs(e[hold])))
        viol = maxF > F_SAFE + 1e-3
        res[name] = (t, e, F, y, rms_appr, maxF, hold_e, viol)
        print(f"\n  {name}")
        print(f"    approach RMS = {rms_appr:.3f} mm | peak contact F = {maxF:.3f} N "
              f"({'VIOLATES' if viol else 'safe'}) | hold err = {hold_e:.3f} mm")

    print("\n" + "-"*68)
    cF = res["MPC + force constraint"][5]; uF = res["MPC, no constraint"][5]
    print(f"  Force-safety: constrained peak {cF:.3f} N vs unconstrained {uF:.3f} N "
          f"(F_safe={F_SAFE}).")
    print("  Story matches Table II: offset-free drive pushes the contact force up to hit")
    print("  the penetrating target; only the hard |F_mpc|<=F_safe bound (via T=-F_mpc/J_k)")
    print("  keeps it safe. Tracking is on the physics plant with Lambda-rescaled gains.")

    # ── Cardiac re-validation (Table V) on the physics plant, force-constrained mode ──
    print("\n" + "="*68)
    print("Cardiac safety-mode re-validation (Table V) — force-constrained MPC")
    print(f"   {'wall condition':<40}{'apprRMS':>8}{'peakF':>8}{'viol':>6}")
    for label, amp, fh, noise in [
            ("static wall (baseline)",                 0.0,    1.0, 0.0),
            ("0.3 mm @ 1 Hz",                          0.3e-3, 1.0, 0.0),
            ("0.5 mm @ 1.2 Hz + 0.2 mm noise",         0.5e-3, 1.2, 0.2e-3)]:
        t, e, F, y = run_mpc(model, data, lam, J_k, k_eff, constrain=True,
                             cardiac_amp=amp, f_heart=fh, pos_noise=noise)
        appr = (t >= 0.3) & (t < 0.95)
        rms = float(np.sqrt(np.mean(e[appr]**2))); maxF = float(np.max(F))
        print(f"   {label:<40}{rms:>8.3f}{maxF:>8.3f}"
              f"{('YES' if maxF > F_SAFE+1e-3 else 'no'):>6}")
    print("   -> NOTE: the bound holds, but mainly because the catheter (k_eff~8 N/m) is")
    print("      ~600x softer than the tissue (5 kN/m): the compliant tip RIDES with the")
    print("      moving wall, so a 0.5 mm excursion barely changes the contact force")
    print("      (peak ~unchanged from static). The cardiac motion is BUFFERED by catheter")
    print("      compliance, not stress-testing the constraint -- the opposite regime from")
    print("      the stiffer RK4 scalar plant. A refinement of, not a match to, RK4 Table V.")
    print("="*68)

    friction_test(lam, J_k, k_eff)


if __name__ == "__main__":
    main()
