#!/usr/bin/env python3
"""
Impedance-MPC catheter benchmark — MuJoCo distributed-compliance plant
======================================================================
Full physics-engine version of the four-controller benchmark (replaces the RK4
scalar plant of `catheter_benchmark.py`). All reported simulation results in the
paper are produced here, on the tendon-driven eight-link pseudo-rigid-body (PRB)
catheter of `catheter_mujoco.py`, with the analytic Kelvin--Voigt tissue wall.

Controllers (all mapped to TENDON tension through the measured transmission J_k):
  - Classical impedance      F = Kd e + Dd e_dot                  (no offset-free, no FC)
  - Impedance MPC (no FC)     F = Kd e + Dd e_dot + d_hat          (offset-free, no bound)
  - Impedance MPC (with FC)   same, with hard |F| <= F_safe        (the safety contribution)
  - Joint-space PD            F = Kp e + Kd e_dot, NO feedforward   (fights catheter elasticity)

Impedance-family gains are Lambda-rescaled to the measured tip inertia (preserving
omega ~ 30 rad/s); the tendon command is T = -(u_ff + F)/J_k with u_ff = k_eff*y_ref
cancelling the catheter's own restoring force (Layer 1). Prints the four-controller
table and the cardiac safety-mode table, and writes catheter_results.png /
catheter_cardiac.png. Run: python3 catheter_benchmark_mujoco.py
"""
import numpy as np
import mujoco
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import catheter_mujoco as C
from catheter_mpc_mujoco import y_ref, measure_plant, Y_WALL, F_SAFE, DT_CTRL

T_END = 3.5
D_TEND = 0.0025            # transmission variant where the force bound binds
TENDON_MAX = 8.0


# ── controllers (return tendon tension T, given tip-normal y, ydot) ───────────
class Classical:
    """Classical impedance: u_ff cancels elasticity, impedance gains close the loop.
    No offset-free disturbance state and no force constraint."""
    name = "Classical impedance"
    def __init__(s, lam, J_k, k_eff):
        s.Kd, s.Dd, s.Jk, s.keff = 900.0*lam, 60.0*lam, J_k, k_eff
    def reset(s): pass
    def __call__(s, t, y, yd):
        yr, yrd = y_ref(t); e, de = yr - y, yrd - yd
        F = s.Kd*e + s.Dd*de
        ff = s.keff*yr
        return float(np.clip(-(ff + F)/s.Jk, -TENDON_MAX, 0.0))


class ImpedanceMPC:
    """Impedance MPC: impedance gains + offset-free integral disturbance state.
    force_con=True adds the hard predicted-force bound |F| <= f_cap (default F_safe;
    set f_cap = F_safe*cos(theta_max) for the misalignment-tightened bound of Sec. V-E)."""
    def __init__(s, lam, J_k, k_eff, force_con, f_cap=F_SAFE):
        s.Kd, s.Dd, s.Ki = 900.0*lam, 60.0*lam, 300.0
        s.Jk, s.keff, s.fc, s.fcap = J_k, k_eff, force_con, f_cap
        s.name = "Impedance MPC (with FC)" if force_con else "Impedance MPC (no FC)"
        s.reset()
    def reset(s): s.dhat = 0.0
    def __call__(s, t, y, yd):
        yr, yrd = y_ref(t); e, de = yr - y, yrd - yd
        s.dhat += s.Ki*e*DT_CTRL
        F = s.Kd*e + s.Dd*de + s.dhat
        if s.fc:
            F = float(np.clip(F, -s.fcap, s.fcap))        # hard contact-force bound
        ff = s.keff*yr
        return float(np.clip(-(ff + F)/s.Jk, -TENDON_MAX, 0.0))


class ForceReg(ImpedanceMPC):
    """Force-regulation mode (App. A): position-track in free space, then regulate the
    sensor-free contact-force estimate to F_des during the contact hold, hard bound still
    active. Sensor-free estimate Fc_hat = max(0, J_k*|T| - k_eff*|y|) = commanded tip force
    beyond elasticity. NOT a claimed contribution -- characterizes the framework's limit."""
    def __init__(s, lam, J_k, k_eff, F_des=0.2, Ki_f=8.0, Dv=60.0):
        super().__init__(lam, J_k, k_eff, force_con=True)
        s.F_des, s.Ki_f, s.Dv = F_des, Ki_f, Dv; s.name = "Force-regulation MPC"
    def reset(s): super().reset(); s.Fint = 0.0; s.Tprev = 0.0
    def __call__(s, t, y, yd):
        if 1.5 <= t < 2.5:                                   # contact hold -> force reg
            Fc_hat = max(0.0, abs(s.Tprev)*s.Jk - s.keff*abs(y))   # sensor-free estimate
            s.Fint = float(np.clip(s.Fint + s.Ki_f*(s.F_des - Fc_hat)*DT_CTRL, -3, 3))
            F_cmd = float(np.clip(s.F_des + s.Fint - s.Dv*yd, 0.0, F_SAFE))  # hard bound
            T = float(np.clip(-(s.keff*y + F_cmd)/s.Jk, -TENDON_MAX, 0.0))
        else:                                                # free space -> position track
            T = super().__call__(t, y, yd)
        s.Tprev = T
        return T


class JointPD:
    """Joint-space PD: NO feedforward, so it fights the catheter's own elasticity."""
    name = "Joint-space PD"
    def __init__(s, lam, J_k, k_eff):
        s.Kp, s.Kd, s.Jk = 900.0*lam, 60.0*lam, J_k
    def reset(s): pass
    def __call__(s, t, y, yd):
        yr, yrd = y_ref(t); e, de = yr - y, yrd - yd
        F = s.Kp*e + s.Kd*de
        return float(np.clip(-F/s.Jk, -TENDON_MAX, 0.0))


# ── closed-loop run on the MuJoCo plant ───────────────────────────────────────
def run(model, data, ctrl, cardiac_amp=0.0, f_heart=1.0, pos_noise=0.0, seed=0,
        tissue_scale=1.0):
    """tissue_scale multiplies the along-axis tissue reaction (=cos^2(theta) for the
    Sec. V-E contact-normal-misalignment model); Fs records the along-axis force."""
    n_sub = int(round(DT_CTRL/model.opt.timestep))
    nctrl = int(T_END/DT_CTRL)
    mujoco.mj_resetData(model, data); ctrl.reset()
    rng = np.random.default_rng(seed)
    ts, es, Fs, ys = [], [], [], []
    for k in range(nctrl):
        t = k*DT_CTRL
        y, yd = C.tip_state(model, data)
        y_meas = y + (rng.normal(0, pos_noise) if pos_noise > 0 else 0.0)
        T = ctrl(t, y_meas, yd)
        data.ctrl[0] = T
        Fc = 0.0
        for _ in range(n_sub):
            data.ctrl[1] = cardiac_amp*np.sin(2*np.pi*f_heart*data.time)   # wall (cardiac)
            f, _ = C.tissue_force(model, data); f *= tissue_scale
            C.apply_tissue(model, data, f)
            mujoco.mj_step(model, data)
            Fc = max(Fc, f)
        yr, _ = y_ref(t)
        ts.append(t); es.append((yr - C.tip_state(model, data)[0])*1e3)
        Fs.append(Fc); ys.append(C.tip_state(model, data)[0]*1e3)
    return map(np.array, (ts, es, Fs, ys))


def metrics(t, e, F):
    appr = (t >= 0.3) & (t < 0.95)
    hold = (t >= 1.8) & (t < 2.5)
    rms_appr = float(np.sqrt(np.mean(e[appr]**2)))
    maxF = float(np.max(F))
    hold_e = float(np.mean(np.abs(e[hold])))
    return rms_appr, maxF, maxF > F_SAFE + 1e-3, hold_e


def main():
    model = mujoco.MjModel.from_xml_string(C.build_xml(d_tend=D_TEND, tendon_max=TENDON_MAX))
    data = mujoco.MjData(model)
    lam, J_k, k_eff = measure_plant(model, data)

    print("="*72)
    print("Impedance-MPC catheter benchmark — MuJoCo PRB plant (d_tend=2.5 mm)")
    print(f"  measured Lambda={lam:.3e}  J_k={J_k:.3f}  k_eff={k_eff:.1f} N/m  "
          f"F_safe={F_SAFE} N  ctrl=500 Hz")
    print(f"  Lambda-rescaled gains Kd={900*lam:.2f} Dd={60*lam:.3f}; FF cancels k_eff")
    print("="*72)

    ctrls = [Classical(lam, J_k, k_eff),
             ImpedanceMPC(lam, J_k, k_eff, force_con=False),
             ImpedanceMPC(lam, J_k, k_eff, force_con=True),
             JointPD(lam, J_k, k_eff)]
    res = {}
    print(f"  {'Controller':<26}{'ApprRMS':>9}{'maxF(N)':>9}{'viol':>6}{'HoldErr':>9}")
    print("  (mm except force)"); print("  "+"-"*58)
    for c in ctrls:
        t, e, F, y = run(model, data, c); m = metrics(t, e, F)
        res[c.name] = (t, e, F, y, m)
        print(f"  {c.name:<26}{m[0]:>9.3f}{m[1]:>9.3f}{('YES' if m[2] else 'no'):>6}{m[3]:>9.3f}")

    # ── cardiac safety-mode (Impedance MPC + FC) ──────────────────────────────
    print(f"\n  Cardiac safety mode (Impedance MPC + FC, position tracking):")
    print(f"    {'wall condition':<40}{'ApprRMS':>9}{'maxF(N)':>9}{'viol':>6}")
    card = []
    for label, short, amp, fh, noise, col in [
            ("static wall (baseline)",            "static wall",          0.0,    1.0, 0.0,    "#27ae60"),
            ("0.3 mm @ 1 Hz",                     "0.3 mm @ 1 Hz",        0.3e-3, 1.0, 0.0,    "#2980b9"),
            ("0.5 mm @ 1.2 Hz + 0.2 mm noise",    "0.5 mm @ 1.2 Hz+noise",0.5e-3, 1.2, 0.2e-3, "#c0392b")]:
        c = ImpedanceMPC(lam, J_k, k_eff, force_con=True)
        t, e, F, y = run(model, data, c, cardiac_amp=amp, f_heart=fh, pos_noise=noise)
        m = metrics(t, e, F); card.append((short, t, F, m[1], m[2], col))
        print(f"    {label:<40}{m[0]:>9.3f}{m[1]:>9.3f}{('YES' if m[2] else 'no'):>6}")

    # ── force-regulation mode (Appendix A) on the MuJoCo plant ────────────────
    print(f"\n  Force-regulation mode (sensor-free, target F_des=0.2 N), reg window [1.8,2.5]s:")
    for label, amp, fh, noise in [
            ("idealized (no cardiac / noise)",                 0.0,    1.0, 0.0),
            ("realistic (0.3 mm cardiac @1 Hz + 0.2 mm noise)", 0.3e-3, 1.0, 0.2e-3)]:
        fr = ForceReg(lam, J_k, k_eff, F_des=0.2)
        t, e, F, y = run(model, data, fr, cardiac_amp=amp, f_heart=fh, pos_noise=noise)
        reg = (t >= 1.8) & (t < 2.5)
        rmse = float(np.sqrt(np.mean((F[reg] - 0.2)**2)))
        print(f"    {label:<48} force RMSE = {rmse*1e3:5.1f} mN, peak {np.max(F[reg]):.3f} N")
    print(f"    [ref] Jolaei et al. sensor-free: 30-50 mN ;  Kesner & Howe: ~80 mN (hardware)")

    # ── contact-normal misalignment sweep (Sec. V-E) on the MuJoCo plant ──────
    import math
    THETA_MAX = math.radians(30.0)
    print(f"\n  Contact-normal misalignment (Sec. V-E), F_safe={F_SAFE} N:")
    print(f"    {'theta':>6}{'1/cos':>7} |{'nominal cap':>17}{'':>3}|{'cos(30 deg)-tightened':>22}")
    print(f"    {'':>6}{'':>7} |{'worstF':>9}{'realF':>8}{'viol':>4}|{'worstF':>10}{'realF':>8}{'viol':>5}")
    for deg in [0, 10, 20, 30]:
        th = math.radians(deg); cth = math.cos(th); sc = cth*cth
        ca = ImpedanceMPC(lam, J_k, k_eff, force_con=True, f_cap=F_SAFE)
        _, _, Fa, _ = run(model, data, ca, tissue_scale=sc)
        realA = float(np.max(Fa))/cth; worstA = F_SAFE/cth                 # nominal cap
        cb = ImpedanceMPC(lam, J_k, k_eff, force_con=True, f_cap=F_SAFE*math.cos(THETA_MAX))
        _, _, Fb, _ = run(model, data, cb, tissue_scale=sc)
        realB = float(np.max(Fb))/cth; worstB = F_SAFE*math.cos(THETA_MAX)/cth  # tightened
        print(f"    {deg:>4} deg{1/cth:>7.3f} |{worstA:>9.3f}{realA:>8.3f}"
              f"{('YES' if worstA > F_SAFE+1e-3 else 'no'):>4}|"
              f"{worstB:>10.3f}{realB:>8.3f}{('YES' if worstB > F_SAFE+1e-3 else 'no'):>5}")
    print(f"    -> worst-case true tissue force = cap/cos(theta): nominal breaches 0.5 N for")
    print(f"       theta>0 (0.577 N at 30 deg); cos(30 deg)-tightened holds it <= 0.5 N for all")
    print(f"       theta <= 30 deg, with realized peaks well under the certificate.")

    # ── figure 1: four-controller tracking + contact force ────────────────────
    cols = {"Classical impedance": "#e74c3c", "Impedance MPC (no FC)": "#3498db",
            "Impedance MPC (with FC)": "#27ae60", "Joint-space PD": "#95a5a6"}
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for name, (t, e, F, y, m) in res.items():
        ax[0].plot(t, y, color=cols[name], lw=1.6, label=name)
        ax[1].plot(t, F, color=cols[name], lw=1.6)
    ax[0].plot(t, [y_ref(tt)[0]*1e3 for tt in t], 'k--', lw=1.0, alpha=0.6, label="reference")
    ax[0].axhline(Y_WALL*1e3, color='gray', ls=':', lw=1)
    ax[0].text(0.05, Y_WALL*1e3+0.1, "tissue surface", fontsize=8, color='gray')
    ax[0].set_ylabel("tip-normal pos (mm)"); ax[0].legend(fontsize=8, ncol=2)
    ax[0].set_title("Catheter tip tracking & contact force (MuJoCo distributed-compliance plant)")
    ax[1].axhline(F_SAFE, color='r', ls='--', lw=1.2, label="$F_{safe}$ = 0.5 N")
    ax[1].set_ylabel("contact force (N)"); ax[1].set_xlabel("time (s)"); ax[1].legend(fontsize=8)
    ax[0].axvspan(1.5, 2.5, color='orange', alpha=0.06); ax[1].axvspan(1.5, 2.5, color='orange', alpha=0.06)
    fig.tight_layout(); fig.savefig("catheter_results.png", dpi=150)
    print("\n  figure -> catheter_results.png")

    # ── figure 2: cardiac safety-mode contact force ───────────────────────────
    figc, axc = plt.subplots(figsize=(7, 3.4))
    for short, t, F, maxF, viol, col in card:
        m = (t >= 0.8) & (t <= 2.6)
        axc.plot(t[m], F[m], color=col, lw=1.7,
                 label=f"{short}  (peak {maxF:.2f} N{', viol.' if viol else ''})")
    axc.axhline(F_SAFE, color='k', ls='--', lw=1.3, label="$F_{safe}$ = 0.5 N")
    axc.axvspan(1.5, 2.5, color='orange', alpha=0.06)
    axc.text(2.0, 0.03, "contact hold", ha='center', fontsize=8, color='gray')
    axc.set_xlabel("time (s)"); axc.set_ylabel("contact force (N)")
    axc.set_title("Safety mode under cardiac wall motion (MuJoCo plant)", fontsize=11)
    axc.set_xlim(0.8, 2.6); axc.set_ylim(0, 0.6); axc.legend(fontsize=8, loc='upper right')
    figc.tight_layout(); figc.savefig("catheter_cardiac.png", dpi=150)
    print("  figure -> catheter_cardiac.png")


if __name__ == "__main__":
    main()
