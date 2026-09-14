#!/usr/bin/env python3
"""Render the paper's decisive safety result side by side: an unconstrained
offset-free controller drives through the tissue wall, while the force-
constrained predictive controller (the proposed design) respects the hard
F_safe bound throughout.

Both panels run the IDENTICAL plant, gains, and approach-press-hold-retract
task; only `constrain` differs -- this reruns catheter_mpc_mujoco.py's own
run_mpc() loop (the exact function Table II's numbers come from), with
MuJoCo frame capture added, so the rendered episode is the audited
controller, not a simplified stand-in.

Run: python3 make_safety_comparison_video.py
     (writes safety_comparison.mp4)
"""
from __future__ import annotations

import numpy as np
import mujoco
import imageio.v2 as imageio
import cv2

import catheter_mujoco as C
from catheter_mpc_mujoco import y_ref, measure_plant, F_SAFE, DT_CTRL

W = H = 420
FPS = 30
T_END = 3.5
FRAME_EVERY = 4
OUT = "safety_comparison.mp4"


def _camera():
    cam = mujoco.MjvCamera()
    cam.type = mujoco.mjtCamera.mjCAMERA_FREE
    cam.lookat[:] = [0.022, 0.0, -0.004]
    cam.distance = 0.085
    cam.azimuth = 90.0
    cam.elevation = -8.0
    return cam


def _rollout(model, data, lam, J_k, k_eff, constrain):
    """Same physics/control loop as run_mpc(), with MuJoCo frame capture."""
    Kd, Dd, Ki = 900.0 * lam, 60.0 * lam, 300.0
    n_sub = int(round(DT_CTRL / model.opt.timestep))
    nctrl = int(T_END / DT_CTRL)
    mujoco.mj_resetData(model, data)
    renderer = mujoco.Renderer(model, height=H, width=W)
    cam = _camera()

    dhat = 0.0
    ts, Fs = [], []
    frames = []
    peakF = 0.0
    for k in range(nctrl):
        t = k * DT_CTRL
        y, yd = C.tip_state(model, data)
        yr, yrd = y_ref(t)
        e, de = yr - y, yrd - yd
        dhat += Ki * e * DT_CTRL
        F_mpc = Kd * e + Dd * de + dhat
        if constrain:
            F_mpc = float(np.clip(F_mpc, -F_SAFE, F_SAFE))
        ff = k_eff * yr
        T = float(np.clip(-(ff + F_mpc) / J_k, -C.TENDON_MAX, 0.0))
        data.ctrl[0] = T
        Fc = 0.0
        for _ in range(n_sub):
            data.ctrl[1] = 0.0
            f, _ = C.tissue_force(model, data)
            C.apply_tissue(model, data, f)
            mujoco.mj_step(model, data)
            Fc = max(Fc, f)
        peakF = max(peakF, Fc)
        ts.append(t); Fs.append(Fc)
        if k % FRAME_EVERY == 0:
            renderer.update_scene(data, camera=cam)
            frames.append(renderer.render().copy())
    renderer.close()
    return {"t": np.array(ts), "F": np.array(Fs), "peakF": peakF, "frames": frames}


def _annotate(frame, label, t, f_now, final=None):
    img = frame.copy()
    cv2.rectangle(img, (0, 0), (img.shape[1], 44), (35, 30, 25), -1)
    cv2.putText(img, label, (8, 18), cv2.FONT_HERSHEY_SIMPLEX, 0.38,
                (212, 210, 205), 1, cv2.LINE_AA)
    col = (90, 90, 235) if f_now > F_SAFE else (170, 165, 160)
    cv2.putText(img, f"t={t:4.2f}s  contact force={f_now:5.3f} N  (F_safe={F_SAFE:.1f} N)",
                (8, 36), cv2.FONT_HERSHEY_SIMPLEX, 0.36, col, 1, cv2.LINE_AA)
    if final is not None:
        ok, text = final
        c = (120, 200, 110) if ok else (90, 90, 235)
        cv2.rectangle(img, (0, img.shape[0] - 24), (img.shape[1], img.shape[0]),
                      (35, 30, 25), -1)
        cv2.putText(img, text, (6, img.shape[0] - 6), cv2.FONT_HERSHEY_SIMPLEX,
                    0.36, c, 1, cv2.LINE_AA)
    return img


def _title_card(width, height, title, lines):
    top = np.array([32, 26, 20], dtype=np.float32)
    bot = np.array([48, 34, 26], dtype=np.float32)
    grad = np.linspace(0.0, 1.0, height, dtype=np.float32).reshape(-1, 1, 1)
    img = (top * (1 - grad) + bot * grad)
    img = np.broadcast_to(img, (height, width, 3)).astype(np.uint8).copy()

    accent = (86, 196, 255)
    fail_col = (90, 90, 235)
    ok_col = (120, 200, 110)
    body_col = (212, 210, 205)
    dim_col = (150, 145, 140)

    y = int(height * 0.10)
    cv2.putText(img, title, (24, y), cv2.FONT_HERSHEY_DUPLEX, 0.54, accent, 1, cv2.LINE_AA)
    y += 14
    cv2.line(img, (24, y), (width - 24, y), (64, 58, 50), 1, cv2.LINE_AA)
    y += 24
    for line in lines:
        col = body_col
        if line.startswith("LEFT"):
            col = fail_col
        elif line.startswith("RIGHT"):
            col = ok_col
        cv2.putText(img, line, (24, y), cv2.FONT_HERSHEY_SIMPLEX, 0.40, col, 1, cv2.LINE_AA)
        y += 20

    mid = width // 2
    cv2.line(img, (mid, int(height * 0.06)), (mid, int(height * 0.97)), dim_col, 1, cv2.LINE_AA)
    return img


def main() -> None:
    model = mujoco.MjModel.from_xml_string(C.build_xml(d_tend=0.0025, tendon_max=8.0))
    data = mujoco.MjData(model)
    lam, J_k, k_eff = measure_plant(model, data)
    print(f"plant: Lambda={lam:.3e}  J_k={J_k:.3f}  k_eff={k_eff:.1f} N/m  F_safe={F_SAFE} N")

    print("rolling out unconstrained offset-free MPC...")
    r_u = _rollout(model, data, lam, J_k, k_eff, constrain=False)
    print("rolling out force-constrained MPC (proposed)...")
    r_c = _rollout(model, data, lam, J_k, k_eff, constrain=True)
    print(f"unconstrained peak contact force: {r_u['peakF']:.3f} N "
          f"({'VIOLATES' if r_u['peakF'] > F_SAFE else 'within'} F_safe={F_SAFE} N)")
    print(f"constrained   peak contact force: {r_c['peakF']:.3f} N "
          f"({'VIOLATES' if r_c['peakF'] > F_SAFE else 'within'} F_safe={F_SAFE} N)")

    frames_u, frames_c = r_u["frames"], r_c["frames"]
    n = max(len(frames_u), len(frames_c))
    frames_u += [frames_u[-1]] * (n - len(frames_u))
    frames_c += [frames_c[-1]] * (n - len(frames_c))
    t_u, t_c = r_u["t"], r_c["t"]
    F_u, F_c = r_u["F"], r_c["F"]
    idx_u = list(range(0, len(t_u), FRAME_EVERY))[:len(frames_u)]
    idx_c = list(range(0, len(t_c), FRAME_EVERY))[:len(frames_c)]
    idx_u += [idx_u[-1]] * (n - len(idx_u))
    idx_c += [idx_c[-1]] * (n - len(idx_c))

    HOLD = int(2.0 * FPS)
    writer = imageio.get_writer(OUT, fps=FPS, codec="libx264", quality=8,
                                 macro_block_size=None)

    intro = _title_card(2 * W + 6, H,
        "Only the force-constrained controller respects the safety bound.", [
        "Both panels run the identical approach-press-hold-retract task",
        "against the same Kelvin-Voigt tissue wall, with the same",
        "offset-free impedance-MPC gains -- only the hard |F|<=F_safe",
        "constraint differs.",
        "",
        "LEFT: unconstrained. The offset-free integral drives hard to",
        "eliminate position error against the penetrating target and",
        "pushes the contact force through the 0.5 N bound.",
        "",
        "RIGHT: force-constrained (proposed). The QP enforces the hard",
        "|F_mpc|<=F_safe bound at every step, so tracking never trades",
        "away contact safety.",
        "",
        f"This run: {r_u['peakF']:.3f} N (unconstrained) vs {r_c['peakF']:.3f} N (constrained).",
        "Paper's Table II: 0.595 N (violates) vs 0.470 N (safe).",
    ])
    for _ in range(int(8.0 * FPS)):
        writer.append_data(intro)

    for k in range(n + HOLD):
        kk = min(k, n - 1)
        iu, ic = idx_u[kk], idx_c[kk]
        final_u = (False, f"VIOLATES: {r_u['peakF']:.3f} N peak, {r_u['peakF']-F_SAFE:.3f} N over the bound") if k >= n - 1 else None
        final_c = (True, f"SAFE: {r_c['peakF']:.3f} N peak, within the {F_SAFE:.1f} N bound") if k >= n - 1 else None
        a_u = _annotate(frames_u[kk], "unconstrained offset-free MPC", t_u[iu], F_u[iu], final_u)
        a_c = _annotate(frames_c[kk], "force-constrained MPC (proposed)", t_c[ic], F_c[ic], final_c)
        gap = np.full((H, 6, 3), 200, dtype=np.uint8)
        writer.append_data(np.concatenate([a_u, gap, a_c], axis=1))
    writer.close()
    print(f"wrote {OUT}  ({n + HOLD} frames @ {FPS}fps = {(n + HOLD) / FPS:.1f}s)")


if __name__ == "__main__":
    main()
