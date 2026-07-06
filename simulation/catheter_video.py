#!/usr/bin/env python3
"""
Benchmark-simulation VIDEO — force-constrained Impedance MPC on the MuJoCo plant
================================================================================
Renders the headline benchmark task (approach -> press -> hold -> retract against a
Kelvin-Voigt tissue wall) running the force-constrained Impedance MPC on the
distributed-compliance MuJoCo PRB catheter of `catheter_mujoco.py` -- the same plant
and controller as `catheter_mpc_mujoco.py` (Table tab:mujoco). Writes an mp4 with the
MuJoCo scene on the left and a live contact-force / tip-position panel on the right,
showing the hard F_safe=0.5 N bound being respected throughout the press.

This is the renderable, physics-engine version of the benchmark: the RK4 scalar plant
(`catheter_benchmark.py`) has no geometry to render, so the video is produced on the
MuJoCo plant, which carries the catheter's distributed bending and tendon transmission.

Run: python3 catheter_video.py   (writes catheter_benchmark.mp4)
"""
import numpy as np
import mujoco
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
import imageio.v2 as imageio
import catheter_mujoco as C
from catheter_mpc_mujoco import y_ref, measure_plant, Y_WALL, F_SAFE, DT_CTRL

# ── video / render settings ───────────────────────────────────────────────────
H = W = 480                 # mujoco frame size (square; default offscreen framebuffer)
FPS = 30                    # output frame rate
FRAME_EVERY = 8             # capture one frame per N control steps (2 ms each)
T_END = 3.5
OUT = "catheter_benchmark.mp4"


def make_camera():
    cam = mujoco.MjvCamera()
    cam.type = mujoco.mjtCamera.mjCAMERA_FREE
    cam.lookat[:] = [0.022, 0.0, -0.004]   # center on the bending segment + wall
    cam.distance = 0.085
    cam.azimuth = 90.0                      # view the x-z bending plane from +y
    cam.elevation = -8.0
    return cam


STAGES = ["approach", "press", "hold", "retract"]


def annotate_scene(scene, t, phase):
    """Overlay component names on the rendered MuJoCo scene (camera is fixed, so the
    anchor pixels are stable). Labels the catheter links, tendon, tip, tissue wall, base,
    and draws the task-stage progression bar (approach->press->hold->retract) at the top
    with the current stage highlighted."""
    h, w = scene.shape[:2]
    fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(scene); ax.axis("off")
    ax.set_xlim(0, w); ax.set_ylim(h, 0)

    # task-stage progression bar (top of the scene)
    ax.text(8, 16, "stage:", fontsize=7.5, color="white", weight="bold", va="center")
    centers = [78, 150, 210, 280]
    for i, st in enumerate(STAGES):
        active = (st == phase)
        ax.text(centers[i], 16, st, fontsize=7.5, ha="center", va="center",
                color="black" if active else "#c8c8c8", weight="bold",
                bbox=dict(boxstyle="round,pad=0.25", lw=0.5,
                          fc="#e67e22" if active else (0, 0, 0, 0.6),
                          ec="white" if active else "#888888"))
        if i < len(STAGES) - 1:
            ax.text((centers[i] + centers[i + 1]) / 2, 16, "›", fontsize=9,
                    ha="center", va="center", color="white")
    ax.text(w - 8, 16, f"t = {t:4.2f} s", fontsize=7.5, color="white", weight="bold",
            ha="right", va="center")

    # component names
    tp = dict(fontsize=7.5, color="white", weight="bold", ha="left", va="center",
              bbox=dict(boxstyle="round,pad=0.25", fc=(0, 0, 0, 0.68), ec="white", lw=0.5))
    ar = dict(arrowstyle="->", color="white", lw=1.1, shrinkA=0, shrinkB=2)
    ax.annotate("catheter (8 PRB links)", xy=(195, 200), xytext=(60, 58), arrowprops=ar, **tp)
    ax.annotate("tendon (actuator)",      xy=(168, 221), xytext=(6, 292),  arrowprops=ar, **tp)
    ax.annotate("tip (contact sphere)",   xy=(420, 230), xytext=(250, 112), arrowprops=ar, **tp)
    ax.annotate("tissue wall (Kelvin-Voigt)", xy=(245, 315), xytext=(110, 446), arrowprops=ar, **tp)
    ax.annotate("base (fixed)",           xy=(98, 208),  xytext=(4, 150),  arrowprops=ar, **tp)
    canvas = FigureCanvasAgg(fig); canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[:, :, :3].copy()
    plt.close(fig)
    return buf


def panel(times, refs, ys, Fs, t_now, phase):
    """Right-hand live panel: contact force (top) + tip-normal pos vs ref (bottom)."""
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
    ax1 = fig.add_axes([0.17, 0.57, 0.78, 0.35])
    ax2 = fig.add_axes([0.17, 0.11, 0.78, 0.35])

    ax1.axhline(F_SAFE, color="r", ls="--", lw=1.4, label="$F_{safe}$ = 0.5 N")
    ax1.plot(times, Fs, color="#27ae60", lw=2.0)
    if Fs:
        ax1.plot(times[-1], Fs[-1], "o", color="#27ae60", ms=5)
    ax1.set_xlim(0, T_END); ax1.set_ylim(0, 0.6)
    ax1.set_ylabel("contact force (N)")
    ax1.legend(loc="upper left", fontsize=8)
    ax1.set_title("Force-constrained Impedance MPC (MuJoCo plant)", fontsize=10)

    ax2.axhline(Y_WALL * 1e3, color="black", ls=":", lw=1.2)
    ax2.text(0.05, Y_WALL * 1e3 + 0.05, "tissue surface", fontsize=7, color="black")
    ax2.plot(times, refs, "k--", lw=1.0, label="reference")
    ax2.plot(times, ys, color="#2980b9", lw=2.0, label="tip")
    if ys:
        ax2.plot(times[-1], ys[-1], "o", color="#2980b9", ms=5)
    ax2.set_xlim(0, T_END); ax2.set_ylim(-0.3, 6.2)
    ax2.set_xlabel("time (s)"); ax2.set_ylabel("tip-normal pos (mm)")
    ax2.legend(loc="upper left", fontsize=8)

    fig.text(0.17, 0.965, f"t = {t_now:4.2f} s", fontsize=11, weight="bold")
    fig.text(0.55, 0.965, phase, fontsize=11, color="#b9770e", weight="bold")

    canvas = FigureCanvasAgg(fig); canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[:, :, :3].copy()
    plt.close(fig)
    return buf


def phase_of(t):
    if t < 1.0:   return "approach"
    if t < 1.5:   return "press"
    if t < 2.5:   return "hold"
    if t < 3.5:   return "retract"
    return "done"


def main():
    model = mujoco.MjModel.from_xml_string(C.build_xml(d_tend=0.0025, tendon_max=8.0))
    data = mujoco.MjData(model)
    lam, J_k, k_eff = measure_plant(model, data)
    print(f"plant: Lambda={lam:.3e}  J_k={J_k:.3f}  k_eff={k_eff:.1f} N/m  F_safe={F_SAFE} N")

    Kd, Dd, Ki = 900.0 * lam, 60.0 * lam, 300.0
    n_sub = int(round(DT_CTRL / model.opt.timestep))
    nctrl = int(T_END / DT_CTRL)

    mujoco.mj_resetData(model, data)
    renderer = mujoco.Renderer(model, height=H, width=W)
    cam = make_camera()

    times, refs, ys, Fs = [], [], [], []
    frames = []
    dhat = 0.0
    peakF = 0.0

    for k in range(nctrl):
        t = k * DT_CTRL
        y, yd = C.tip_state(model, data)
        yr, yrd = y_ref(t)
        e, de = yr - y, yrd - yd
        dhat += Ki * e * DT_CTRL
        F_mpc = float(np.clip(Kd * e + Dd * de + dhat, -F_SAFE, F_SAFE))   # hard bound
        ff = k_eff * yr
        T = float(np.clip(-(ff + F_mpc) / J_k, -C.TENDON_MAX, 0.0))        # J_k remap
        data.ctrl[0] = T
        Fc = 0.0
        for _ in range(n_sub):
            data.ctrl[1] = 0.0                                            # static wall
            f, _ = C.tissue_force(model, data); C.apply_tissue(model, data, f)
            mujoco.mj_step(model, data)
            Fc = max(Fc, f)
        ytn, _ = C.tip_state(model, data)
        peakF = max(peakF, Fc)
        times.append(t); refs.append(yr * 1e3); ys.append(ytn * 1e3); Fs.append(Fc)

        if k % FRAME_EVERY == 0:
            renderer.update_scene(data, camera=cam)
            scene = annotate_scene(renderer.render(), t, phase_of(t))
            pan = panel(times, refs, ys, Fs, t, phase_of(t))
            frames.append(np.hstack([scene, pan]))

    renderer.close()
    print(f"peak contact force over run = {peakF:.3f} N  "
          f"({'VIOLATES' if peakF > F_SAFE + 1e-3 else 'within'} F_safe={F_SAFE} N)")
    print(f"writing {len(frames)} frames -> {OUT} at {FPS} fps")
    imageio.mimsave(OUT, frames, fps=FPS, codec="libx264", quality=8,
                    macro_block_size=8)
    print(f"done -> {OUT}")


if __name__ == "__main__":
    main()
