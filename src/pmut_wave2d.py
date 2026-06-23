"""
pmut_wave2d.py
==============
2-D finite-difference time-domain (FDTD) acoustic wave simulation of a PMUT
pulse-echo inside the oil-filled steel tube, rendered to an animated MP4/GIF.

Approach follows WaveSimulator2D (github.com/0x23/WaveSimulator2D): an explicit
time-stepped grid solver with a spatial material map and absorbing (sponge)
borders.  We use the *variable-density* staggered velocity-pressure scheme so
reflections are governed by the true acoustic impedance Z = rho*c -- essential
to get the steel<->oil echoes right.

The rendered animation has THREE stacked panels:
  1. the 2-D wave field propagating in the tube,
  2. the PMUT TRANSMITTER drive signal vs time, with a cursor synced to the wave,
  3. the PMUT RECEIVER (measured) signal vs time + noise, with the same cursor,
plus a live read-out of the physical simulation time.

Geometry (longitudinal slice through the tube axis):
    x = axial (3 m oil column + thick steel end block); y = transverse (diameter)
    steel walls top & bottom; near-end 10 mm steel cap with the PMUT behind it;
    far end = thick steel block (near-total reflector).

Honest caveats (also in code):
* This 2-D slice is for VISUALISATION; quantitative SNR lives in pmut_echo_sim.py.
* Oil sound speed (1440 m/s) is exact => round-trip time-of-flight ~4.17 ms is
  physical.  Steel speed is reduced (c_steel_sim) only to relax the CFL step;
  steel DENSITY is kept high so the impedance contrast / reflection is preserved.
* Oil attenuation is applied as a real per-cell loss matching the analytical
  model at the 250 kHz operating point (20 dB/m @ 1 MHz -> 1.25 dB/m @ 250 kHz).
* The receiver-noise level is set for visibility; the true link SNR (see
  pmut_echo_sim.py) is far higher.
"""

from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import imageio.v2 as imageio


# ----------------------------------------------------------------------------
# PARAMETERS
# ----------------------------------------------------------------------------
class P:
    # grid / geometry
    dx = 1.0e-3                 # cell size [m] (>=~6 cells per oil wavelength @250kHz)
    oil_len = 3.0               # oil column length [m]
    outer_d = 300.0e-3          # tube outer diameter [m]
    wall = 10.0e-3              # side wall thickness [m]
    end_block = 60.0e-3         # far-end steel block thickness [m]
    pml = 20                    # sponge layer thickness [cells]

    # media
    c_oil = 1440.0;  rho_oil = 850.0
    c_steel_sim = 2500.0;  rho_steel = 7850.0     # speed reduced for dt; see header

    # oil attenuation (matches analytical model: 20 dB/m AT 250 kHz, ~f^2)
    # => 120 dB round-trip absorption over 6 m: the wave dies in transit.
    att_oil_db_m_ref = 20.0
    att_ref_f = 250.0e3
    att_n = 2.0

    # source / operating point
    f_sim = 250.0e3             # excitation centre frequency [Hz]
    n_cycles = 4.0              # Gaussian-burst width (cycles)
    src_radius = 30.0e-3        # source half-height [m] (visual clarity)

    # receiver noise (set for visibility; see header)
    rx_display_snr_db = 22.0    # echo-to-noise on the displayed receiver trace

    # time / output
    sim_time = 4.7e-3           # total simulated time [s] (> round trip 4.17 ms)
    cfl = 0.5
    n_frames = 220
    gif_path = "wave_tube.gif"
    mp4_path = "wave_tube.mp4"


def build_model(p: P):
    """Construct density / bulk-modulus / damping maps and grid indices."""
    Lx = p.oil_len + p.end_block + p.dx * (p.pml + 2)
    Ly = p.outer_d + p.dx * 2 * (p.pml + 2)
    Nx = int(round(Lx / p.dx))
    Ny = int(round(Ly / p.dx))

    rho = np.full((Nx, Ny), p.rho_oil, dtype=np.float32)
    c = np.full((Nx, Ny), p.c_oil, dtype=np.float32)

    ymid = Ny // 2
    half_in = int(round((p.outer_d / 2 - p.wall) / p.dx))
    half_out = int(round((p.outer_d / 2) / p.dx))
    y0, y1 = ymid - half_in, ymid + half_in
    yw0, yw1 = ymid - half_out, ymid + half_out

    steel = np.zeros((Nx, Ny), dtype=bool)
    steel[:, yw0:y0] = True                                 # bottom wall
    steel[:, y1:yw1] = True                                 # top wall
    x_cap = p.pml + int(round(p.wall / p.dx))
    steel[p.pml:x_cap, yw0:yw1] = True                      # near cap
    x_oilend = x_cap + int(round(p.oil_len / p.dx))
    steel[x_oilend:x_oilend + int(round(p.end_block / p.dx)), yw0:yw1] = True  # far block

    rho[steel] = p.rho_steel
    c[steel] = p.c_steel_sim
    K = rho * c ** 2

    # damping field = sponge borders * oil bulk attenuation
    damp = np.ones((Nx, Ny), dtype=np.float32)
    ramp = np.linspace(0.0, 1.0, p.pml, dtype=np.float32)
    d_edge = 0.06
    for i in range(p.pml):
        f = d_edge * (1 - ramp[i]) ** 2
        damp[i, :] = np.minimum(damp[i, :], 1 - f)
        damp[-1 - i, :] = np.minimum(damp[-1 - i, :], 1 - f)
        damp[:, i] = np.minimum(damp[:, i], 1 - f)
        damp[:, -1 - i] = np.minimum(damp[:, -1 - i], 1 - f)

    return dict(Nx=Nx, Ny=Ny, rho=rho, c=c, K=K, steel=steel, damp=damp,
                x_cap=x_cap, x_oilend=x_oilend, y0=y0, y1=y1, yw0=yw0, yw1=yw1,
                ymid=ymid)


def run(p: P):
    m = build_model(p)
    Nx, Ny = m["Nx"], m["Ny"]
    rho, K = m["rho"], m["K"]
    steel = m["steel"]
    cmax = float(m["c"].max())
    dt = p.cfl * p.dx / (cmax * np.sqrt(2.0))
    nsteps = int(round(p.sim_time / dt))

    # oil attenuation at the operating frequency -> per-step decay in oil cells
    alpha_db_m = p.att_oil_db_m_ref * (p.f_sim / p.att_ref_f) ** p.att_n
    alpha_np = alpha_db_m / 8.686
    oil_decay = float(np.exp(-alpha_np * p.c_oil * dt))
    damp = m["damp"].copy()
    damp[~steel] *= oil_decay                               # bulk oil loss

    print(f"  grid {Nx}x{Ny} = {Nx*Ny/1e3:.0f}k cells | dt={dt*1e9:.1f} ns | "
          f"{nsteps} steps | sim {p.sim_time*1e3:.2f} ms")
    print(f"  oil attenuation @ {p.f_sim/1e3:.0f} kHz = {alpha_db_m:.2f} dB/m "
          f"(=> {alpha_db_m*2*p.oil_len:.1f} dB round trip)")

    p_f = np.zeros((Nx, Ny), dtype=np.float32)
    vx = np.zeros((Nx - 1, Ny), dtype=np.float32)
    vy = np.zeros((Nx, Ny - 1), dtype=np.float32)

    inv_rho_x = (dt / p.dx) / (0.5 * (rho[1:, :] + rho[:-1, :]))
    inv_rho_y = (dt / p.dx) / (0.5 * (rho[:, 1:] + rho[:, :-1]))
    div = K * dt
    damp_vx = 0.5 * (damp[1:, :] + damp[:-1, :])
    damp_vy = 0.5 * (damp[:, 1:] + damp[:, :-1])

    sx = m["x_cap"] + 2
    srad = int(round(p.src_radius / p.dx))
    sy0, sy1 = m["ymid"] - srad, m["ymid"] + srad
    t0 = p.n_cycles / p.f_sim
    tau = p.n_cycles / (2.0 * p.f_sim)
    src_end = 2 * t0 + 4 * tau

    frame_every = max(1, nsteps // p.n_frames)
    steel_mask = steel.T[::-1]
    field_frames, frame_steps = [], []
    tx_sig = np.zeros(nsteps, np.float32)
    rx_sig = np.zeros(nsteps, np.float32)

    for n in range(nsteps):
        t = n * dt
        vx -= inv_rho_x * (p_f[1:, :] - p_f[:-1, :])
        vy -= inv_rho_y * (p_f[:, 1:] - p_f[:, :-1])
        dvx = np.zeros_like(p_f); dvy = np.zeros_like(p_f)
        dvx[1:-1, :] = (vx[1:, :] - vx[:-1, :]) / p.dx
        dvy[:, 1:-1] = (vy[:, 1:] - vy[:, :-1]) / p.dx
        p_f -= div * (dvx + dvy)
        s = 0.0
        if t < src_end:
            s = np.sin(2 * np.pi * p.f_sim * (t - t0)) * np.exp(-((t - t0) / tau) ** 2)
            p_f[sx, sy0:sy1] += np.float32(s)
        p_f *= damp; vx *= damp_vx; vy *= damp_vy

        tx_sig[n] = s
        rx_sig[n] = p_f[sx, sy0:sy1].mean()

        if n % frame_every == 0:
            field_frames.append(field_to_rgb(p_f, steel_mask))
            frame_steps.append(n)
            if n % (frame_every * 22) == 0:
                print(f"    step {n}/{nsteps}  t={t*1e3:.2f} ms")

    t_axis = np.arange(nsteps) * dt
    return dict(field=field_frames, fstep=np.array(frame_steps),
                tx=tx_sig, rx=rx_sig, t=t_axis, dt=dt, nsteps=nsteps, m=m)


def field_to_rgb(p_f, steel_mask):
    img = p_f.T[::-1]
    norm = np.clip(img / 0.06, -1, 1) * 0.5 + 0.5
    rgb = (cm.RdBu_r(norm)[..., :3] * 255).astype(np.uint8)
    rgb[steel_mask] = np.array([90, 90, 95], np.uint8)
    rgb = rgb[:, ::3, :]                                    # axial downsample
    return rgb


def compose(out, p: P):
    """Build the 3-panel synced animation and write MP4 + GIF.

    Wave-field panel: the real FDTD propagation (with true attenuation).
    Receiver panel: driven by the CALIBRATED, reciprocity-consistent 1-D link
    budget (pmut_echo_sim), which is passive (no gain>1) and gives an absolute
    echo voltage + modeled receiver noise floor + absolute SNR at every
    frequency.  We plot the echo burst at the round-trip time on top of the
    modeled noise floor (transmit leakage clipped off-scale), so the panel shows
    the true SNR -- clean above 135 kHz... below it, buried above it.
    """
    from pmut_echo_sim import Config, link_budget

    field = out["field"]; fstep = out["fstep"]
    t = out["t"]; tx = out["tx"]
    tof = 2 * p.oil_len / p.c_oil
    tms = t * 1e3
    f = p.f_sim
    Vin = 10.0
    tx_v = tx / (np.abs(tx).max() + 1e-12) * Vin            # transmitter drive [V]
    rng = np.random.default_rng(0)
    tb = p.n_cycles / f

    # calibrated link budget at this frequency
    cfg = Config(f0=f, geometry="tube", att_oil_db_m=p.att_oil_db_m_ref,
                 att_ref_f=p.att_ref_f, att_n=p.att_n)
    lb = link_budget(cfg)
    echo_amp = lb["V_echo"] * 1e6                           # echo [µV]
    noise = lb["V_noise"] * 1e6                             # modeled noise floor [µV]
    snr = lb["snr_model"]
    verdict = "DETECTABLE" if snr >= 12 else "buried in noise"

    # receiver trace in absolute µV: noise + echo burst at TOF + transmit leakage
    rx_uv = rng.normal(0, noise, size=t.shape).astype(np.float64)
    te = t - tof
    rx_uv += np.where((te >= 0) & (te < tb),
                      np.sin(2*np.pi*f*te) * np.sin(np.pi*np.clip(te/tb,0,1))**2, 0.0) * echo_amp
    rx_uv += np.where(t < tb,
                      np.sin(2*np.pi*f*t) * np.sin(np.pi*np.clip(t/tb,0,1))**2, 0.0) * max(noise, echo_amp) * 200
    # NORMALISE to the noise floor so 100/135/200 kHz share ONE fixed y-scale:
    # noise -> +/-1 band, the 12 dB detection threshold -> fixed +/-4 lines.
    rx_sig = rx_uv / noise
    rx_unit = "RX / noise floor  (×)"
    noise_lbl = f"±1 (= {noise/1e3:.1f} mV modeled)"
    thr = 10 ** (12 / 20)                                   # 12 dB in amplitude = ~3.98
    tx_lim = 1.1 * Vin
    rx_lim = 6.0                                            # fixed across all clips

    # ---- persistent figure; update image/cursors per frame ------------------
    fig = plt.figure(figsize=(11, 7.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.0, 1.0, 1.0], hspace=0.55)
    axw = fig.add_subplot(gs[0]); axt = fig.add_subplot(gs[1]); axr = fig.add_subplot(gs[2])

    im = axw.imshow(field[0], aspect="auto",
                    extent=[0, p.oil_len + p.end_block, -p.outer_d / 2 * 1e3, p.outer_d / 2 * 1e3])
    axw.set_title(f"Acoustic wave in oil-filled steel tube @ {f/1e3:.0f} kHz "
                  f"(SNR {snr:+.0f} dB, {verdict})  —  red/blue = pressure, gray = steel",
                  fontsize=10)
    axw.set_xlabel("axial position [m]"); axw.set_ylabel("radius [mm]")
    time_txt = axw.text(0.015, 0.88, "", transform=axw.transAxes, fontsize=11,
                        color="k", fontweight="bold",
                        bbox=dict(fc="white", ec="0.5", alpha=0.85, pad=3))

    axt.plot(tms, tx_v, color="C1", lw=0.8)
    axt.set_ylim(-tx_lim, tx_lim); axt.set_xlim(0, tms[-1])
    axt.set_ylabel("drive [V]")
    axt.set_title(f"PMUT transmitter drive signal  ({f/1e3:.0f} kHz burst, {Vin:.0f} V)", fontsize=9)
    axt.grid(alpha=0.3)
    cur_t = axt.axvline(0, color="k", lw=1.4)

    axr.plot(tms, rx_sig, color="C0", lw=0.5, label="measured (echo + noise)")
    axr.axhspan(-1, 1, color="0.5", alpha=0.18, label=f"noise floor {noise_lbl}")
    axr.axhline(thr, color="green", ls="--", lw=1.0)
    axr.axhline(-thr, color="green", ls="--", lw=1.0, label="±12 dB threshold")
    axr.axvline(tof * 1e3, color="red", ls="--", lw=1.1, label=f"echo @ {tof*1e3:.2f} ms")
    axr.set_ylim(-rx_lim, rx_lim); axr.set_xlim(0, tms[-1])
    axr.set_ylabel(rx_unit); axr.set_xlabel("time [ms]")
    axr.set_title(f"PMUT receiver signal  (echo SNR = {snr:+.1f} dB — {verdict}; "
                  f"shared scale, transmit leakage off-scale)", fontsize=9)
    axr.grid(alpha=0.3); axr.legend(fontsize=7, loc="upper left", ncol=2)
    cur_r = axr.axvline(0, color="k", lw=1.4)

    fig.canvas.draw()
    frames_rgb = []
    for k, fs in enumerate(fstep):
        tt = t[fs]
        im.set_data(field[k])
        cur_t.set_xdata([tt * 1e3, tt * 1e3])
        cur_r.set_xdata([tt * 1e3, tt * 1e3])
        time_txt.set_text(f"sim time  t = {tt*1e3:6.3f} ms   ({100*tt/t[-1]:4.1f} %)")
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
        buf = buf[:buf.shape[0] // 2 * 2, :buf.shape[1] // 2 * 2, :]
        frames_rgb.append(buf)
    plt.close(fig)

    print(f"  writing {len(frames_rgb)} frames -> {p.mp4_path} / {p.gif_path}")
    w = imageio.get_writer(p.mp4_path, fps=30, codec="libx264", quality=8,
                           macro_block_size=None)
    for fr in frames_rgb:
        w.append_data(fr)
    w.close()
    imageio.mimsave(p.gif_path, frames_rgb[::2], fps=15, loop=0)   # lighter gif


if __name__ == "__main__":
    import sys
    p = P()
    # optional CLI: frequency in kHz (default 250).  Grid auto-set to ~9 cells
    # per oil wavelength; outputs tagged by frequency.
    f_khz = float(sys.argv[1]) if len(sys.argv) > 1 else 250.0
    p.f_sim = f_khz * 1e3
    dx = p.c_oil / p.f_sim / 9.0
    p.dx = float(min(3.0e-3, max(1.0e-3, round(dx / 0.5e-3) * 0.5e-3)))
    tag = f"{int(round(f_khz)):03d}khz"
    p.gif_path = f"{tag}.gif"
    p.mp4_path = f"{tag}.mp4"
    print(f"2-D FDTD acoustic wave simulation (oil-filled steel tube) @ {f_khz:.0f} kHz "
          f"(dx={p.dx*1e3:.1f} mm)")
    out = run(p)
    compose(out, p)
    print("done.")
