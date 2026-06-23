"""
ranging_precision.py
Closes the loop from SNR to ranging precision sigma_d(f) for single-reflector
pulse-echo (the hydraulic-cylinder piston face is one strong specular echo).

Classical axial resolution ~ c/(2B) with B ~ f0/Q is poor at low frequency, but
that measures separating two reflectors. For ONE isolated strong echo the relevant
figure is the Cramer-Rao bound on time-of-flight, reachable by carrier-phase (I/Q)
estimation. We verify the phase-precision formula by Monte-Carlo: synthesize the
modeled echo + modeled receiver noise, estimate the carrier phase over the known
burst window with the maximum-likelihood single-bin estimator, and convert phase
jitter to a one-way distance error via  delta_d = c/(4 pi f) * sigma_phi.
"""
import dataclasses
import numpy as np
from pmut_echo_sim import Config, link_budget


def precision(cfg, f, n_trials=20000, n_cycles=10, seed=0):
    lb = link_budget(dataclasses.replace(cfg, f0=f), f=f)
    A = lb["V_echo"]            # echo voltage amplitude [V]
    sigma = lb["V_noise"]       # modeled receiver noise rms [V]
    c = cfg.c_oil
    fs = 40.0 * f               # generous oversampling
    tb = n_cycles / f
    t = np.arange(0, tb, 1.0 / fs)
    env = np.sin(np.pi * t / tb) ** 2                 # Hann^2 burst envelope
    s = A * env * np.cos(2 * np.pi * f * t)           # received echo (true phase 0)
    ref = env * np.exp(-1j * 2 * np.pi * f * t)       # ML matched reference
    rng = np.random.default_rng(seed)

    ph = np.empty(n_trials)
    for k in range(n_trials):
        r = s + rng.normal(0, sigma, t.size)
        ph[k] = np.angle(np.sum(r * ref))             # single-bin phase estimate
    # circular standard deviation (handles phase wrap at low SNR)
    R = np.abs(np.mean(np.exp(1j * ph)))
    sphi_emp = np.sqrt(-2.0 * np.log(max(R, 1e-12)))
    sd_emp = c / (4 * np.pi * f) * sphi_emp

    snr_lin = 10 ** (lb["snr_model"] / 10)
    sphi_th = 1.0 / np.sqrt(2 * snr_lin)              # high-SNR CRLB
    sd_th = c / (4 * np.pi * f) * sphi_th
    lam = c / f
    return dict(f=f, snr_db=lb["snr_model"], snr_lin=snr_lin,
                lam_mm=lam * 1e3, sphi_th=sphi_th, sphi_emp=sphi_emp,
                sd_th_mm=sd_th * 1e3, sd_emp_mm=sd_emp * 1e3)


if __name__ == "__main__":
    cfg = Config()
    B = cfg.f0 / cfg.Q
    print(f"resonant bandwidth B = f0/Q (at 135 kHz) ~ {135e3/cfg.Q/1e3:.2f} kHz")
    print(f"classical axial resolution c/(2B) at 135 kHz ~ "
          f"{cfg.c_oil/(2*135e3/cfg.Q)*1e3:.0f} mm\n")
    hdr = f"{'f[kHz]':>7} {'SNR[dB]':>8} {'lam[mm]':>8} {'sig_phi[rad]':>12} " \
          f"{'sd_theory[mm]':>14} {'sd_MC[mm]':>11}"
    print(hdr); print("-" * len(hdr))
    for fk in [75, 100, 135, 200, 250]:
        r = precision(cfg, fk * 1e3)
        print(f"{fk:7d} {r['snr_db']:8.1f} {r['lam_mm']:8.2f} "
              f"{r['sphi_emp']:12.4f} {r['sd_th_mm']:14.4f} {r['sd_emp_mm']:11.4f}")
