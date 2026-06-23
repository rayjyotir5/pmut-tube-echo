"""
pmut_echo_sim.py
================
Pulse-echo link-budget + SNR simulation for a PMUT transmitting through a
10 mm steel wall into synthetic mineral oil, ranging a target at 3 m.

WHY THIS EXISTS
---------------
The reference paper (Dangi & Pratap, "System level modeling and design maps of
PMUTs with residual stresses", Sensors & Actuators A 262 (2017) 18-28) models
the TRANSDUCER only: its transmit pressure per volt, receive charge per pascal,
quality factor Q, and the near-invariant mode-shape participation factors
Lambda', Lambda1, Lambda2.  It validates radiation into AIR out to ~0.5 m.

It says NOTHING about:
  * transmission through a layered steel/oil medium,
  * propagation + beam spreading to 3 m,
  * pulse-echo (round trip) off a target,
  * any noise model or SNR.

This script keeps the paper's transducer relations and adds the missing ACOUSTIC
CHANNEL so we can answer: at 1 MHz, through 10 mm steel into mineral oil, does a
target at 3 m return an echo with SNR >= 12 dB?

MODEL OVERVIEW (1-D, normal incidence, on-axis)
-----------------------------------------------
  PMUT  ->  10 mm steel wall  ->  oil (3 m)  ->  target  ->  oil (3 m)  ->  wall  ->  PMUT

  Transmit  : paper transducer -> face pressure into steel
  Wall      : single-layer acoustic transfer-matrix (intensity transmission tau),
              handles the steel thickness resonances (n*c/2d windows)
  Oil path  : on-axis baffled-piston far-field spreading g(z)=pi*a^2/(lambda*z),
              plus frequency-dependent attenuation; monostatic flat-target uses
              reciprocal aperture so the round-trip geometric factor is g(z)^2.
  Target    : pressure reflection coefficient Gamma
  Receive   : paper charge sensitivity Qrx/Pin (Eq 17) -> charge-amp voltage
  SNR       : swept against a range of receiver noise floors; >=12 dB region shaded

Everything physical is in CONFIG below and documented.  All channel additions
are clearly separated from the paper-derived transducer block.
"""

from __future__ import annotations

import dataclasses
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------------
# CONFIGURATION  (every number here is an explicit, documented assumption)
# ----------------------------------------------------------------------------
@dataclasses.dataclass
class Config:
    # --- Operating point -----------------------------------------------------
    # NOTE: 1 MHz is NOT assumed optimal.  In the oil-filled tube the SNR-optimal
    # frequency is far lower (oil attenuation ~f^2 AND receive sensitivity ~1/f^2
    # both reward low f); default re-centred to a practical value, swept below.
    f0: float = 250.0e3        # centre frequency [Hz]
    Vin: float = 10.0          # transmit drive amplitude [V]
    range_m: float = 3.0       # one-way distance to reflector in oil [m]

    # --- Geometry ------------------------------------------------------------
    # 'tube' : oil-filled steel tube waveguide (this task).  PMUT at the near end
    #          (behind the 10 mm steel end-cap), far end is a thick steel
    #          reflector (~total reflection).  The wave is laterally confined, so
    #          free-field spherical spreading is replaced by a bounded cross-
    #          section "fill" factor (a_ap/R_in)^2.
    # 'free' : original unbounded free-field pulse-echo (kept for comparison).
    geometry: str = "tube"
    tube_outer_d: float = 300.0e-3   # outer diameter [m]
    tube_len: float = 3.0e-3         # (unused placeholder; range_m is the length)
    gamma_end: float = 0.99          # far-end (thick steel) pressure reflection
    wg_loss_db_m: float = 0.5        # waveguide wall-bounce loss [dB/m] (one-way)

    # --- Transducer aperture (single PMUT is far too small at 3 m; a practical
    #     system uses an array. a_ap is the effective radiating aperture radius) -
    a_ap: float = 5.0e-3       # effective aperture radius [m] (10 mm aperture)

    # --- Transmit drive ------------------------------------------------------
    # We anchor transmit on the ACOUSTIC POWER actually launched into the steel
    # wall, NOT on an imposed surface displacement.  Reason: steel's impedance
    # (~46 MRayl) is enormous, so a "free" PMUT displacement would imply absurd
    # face pressures/powers (P = Z*u).  In reality the heavy steel load clamps
    # the plate and the delivered power is what matters -- and getting power
    # THROUGH steel is the hard part.  W_tx is the documented, tunable knob.
    W_tx: float = 20.0e-3      # acoustic power launched into the wall [W]
    # (For reference the script reports the equivalent uniform face velocity /
    #  displacement so you can sanity-check against a real device.)

    # --- Paper transducer parameters (mode-shape factors ~invariant, Sec 3-4) -
    Q: float = 55.0            # quality factor (paper's design-map default)
    Lambda_p: float = 0.99     # Lambda'  (moment participation, plate regime)
    Lambda1: float = 0.33      # uniform-pressure participation
    Lambda2: float = 0.20      # inertial participation
    e31f: float = -24.24       # eff. transverse piezo coeff [C/m^2] (p-PMUT)
    zp: float = 5.4e-6         # mid-plane of PZT to neutral plane [m] (~stack)
    rho_h: float = 0.0265      # mass per unit area [kg/m^2] (p-PMUT stack, derived)
    Cf: float = 10.0e-12       # charge-amplifier feedback capacitance [F]

    # --- RECIPROCITY-CONSISTENT transducer anchor ----------------------------
    # Transmit and receive both come from ONE coupling so the two-way transfer is
    # passive (no gain > 1) and area-scaling is self-consistent.  We anchor the
    # transmit to the paper's MEASURED deflection sensitivity (Fig 9: a 1 mm,
    # 126 kHz p-PMUT deflects ~3 nm/V at resonance) and DERIVE receive from it via
    # the lumped-model reciprocity relation  (qdot/V) = w00 * (Q/Pin).
    # The aperture is treated as an array of resonant elements (each sized so its
    # first mode = the operating frequency); element count N cancels in SNR.
    defl_sens_ref: float = 3.0e-9   # resonant centre deflection per volt [m/V] @ ref
    a_ref: float = 500.0e-6         # reference element radius [m] (1 mm device)
    f_ref: float = 126.0e3          # reference resonance [Hz]
    De_flex: float = 1.0e-5         # equiv. flexural rigidity [N m] (calib. from ref)
    alpha00: float = 3.19           # 1st clamped-plate eigenvalue (plate regime)
    elec_frac: float = 0.8          # electrode radius / element radius

    # --- Receiver noise model (replaces the arbitrary noise-floor sweep) ------
    # A piezo PMUT into a charge amplifier is dominated by amplifier voltage
    # noise multiplied by the "noise gain" (Cin/Cf), where Cin is set by the
    # PMUT clamped capacitance Cs.  Thin-film PZT has high Cs/area, so this term
    # is large.  Bandwidth is the resonant bandwidth f0/Q (matched receiver).
    temp_K: float = 300.0        # temperature [K]
    amp_en: float = 4.0e-9       # amplifier input voltage noise [V/sqrt(Hz)]
    amp_in: float = 1.0e-15      # amplifier input current noise [A/sqrt(Hz)] (FET)
    eps_r_pzt: float = 1000.0    # PZT relative permittivity
    h_pzt: float = 650.0e-9      # PZT layer thickness [m] (p-PMUT stack)
    fill_factor: float = 0.7     # electrode area / aperture area
    C_stray: float = 5.0e-12     # stray + amplifier input capacitance [F]
    tan_delta: float = 0.02      # PZT dielectric loss tangent (thermal noise)

    # --- Steel wall ----------------------------------------------------------
    d_steel: float = 10.0e-3   # wall thickness [m]
    rho_steel: float = 7850.0  # [kg/m^3]
    c_steel: float = 5900.0    # longitudinal speed [m/s]
    att_steel_db_m: float = 5.0   # steel attenuation at f0 [dB/m] (small)

    # --- Synthetic mineral oil ----------------------------------------------
    rho_oil: float = 850.0     # [kg/m^3]
    c_oil: float = 1440.0      # speed of sound [m/s]
    # Oil attenuation:  alpha(f) = att_oil_db_m * (f/att_ref_f)^att_n   [dB/m]
    # Referenced to a FIXED frequency (not the operating f0) so the frequency
    # sweep is physical.  Mineral/transformer oils vary widely (~1-10 dB/m @ 1 MHz).
    att_oil_db_m: float = 20.0  # [dB/m] at att_ref_f (measured AT 250 kHz)
    att_n: float = 2.0         # frequency exponent
    att_ref_f: float = 250.0e3 # reference frequency for attenuation specs [Hz]

    # --- Transducer coupling to the steel wall ------------------------------
    # 'matched' : transducer well bonded/couplant-matched to steel (Z1 = Z_steel)
    #             -> single steel/oil interface loss, 10 mm adds only delay+atten.
    # 'air'     : dry / air-backed transducer (Z1 = Z_air) -> strong thickness
    #             resonance windows at n*c_steel/(2 d).
    coupling: str = "matched"
    rho_air: float = 1.21
    c_air: float = 343.0

    # --- Target --------------------------------------------------------------
    # 'flat' : large specular reflector (far vessel wall / liquid interface) at
    #          range z.  Mirror -> image source at 2z, round-trip pressure
    #          factor g(2z); the beam is only folded, not re-spread.
    # 'point': small discrete target of radius target_b.  Intercepts the beam
    #          and backscatters as a secondary source -> g(z) out, then a
    #          separate spherical return.  Far more loss; conservative.
    target_type: str = "flat"
    gamma_target: float = 0.95   # pressure reflection coefficient (flat target)
    target_b: float = 5.0e-3     # effective target radius [m] (point target)

    # --- SNR threshold -------------------------------------------------------
    snr_threshold_db: float = 12.0


# ----------------------------------------------------------------------------
# DERIVED ACOUSTIC IMPEDANCES
# ----------------------------------------------------------------------------
def impedances(cfg: Config):
    Z_steel = cfg.rho_steel * cfg.c_steel
    Z_oil = cfg.rho_oil * cfg.c_oil
    Z_air = cfg.rho_air * cfg.c_air
    return Z_steel, Z_oil, Z_air


def tube_inner_radius(cfg: Config):
    """Inner radius of the oil column = (outer_d - 2*wall)/2."""
    return (cfg.tube_outer_d - 2 * cfg.d_steel) / 2.0


def db(x):
    """Power/intensity ratio -> dB (10 log10).  Use for intensity terms."""
    return 10.0 * np.log10(np.clip(x, 1e-300, None))


# ----------------------------------------------------------------------------
# CHANNEL: STEEL WALL  (single-layer acoustic transfer matrix, intensity tau)
# ----------------------------------------------------------------------------
def wall_transmission(f, cfg: Config):
    """
    Intensity transmission coefficient through the steel wall, into the oil.

    Single layer (impedance Z, thickness d) between source medium Z1 and load Z2:

        tau = 4 Z1 Z2 / [ (Z1+Z2)^2 cos^2(kd) + (Z + Z1 Z2 / Z)^2 sin^2(kd) ]

    Limits: kd->0 gives the bare interface 4 Z1 Z2/(Z1+Z2)^2; half-wave (kd=n*pi)
    is transparent; quarter-wave matches if Z=sqrt(Z1 Z2).  Steel attenuation is
    folded in by letting kd carry a small imaginary part.
    """
    Z_steel, Z_oil, Z_air = impedances(cfg)
    f = np.asarray(f, dtype=float)

    Z = Z_steel
    Z2 = Z_oil
    Z1 = Z_steel if cfg.coupling == "matched" else Z_air

    k = 2 * np.pi * f / cfg.c_steel
    # attenuation (Np/m) for the steel layer at this frequency
    a_np = (cfg.att_steel_db_m * (f / cfg.att_ref_f) ** 2) / 8.686
    kd = (k - 1j * a_np) * cfg.d_steel

    cos2 = np.abs(np.cos(kd)) ** 2
    sin2 = np.abs(np.sin(kd)) ** 2
    denom = (Z1 + Z2) ** 2 * cos2 + (Z + Z1 * Z2 / Z) ** 2 * sin2
    tau = 4.0 * Z1 * Z2 / denom
    return tau


# ----------------------------------------------------------------------------
# CHANNEL: OIL PROPAGATION  (spreading + attenuation)
# ----------------------------------------------------------------------------
def oil_alpha_np(f, cfg: Config):
    """Oil pressure attenuation [Np/m] at frequency f (referenced to att_ref_f)."""
    return (cfg.att_oil_db_m * (np.asarray(f, float) / cfg.att_ref_f) ** cfg.att_n) / 8.686


def spreading_pressure(f, z, cfg: Config):
    """
    On-axis baffled-piston far-field pressure spreading factor g(z) = pi a^2/(lambda z)
    (one-way).  Valid in the far field z >> a^2/lambda; we assert that below.
    """
    lam = cfg.c_oil / np.asarray(f, float)
    return np.pi * cfg.a_ap ** 2 / (lam * z)


# ----------------------------------------------------------------------------
# TRANSDUCER  (paper-derived: transmit pressure, receive charge sensitivity)
# ----------------------------------------------------------------------------
def transmit_face_pressure(cfg: Config):
    """
    Acoustic pressure amplitude launched by the transducer into the steel,
    anchored on launched acoustic power W_tx (see Config).

        I_face = W_tx / aperture_area
        P_face = sqrt(2 * Z_steel * I_face)

    Also returns the equivalent uniform face velocity u0 = P_face / Z_steel and
    displacement w0 = u0 / omega for a physical-plausibility sanity check.
    """
    Z_steel, _, _ = impedances(cfg)
    area = np.pi * cfg.a_ap ** 2
    I_face = cfg.W_tx / area
    P_face = np.sqrt(2 * Z_steel * I_face)
    u0 = P_face / Z_steel
    w0 = u0 / (2 * np.pi * cfg.f0)
    return P_face, w0, u0


def receive_charge_sensitivity(cfg: Config):
    """
    Charge generated per unit incident pressure |Qrx/Pin|  [C/Pa], paper Eq (17):

        |Qrx/Pin| = (2 pi Q / omega00^2) (e31f zp / rho_h) (Lambda' Lambda1 / Lambda2)

    omega00 = 2 pi f0 (operate at first resonance).
    """
    w00 = 2 * np.pi * cfg.f0
    s = (2 * np.pi * cfg.Q / w00 ** 2) \
        * (abs(cfg.e31f) * cfg.zp / cfg.rho_h) \
        * (cfg.Lambda_p * cfg.Lambda1 / cfg.Lambda2)
    return s


def transducer(cfg: Config, f=None):
    """
    Reciprocity-consistent transducer model for the aperture at frequency f.

    The aperture (radius a_ap) is an array of N resonant PMUT elements, each
    sized so its first mode equals f.  Transmit and receive both follow from the
    measured deflection sensitivity + the lumped-model reciprocity relation, so
    the two-way transfer is passive and area-scaling cancels in SNR.

    Returns dict with: a_el, N, w0_per_V (resonant), qdot_per_V (total volume
    velocity per volt), s_q (total charge per pascal), Cs (total capacitance).
    """
    if f is None:
        f = cfg.f0
    eps0 = 8.854e-12
    w00 = 2 * np.pi * f
    # element radius so the clamped-plate fundamental = f  (plate regime)
    #   w00 = alpha00^2 * sqrt(De/(rho_h a^4))  ->  a^2 = alpha00^2 sqrt(De/rho_h)/w00
    a_el = np.sqrt(cfg.alpha00 ** 2 * np.sqrt(cfg.De_flex / cfg.rho_h) / w00)
    A_el = np.pi * a_el ** 2
    # resonant deflection per volt scales as a^2 in the plate regime (fixed stack)
    w0_per_V = cfg.defl_sens_ref * (a_el / cfg.a_ref) ** 2
    N = cfg.fill_factor * (cfg.a_ap / a_el) ** 2          # elements tiling aperture
    # transmit: per-element volume velocity per volt at resonance
    qdot_el = w00 * w0_per_V * cfg.Lambda1 * A_el          # [m^3/s/V]
    qdot_per_V = N * qdot_el
    # receive (reciprocity): Q/Pin = (qdot/V)/w00
    s_q = qdot_per_V / w00                                 # [C/Pa] total
    # clamped capacitance (total electrode area)
    Cs_el = eps0 * cfg.eps_r_pzt * np.pi * (cfg.elec_frac * a_el) ** 2 / cfg.h_pzt
    Cs = N * Cs_el
    return {"a_el": a_el, "N": N, "w0_per_V": w0_per_V,
            "qdot_per_V": qdot_per_V, "s_q": s_q, "Cs": Cs}


def receiver_noise(cfg: Config, Cs, f=None):
    """
    Output-referred receiver noise voltage [V rms] of a PMUT + charge amplifier.

    Dominant term: amplifier voltage noise e_n times the capacitive noise gain
    (Cin/Cf), where Cin = Cs + C_stray + Cf and Cs is the (large) PZT clamped
    capacitance.  Plus amplifier current noise through the feedback transimpedance
    and the PZT dielectric-loss thermal (Johnson) noise.  Integrated over the
    resonant equivalent-noise bandwidth ENBW = (pi/2)*(f0/Q).

    Returns (V_noise_out, breakdown_dict, Cs).
    """
    if f is None:
        f = cfg.f0
    kB = 1.380649e-23
    Cin = Cs + cfg.C_stray + cfg.Cf
    noise_gain = Cin / cfg.Cf                      # output/input voltage-noise gain
    enbw = (np.pi / 2) * (f / cfg.Q)               # 2nd-order resonant ENBW [Hz]

    # spectral densities referred to the amplifier output [V/sqrt(Hz)]
    e_v = cfg.amp_en * noise_gain                                  # voltage noise
    e_i = cfg.amp_in / (2 * np.pi * f * cfg.Cf)                    # current noise (transimped.)
    # dielectric-loss (Johnson) noise: parallel loss conductance G = w*Cs*tanδ
    G = 2 * np.pi * f * Cs * cfg.tan_delta
    e_d = np.sqrt(4 * kB * cfg.temp_K * G) / (2 * np.pi * f * cfg.Cf)

    V_v = e_v * np.sqrt(enbw)
    V_i = e_i * np.sqrt(enbw)
    V_d = e_d * np.sqrt(enbw)
    V_tot = np.sqrt(V_v ** 2 + V_i ** 2 + V_d ** 2)
    return V_tot, {"voltage": V_v, "current": V_i, "dielectric": V_d}, Cs


# ----------------------------------------------------------------------------
# FULL ROUND-TRIP LINK BUDGET  (intensity, monostatic flat target)
# ----------------------------------------------------------------------------
def link_budget(cfg: Config, f=None, z=None):
    """
    Returns a dict with the received echo voltage and a per-stage dB breakdown.
    Intensity is tracked through the channel; converted to incident pressure at
    the receiver, then to charge-amp voltage via the paper receive sensitivity.
    """
    if f is None:
        f = cfg.f0
    if z is None:
        z = cfg.range_m

    Z_steel, Z_oil, _ = impedances(cfg)

    # --- far-field sanity ---------------------------------------------------
    lam_oil = cfg.c_oil / f
    rayleigh = cfg.a_ap ** 2 / lam_oil           # near/far boundary
    far_field = z > rayleigh

    # --- transducer (reciprocity-consistent transmit & receive) -------------
    td = transducer(cfg, f)
    u0 = td["qdot_per_V"] * cfg.Vin / (np.pi * cfg.a_ap ** 2)   # face velocity [m/s]
    P_face = Z_steel * u0                        # face pressure into steel [Pa]
    w0 = u0 / (2 * np.pi * f)                     # equiv. face displacement [m]
    I_face = P_face ** 2 / (2 * Z_steel)         # intensity into steel [W/m^2]
    W_tx_eff = I_face * np.pi * cfg.a_ap ** 2     # implied transmit acoustic power [W]

    # --- channel intensity factors -----------------------------------------
    tau = float(wall_transmission(f, cfg))       # per crossing (intensity)
    a_oil = float(oil_alpha_np(f, cfg))
    atten_round_I = np.exp(-2 * a_oil * 2 * z)   # intensity, round trip (2z)

    # Build the round-trip loss as an ORDERED list of (name, intensity factor).
    # I_back is the product of exactly these -> the dB breakdown can never drift
    # out of sync with the computed result (single source of truth).
    factors = [("near-cap out (steel->oil)", tau)]

    if cfg.geometry == "tube":
        # WAVEGUIDE: the oil column confines the beam.  Total transmitted power
        # propagates with attenuation only (no spherical spreading); the loss
        # appears when the small receiver samples the returned, cross-section-
        # filled field -> single fill factor (a_ap/R_in)^2 (=1 if the aperture
        # fills the tube).  Far end is a thick-steel near-total reflector.
        R_in = tube_inner_radius(cfg)
        fill = min(1.0, (cfg.a_ap / R_in) ** 2)
        wg_atten_I = np.exp(-2 * (cfg.wg_loss_db_m / 8.686) * (2 * z))  # wall bounces
        factors.append(("tube fill (a_ap/R_in)^2", fill))
        factors.append(("far-end reflection", cfg.gamma_end ** 2))
        factors.append(("waveguide wall loss", wg_atten_I))
    elif cfg.target_type == "flat":
        # Specular mirror: image source at 2z. Round-trip geometric pressure
        # factor is g(2z); apply its square once as an intensity factor.
        g2z = float(spreading_pressure(f, 2 * z, cfg))
        factors.append(("geom spreading (flat, image @2z)", g2z ** 2))
        factors.append(("target reflection", cfg.gamma_target ** 2))
    else:
        # Point/small target: spread out to z, backscatter, spread back to z.
        # Backscatter intensity cross-section of a rigid target ~ pi b^2;
        # fraction returned to the aperture ~ area / (4 pi z^2).
        g = float(spreading_pressure(f, z, cfg))
        sigma = np.pi * cfg.target_b ** 2                 # backscatter x-section
        recv = (np.pi * cfg.a_ap ** 2) / (4 * np.pi * z ** 2)
        factors.append(("spreading out", g ** 2))
        factors.append(("target backscatter", sigma * recv * cfg.gamma_target ** 2))
        factors.append(("spreading return", g ** 2))

    factors.append(("near-cap return (oil->steel)", tau))
    factors.append(("oil absorption (round trip)", atten_round_I))

    I_back = I_face
    for _, fac in factors:
        I_back *= fac

    # --- receive (same transducer -> reciprocity-consistent) ----------------
    P_inc = np.sqrt(2 * Z_steel * I_back)        # incident pressure in steel [Pa]
    s_q = td["s_q"]                              # [C/Pa] (derived from transmit)
    Q_echo = s_q * P_inc                         # received charge [C]
    V_echo = Q_echo / cfg.Cf                     # charge-amp output [V]

    # --- modeled receiver noise + absolute SNR ------------------------------
    Cs = td["Cs"]
    V_noise, n_parts, _ = receiver_noise(cfg, Cs, f)
    snr_model = 20.0 * np.log10(max(V_echo, 1e-300) / V_noise)
    # passivity flag: a passive 2-way link cannot give voltage gain > 1
    reciprocity_violation = V_echo > cfg.Vin

    # --- per-stage dB breakdown (intensity) ---------------------------------
    breakdown = {"transmit source (ref 0 dB)": 0.0}
    for name, fac in factors:
        breakdown[name] = db(fac)

    return {
        "f": f, "z": z,
        "P_face": P_face, "w0": w0, "u0": u0, "I_face": I_face,
        "tau_wall": tau, "alpha_oil_np": a_oil,
        "rayleigh": rayleigh, "far_field": far_field,
        "I_back": I_back, "P_inc": P_inc,
        "Q_echo": Q_echo, "V_echo": V_echo,
        "charge_sens": s_q,
        "V_noise": V_noise, "noise_parts": n_parts, "Cs": Cs,
        "snr_model": snr_model, "reciprocity_violation": reciprocity_violation,
        "W_tx_eff": W_tx_eff, "a_el": td["a_el"], "N": td["N"],
        "w0_per_V": td["w0_per_V"],
        "tof_round": 2 * (z / cfg.c_oil + cfg.d_steel / cfg.c_steel),
        "breakdown": breakdown,
    }


def snr_db(V_echo, V_noise):
    """SNR in dB from echo and noise voltage amplitudes (20 log10)."""
    return 20.0 * np.log10(np.clip(V_echo / V_noise, 1e-300, None))


# ----------------------------------------------------------------------------
# VISUALISATION
# ----------------------------------------------------------------------------
def make_compare_figure(freqs_khz=(100.0, 135.0, 200.0),
                        outfile="pmut_tube_results.png"):
    """Single 6-panel figure comparing the calibrated model at several
    operating frequencies (defaults match the wave animations)."""
    base = Config()
    colors = plt.cm.viridis(np.linspace(0.12, 0.78, len(freqs_khz)))
    lbs = [link_budget(dataclasses.replace(base, f0=fk * 1e3), f=fk * 1e3)
           for fk in freqs_khz]

    fig = plt.figure(figsize=(16, 9.5), constrained_layout=True)
    fig.suptitle(
        f"PMUT pulse-echo in oil-filled steel tube  "
        f"(L={base.range_m:.0f} m, OD={base.tube_outer_d*1e3:.0f} mm, "
        f"wall={base.d_steel*1e3:.0f} mm)  —  comparison at "
        f"{', '.join(f'{int(f)}' for f in freqs_khz)} kHz\n"
        f"calibrated transducer (Dangi & Pratap 2017) + waveguide channel "
        f"+ modeled receiver noise",
        fontsize=12, fontweight="bold")
    gs = fig.add_gridspec(2, 3)

    # ---- Panel 1: steel-wall transmission (matched + air-backed) ------------
    ax1 = fig.add_subplot(gs[0, 0])
    fsw = np.linspace(0.05e6, 1.2e6, 2000)
    ax1.plot(fsw / 1e3, db(wall_transmission(fsw, dataclasses.replace(base, coupling="matched"))),
             color="C0", lw=1.6, label="matched coupling")
    ax1.plot(fsw / 1e3, db(wall_transmission(fsw, dataclasses.replace(base, coupling="air"))),
             color="0.55", lw=1.0, alpha=0.85, label="air-backed (dry)")
    for n in range(1, 5):
        fr = n * base.c_steel / (2 * base.d_steel)
        if fr < 1.2e6:
            ax1.axvline(fr / 1e3, color="grey", ls=":", lw=0.7)
    for fk, c in zip(freqs_khz, colors):
        ax1.axvline(fk, color=c, lw=1.4)
    ax1.set_xlabel("Frequency [kHz]"); ax1.set_ylabel("Wall transmission [dB]")
    ax1.set_title("Steel-wall transmission\n(dotted = n·c/2d resonances; lines = op. freqs)")
    ax1.legend(fontsize=7, loc="lower right"); ax1.grid(alpha=0.3); ax1.set_ylim(-65, 3)

    # ---- Panel 2: round-trip link-budget breakdown per frequency ------------
    ax2 = fig.add_subplot(gs[0, 1])
    for lb, fk, c in zip(lbs, freqs_khz, colors):
        vals = list(lb["breakdown"].values())
        cum = np.cumsum(vals)
        ax2.step(range(len(cum)), cum, where="mid", color=c, lw=1.8,
                 label=f"{int(fk)} kHz: {cum[-1]:.0f} dB")
    stages = list(lbs[0]["breakdown"].keys())
    ax2.set_xticks(range(len(stages)))
    ax2.set_xticklabels(stages, rotation=40, ha="right", fontsize=6)
    ax2.set_ylabel("Cumulative intensity [dB]")
    ax2.set_title("Round-trip link budget\n(only oil absorption differs)")
    ax2.legend(fontsize=7); ax2.grid(alpha=0.3)

    # ---- Panel 3 (headline): absolute SNR vs frequency ----------------------
    ax3 = fig.add_subplot(gs[0, 2])
    fs2 = np.linspace(0.05e6, 0.4e6, 400)        # focus where it matters
    snr_f = np.array([link_budget(dataclasses.replace(base, f0=ff), f=ff)["snr_model"]
                      for ff in fs2])
    ax3.plot(fs2 / 1e3, snr_f, color="C4", lw=2)
    ax3.axhline(base.snr_threshold_db, color="green", ls="--", lw=1.3, label="12 dB")
    ax3.axhline(0, color="0.6", lw=0.7)
    for lb, fk, c in zip(lbs, freqs_khz, colors):
        ax3.plot(fk, lb["snr_model"], "o", color=c, ms=9, zorder=5)
        ax3.annotate(f"{int(fk)} kHz\n{lb['snr_model']:+.0f} dB",
                     (fk, lb["snr_model"]), fontsize=7, color=c,
                     xytext=(6, -4), textcoords="offset points")
    ax3.set_xlabel("Frequency [kHz]"); ax3.set_ylabel("absolute SNR [dB]")
    ax3.set_ylim(-80, 60)                        # clip the exp(-kf²) cliff tail
    ax3.set_title("SNR vs frequency (modeled noise)\ncurve plunges off-scale >250 kHz (cliff)")
    ax3.legend(fontsize=7); ax3.grid(alpha=0.3)

    # ---- Panel 4: received echo WAVEFORM at TOF (3 stacked, own scales) ------
    # Real receiver voltage (echo burst + modeled noise) around the round-trip
    # time.  Echo amplitudes span ~2000:1, so each frequency gets its own y-scale
    # (absolute, not normalised) -- the only way to keep the actual waveform
    # readable for all three.
    tof = lbs[0]["tof_round"]
    gs4 = gs[1, 0].subgridspec(3, 1, hspace=0.32)
    rng = np.random.default_rng(0)
    fs_t = 20e6
    for i, (lb, fk, c) in enumerate(zip(lbs, freqs_khz, colors)):
        ax = fig.add_subplot(gs4[i])
        f = fk * 1e3
        tb = 5 / f
        echoV = lb["V_echo"]; noiseV = lb["V_noise"]
        tw = np.arange(tof - 12e-6, tof + tb + 12e-6, 1 / fs_t)
        te = tw - tof
        sig = (np.where((te >= 0) & (te < tb),
               np.sin(2*np.pi*f*te) * np.sin(np.pi*np.clip(te/tb, 0, 1))**2, 0.0) * echoV
               + rng.normal(0, noiseV, size=tw.size))
        u, lab = (1e3, "mV") if max(echoV, noiseV) > 1e-3 else (1e6, "µV")
        ax.plot((tw - tof) * 1e6, sig * u, color=c, lw=0.6)
        ax.axvline(0, color="red", ls="--", lw=0.9)
        ax.text(0.015, 0.80, f"{int(fk)} kHz  ·  echo {echoV*u:.3g} {lab}  ·  "
                f"SNR {lb['snr_model']:+.0f} dB", transform=ax.transAxes, fontsize=7.5,
                color=("green" if lb["snr_model"] >= 12 else "firebrick"))
        ax.set_ylabel(lab, fontsize=7); ax.tick_params(labelsize=6); ax.grid(alpha=0.3)
        if i == 0:
            ax.set_title("Received echo waveform at round-trip time "
                         "(own scale each; echo + modeled noise)", fontsize=9)
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("time relative to echo @ 4.17 ms  [µs]", fontsize=8)

    # ---- Panel 5: echo voltage & modeled noise floor vs frequency -----------
    ax5 = fig.add_subplot(gs[1, 1])
    ev, nv = [], []
    for ff in fs2:
        l = link_budget(dataclasses.replace(base, f0=ff), f=ff)
        ev.append(l["V_echo"] * 1e6); nv.append(l["V_noise"] * 1e6)
    ax5.semilogy(fs2 / 1e3, ev, color="C0", lw=2, label="echo |V|")
    ax5.semilogy(fs2 / 1e3, nv, color="C3", lw=2, label="noise floor")
    for lb, fk, c in zip(lbs, freqs_khz, colors):
        ax5.plot(fk, lb["V_echo"] * 1e6, "o", color=c, ms=8, zorder=5)
    ax5.set_xlabel("Frequency [kHz]"); ax5.set_ylabel("Voltage [µV]")
    ax5.set_ylim(1e0, 1e7)                       # keep crossover region readable
    ax5.set_title("Echo vs modeled noise\n(Cs %.0f nF, en=%.0f nV/√Hz)"
                  % (lbs[0]["Cs"] * 1e9, base.amp_en * 1e9))
    ax5.legend(fontsize=7); ax5.grid(alpha=0.3, which="both")

    # ---- Panel 6: SNR vs (range, frequency) feasibility map -----------------
    ax6 = fig.add_subplot(gs[1, 2])
    rr = np.linspace(0.2, base.range_m, 70)
    ff6 = np.linspace(50e3, 300e3, 70)
    Rr, Ff = np.meshgrid(rr, ff6)
    Z = np.zeros_like(Rr)
    for i in range(Ff.shape[0]):
        for j in range(Ff.shape[1]):
            l = link_budget(dataclasses.replace(base, f0=Ff[i, j], range_m=Rr[i, j]),
                            f=Ff[i, j])
            Z[i, j] = l["snr_model"]
    pcm = ax6.pcolormesh(Rr, Ff / 1e3, Z, shading="auto", cmap="viridis",
                         vmin=-40, vmax=60)
    cs = ax6.contour(Rr, Ff / 1e3, Z, levels=[base.snr_threshold_db],
                     colors="white", linewidths=2)
    ax6.clabel(cs, fmt="12 dB")
    for fk, c in zip(freqs_khz, colors):
        ax6.plot(base.range_m, fk, "o", color="red", ms=8, mec="white")
    fig.colorbar(pcm, ax=ax6, label="SNR [dB]")
    ax6.set_xlabel("Target range [m]"); ax6.set_ylabel("Frequency [kHz]")
    ax6.set_title("Feasibility: SNR vs range & frequency\n(dots = 3 m operating points)")

    fig.savefig(outfile, dpi=130, bbox_inches="tight")
    print(f"\nSaved figure -> {outfile}")


def make_figures(cfg: Config, outfile="pmut_echo_results.png"):
    lb = link_budget(cfg)

    R_in = tube_inner_radius(cfg)
    fig = plt.figure(figsize=(15, 10))
    fig.suptitle(
        f"PMUT pulse-echo in oil-filled steel tube  "
        f"(L={cfg.range_m:.0f} m, OD={cfg.tube_outer_d*1e3:.0f} mm, "
        f"wall={cfg.d_steel*1e3:.0f} mm)  @ {cfg.f0/1e3:.0f} kHz\n"
        f"(transducer model from Dangi & Pratap 2017; waveguide channel added here)",
        fontsize=13, fontweight="bold")

    # ---- Panel 1: steel wall transmission spectrum -------------------------
    ax1 = fig.add_subplot(2, 3, 1)
    fsweep = np.linspace(0.05e6, 2.0e6, 4000)
    tau_matched = wall_transmission(fsweep, dataclasses.replace(cfg, coupling="matched"))
    tau_air = wall_transmission(fsweep, dataclasses.replace(cfg, coupling="air"))
    ax1.plot(fsweep / 1e6, db(tau_matched), label="matched coupling", lw=1.8)
    ax1.plot(fsweep / 1e6, db(tau_air), label="air-backed (dry)", lw=1.0, alpha=0.8)
    for n in range(1, 7):
        fr = n * cfg.c_steel / (2 * cfg.d_steel)
        if fr < 2.0e6:
            ax1.axvline(fr / 1e6, color="grey", ls=":", lw=0.8)
    ax1.axvline(cfg.f0 / 1e6, color="red", lw=1.5, label="operating f")
    ax1.set_xlabel("Frequency [MHz]")
    ax1.set_ylabel("Wall intensity transmission [dB]")
    ax1.set_title("Steel-wall transmission\n(dotted = n·c/2d resonances)")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    ax1.set_ylim(-60, 2)

    # ---- Panel 2: link-budget waterfall ------------------------------------
    ax2 = fig.add_subplot(2, 3, 2)
    stages = list(lb["breakdown"].keys())
    vals = list(lb["breakdown"].values())
    cum = np.cumsum(vals)
    ax2.step(range(len(cum)), cum, where="mid", color="C3", lw=2)
    ax2.fill_between(range(len(cum)), cum, cum.min() - 10, step="mid", alpha=0.15, color="C3")
    ax2.set_xticks(range(len(stages)))
    ax2.set_xticklabels(stages, rotation=40, ha="right", fontsize=7)
    ax2.set_ylabel("Cumulative intensity [dB re source]")
    ax2.set_title(f"Round-trip link budget\n(total {cum[-1]:.0f} dB acoustic)")
    ax2.grid(alpha=0.3)

    # ---- Panel 3: echo vs MODELED noise floor, vs frequency ----------------
    # Replaces the arbitrary noise-floor sweep with a physical receiver-noise
    # model.  Shows the echo voltage and the modeled noise floor (both in µV)
    # so the crossing (SNR=0) and 12 dB margin are read directly.
    ax3 = fig.add_subplot(2, 3, 3)
    fmod = np.linspace(0.05e6, 1.2e6, 400)
    e_v, n_v = [], []
    for ff in fmod:
        lbf = link_budget(dataclasses.replace(cfg, f0=ff), f=ff)
        e_v.append(lbf["V_echo"] * 1e6)
        n_v.append(lbf["V_noise"] * 1e6)
    e_v, n_v = np.array(e_v), np.array(n_v)
    ax3.semilogy(fmod / 1e3, e_v, color="C0", lw=2, label="echo |V|")
    ax3.semilogy(fmod / 1e3, n_v, color="C3", lw=2, label="modeled noise floor")
    ax3.semilogy(fmod / 1e3, n_v * 10 ** (cfg.snr_threshold_db / 20), color="C2",
                 ls="--", lw=1, label="noise + 12 dB")
    ax3.axvline(cfg.f0 / 1e3, color="k", lw=1, ls=":", label=f"{cfg.f0/1e3:.0f} kHz")
    ax3.set_xlabel("Frequency [kHz]")
    ax3.set_ylabel("Voltage [µV]")
    ax3.set_title(f"Echo vs modeled noise\n(Cs={lb['Cs']*1e9:.0f} nF, en={cfg.amp_en*1e9:.0f} nV/√Hz)")
    ax3.legend(fontsize=7); ax3.grid(alpha=0.3, which="both")
    ax3.set_ylim(max(1e-3, n_v.min() / 10), max(e_v.max(), n_v.max()) * 10)

    # ---- Panel 4: received echo (zoom at time-of-flight) -------------------
    # The echo is a ~5 us burst arriving ~4 ms after TX, so we zoom on a window
    # around the round-trip time and show the actual echo voltage + noise set so
    # the displayed SNR equals the threshold (illustrates the detection limit).
    ax4 = fig.add_subplot(2, 3, 4)
    fs = 50e6
    tof = lb["tof_round"]
    n_cycles = 5
    tb_len = n_cycles / cfg.f0
    win = 12e-6
    tw = np.arange(tof - win, tof + tb_len + win, 1 / fs)
    te = tw - tof
    echo = np.where((te >= 0) & (te < tb_len),
                    np.sin(2 * np.pi * cfg.f0 * te) *
                    np.sin(np.pi * np.clip(te / tb_len, 0, 1)) ** 2, 0.0) \
        * lb["V_echo"] * 1e6                       # echo in µV
    # noise at exactly the threshold SNR so the plot shows the detection margin
    v_noise = lb["V_echo"] / (10 ** (cfg.snr_threshold_db / 20)) * 1e6  # µV rms
    rng = np.random.default_rng(0)
    noise = rng.normal(0, v_noise, size=tw.size)
    ax4.plot((tw - tof) * 1e6, echo + noise, color="C0", lw=0.7,
             label=f"echo + noise @ {cfg.snr_threshold_db:.0f} dB")
    ax4.plot((tw - tof) * 1e6, echo, color="C3", lw=1.4, label="echo (clean)")
    ax4.set_xlabel(f"Time relative to TX  [µs]   (echo @ {tof*1e3:.2f} ms)")
    ax4.set_ylabel("RX voltage [µV]")
    ax4.set_title("Received echo at target range\n(zoom at round-trip time)")
    ax4.legend(fontsize=8)
    ax4.grid(alpha=0.3)
    # inset: full TX->echo timeline (envelopes) for time-of-flight context
    axin = ax4.inset_axes([0.58, 0.62, 0.40, 0.34])
    axin.plot([0, 0], [0, 1], color="C1", lw=2)
    axin.plot([tof * 1e3, tof * 1e3], [0, lb["V_echo"] / cfg.Vin], color="C3", lw=2)
    axin.set_xlim(-0.2, tof * 1e3 * 1.15)
    axin.set_yscale("log")
    axin.set_ylim(1e-5, 2)
    axin.set_xlabel("ms", fontsize=6)
    axin.set_title("TX → echo (norm.)", fontsize=6)
    axin.tick_params(labelsize=5)

    # ---- Panel 5: ABSOLUTE SNR vs frequency (modeled noise) ----------------
    ax5 = fig.add_subplot(2, 3, 5)
    fs2 = np.linspace(0.05e6, 1.2e6, 400)
    snr_f, recip = [], []
    for ff in fs2:
        lbf = link_budget(dataclasses.replace(cfg, f0=ff), f=ff)
        snr_f.append(lbf["snr_model"]); recip.append(lbf["reciprocity_violation"])
    snr_f = np.array(snr_f); recip = np.array(recip)
    ax5.plot(fs2 / 1e3, snr_f, color="C4", lw=2)
    # shade the low-f band where the (decoupled) signal model over-predicts
    if recip.any():
        fbad = fs2[recip].max()
        ax5.axvspan(fs2[0] / 1e3, fbad / 1e3, color="orange", alpha=0.12,
                    label="signal over-predicted\n(reciprocity)")
    ax5.axhline(cfg.snr_threshold_db, color="green", ls="--", lw=1.5, label="12 dB")
    ax5.axhline(0, color="0.6", lw=0.8)
    ax5.axvline(cfg.f0 / 1e3, color="red", lw=1.2, label=f"operating {cfg.f0/1e3:.0f} kHz")
    ax5.axvline(1000, color="grey", lw=1.0, ls=":", label="1 MHz")
    ax5.set_xlabel("Frequency [kHz]")
    ax5.set_ylabel("absolute SNR [dB]")
    ax5.set_title("SNR vs frequency (MODELED noise)\nsignal collapse ~exp(-k f²)")
    ax5.legend(fontsize=7)
    ax5.grid(alpha=0.3)

    # ---- Panel 6: SNR vs range x aperture (feasibility map) ----------------
    ax6 = fig.add_subplot(2, 3, 6)
    ranges = np.linspace(0.2, cfg.range_m, 100)
    apertures = np.linspace(1e-3, R_in, 100)         # up to full tube radius
    Zr, Za = np.meshgrid(ranges, apertures)
    snr_map = np.zeros_like(Zr)
    V_noise_fixed = lb["V_noise"]                    # fixed design receiver
    for i in range(Za.shape[0]):
        for j in range(Za.shape[1]):
            c2 = dataclasses.replace(cfg, a_ap=Za[i, j], range_m=Zr[i, j])
            lbm = link_budget(c2)
            snr_map[i, j] = 20 * np.log10(max(lbm["V_echo"], 1e-300) / V_noise_fixed)
    pcm = ax6.pcolormesh(Zr, Za * 1e3, snr_map, shading="auto", cmap="viridis")
    cs = ax6.contour(Zr, Za * 1e3, snr_map, levels=[cfg.snr_threshold_db],
                     colors="white", linewidths=2)
    ax6.clabel(cs, fmt=f"{cfg.snr_threshold_db:.0f} dB")
    ax6.plot(cfg.range_m, cfg.a_ap * 1e3, "r*", ms=15, label="design point")
    fig.colorbar(pcm, ax=ax6, label="SNR [dB] (modeled noise)")
    ax6.set_xlabel("Target range [m]")
    ax6.set_ylabel("Aperture radius [mm]")
    ax6.set_title("Feasibility map")
    ax6.legend(fontsize=8, loc="upper left")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outfile, dpi=130)
    print(f"\nSaved figure -> {outfile}")
    return lb


# ----------------------------------------------------------------------------
# TEXT REPORT
# ----------------------------------------------------------------------------
def print_report(cfg: Config, lb):
    Z_steel, Z_oil, Z_air = impedances(cfg)
    print("=" * 70)
    print("PMUT PULSE-ECHO LINK BUDGET")
    print("=" * 70)
    print(f"  Centre frequency      : {cfg.f0/1e6:.3f} MHz")
    print(f"  Drive voltage         : {cfg.Vin:.1f} V")
    print(f"  Aperture radius       : {cfg.a_ap*1e3:.2f} mm")
    print(f"  Target range (1-way)  : {cfg.range_m:.2f} m")
    print(f"  Steel wall            : {cfg.d_steel*1e3:.1f} mm  "
          f"(Z={Z_steel/1e6:.1f} MRayl, c={cfg.c_steel} m/s)")
    print(f"  Oil                   : Z={Z_oil/1e6:.2f} MRayl, c={cfg.c_oil} m/s, "
          f"atten={cfg.att_oil_db_m} dB/m @ {cfg.att_ref_f/1e3:.0f} kHz")
    print(f"  Coupling              : {cfg.coupling}")
    if cfg.geometry == "tube":
        R_in = tube_inner_radius(cfg)
        print(f"  Geometry              : TUBE  outer D={cfg.tube_outer_d*1e3:.0f} mm, "
              f"wall={cfg.d_steel*1e3:.0f} mm, inner R={R_in*1e3:.0f} mm")
        print(f"  Tube fill (a/R_in)^2  : {min(1.0,(cfg.a_ap/R_in)**2):.4f} "
              f"({db(min(1.0,(cfg.a_ap/R_in)**2)):.1f} dB)")
    else:
        print(f"  Geometry              : FREE field, target={cfg.target_type}")
    print("-" * 70)
    print(f"  Transducer (reciprocity-consistent, anchored to {cfg.defl_sens_ref*1e9:.0f} nm/V):")
    print(f"      element radius a_el   : {lb['a_el']*1e6:.0f} µm  (resonant at f0)  "
          f"x {lb['N']:.0f} elements in aperture")
    print(f"      defl. sensitivity     : {lb['w0_per_V']*1e9:.2f} nm/V (resonant)")
    print(f"  TX power into wall (derived): {lb['W_tx_eff']*1e3:.1f} mW "
          f"(u0={lb['u0']*1e3:.2f} mm/s, P_face={lb['P_face']/1e3:.1f} kPa)")
    print(f"  Steel resonances at   : "
          + ", ".join(f"{n*cfg.c_steel/(2*cfg.d_steel)/1e3:.0f}" for n in range(1, 6))
          + " kHz")
    print(f"  Wall transmission     : {db(lb['tau_wall']):.1f} dB / crossing")
    print(f"  Far-field?            : {lb['far_field']} "
          f"(Rayleigh dist {lb['rayleigh']*1e3:.2f} mm)")
    print("-" * 70)
    print("  Round-trip breakdown (intensity, dB):")
    for k, v in lb["breakdown"].items():
        print(f"      {k:<32s}: {v:8.1f}")
    print(f"      {'TOTAL acoustic':<32s}: "
          f"{sum(lb['breakdown'].values()):8.1f}")
    print("-" * 70)
    print(f"  Charge sens (reciprocity): {lb['charge_sens']*1e15:.3f} fC/Pa "
          f"  [paper Eq.17 would give {receive_charge_sensitivity(cfg)*1e15:.1f} fC/Pa "
          f"-- not reciprocity-consistent]")
    print(f"  Incident echo pressure: {lb['P_inc']*1e3:.4g} mPa")
    print(f"  Echo charge           : {lb['Q_echo']*1e15:.4g} fC")
    print(f"  Echo voltage          : {lb['V_echo']*1e6:.4g} µV")
    print(f"  Round-trip time-of-flt: {lb['tof_round']*1e6:.2f} µs")
    print("-" * 70)
    # MODELED receiver noise (replaces the arbitrary floor sweep)
    np_ = lb["noise_parts"]
    print("  Receiver noise model (charge amp + PZT capacitance):")
    print(f"      PMUT capacitance Cs   : {lb['Cs']*1e9:.1f} nF  "
          f"(noise gain ~ {(lb['Cs']+cfg.C_stray+cfg.Cf)/cfg.Cf:.0f})")
    print(f"      noise floor (output)  : {lb['V_noise']*1e6:.3g} µV rms"
          f"   [V={np_['voltage']*1e6:.2g}, I={np_['current']*1e6:.2g}, "
          f"diel={np_['dielectric']*1e6:.2g} µV]")
    s = lb["snr_model"]
    flag = "OK" if s >= cfg.snr_threshold_db else "FAIL"
    print(f"      echo {lb['V_echo']*1e6:.3g} µV vs noise {lb['V_noise']*1e6:.3g} µV"
          f"  -> SNR = {s:.1f} dB   [{flag} vs {cfg.snr_threshold_db:.0f} dB]")
    if lb["reciprocity_violation"]:
        print("      !! V_echo > Vin: signal over-predicted here (reciprocity) — SNR is an upper bound")
    # SNR-optimal frequency (modeled noise)
    fscan = np.linspace(50e3, 1.2e6, 400)
    snrf = np.array([link_budget(dataclasses.replace(cfg, f0=ff), f=ff)["snr_model"]
                     for ff in fscan])
    fopt = fscan[int(np.argmax(snrf))]
    s1M = link_budget(dataclasses.replace(cfg, f0=1e6), f=1e6)["snr_model"]
    print(f"  SNR-optimal frequency : {fopt/1e3:.0f} kHz (SNR {snrf.max():.0f} dB)   "
          f"[1 MHz gives {s1M:.0f} dB]")
    print("=" * 70)


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    import sys
    # operating frequencies to compare [kHz]; default matches the animations
    freqs = [float(x) for x in sys.argv[1:]] or [100.0, 135.0, 200.0]
    for fk in freqs:
        print_report(dataclasses.replace(Config(), f0=fk * 1e3),
                     link_budget(dataclasses.replace(Config(), f0=fk * 1e3), f=fk * 1e3))
    make_compare_figure(tuple(freqs), outfile="pmut_tube_results.png")
