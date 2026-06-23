# Decisions log

Key modeling decisions and why, in order. Each one is a fork where the result
would differ materially under another choice.

1. **Anchor transmit on launched acoustic power, not imposed displacement.**
   Imposing a free PMUT velocity into steel (Z=46 MRayl) implies GPa face pressures
   and kilowatts of power. Power is the physical knob; later replaced by the
   reciprocity-consistent transmit (decision 7).

2. **Treat the geometry as a waveguide, not free field.** The oil-filled tube
   confines the beam, so spherical 1/r^2 spreading (-41 dB at 3 m) is replaced by a
   bounded cross-section fill factor (a_ap/R_in)^2. This is the single biggest
   change to feasibility and the reason the tube case differs from free field.

2b. **Far end = thick steel reflector (Gamma ~ 0.99).** Per the human's spec, the
    end-cap is thick enough to ignore the steel-air interface, so it is a near-total
    reflector; no second medium transition on the return path.

3. **Oil attenuation referenced to 250 kHz, scaling as f^2.** alpha(f) =
   20 dB/m * (f/250kHz)^2. Classical thermoviscous absorption. This produces the
   dominant frequency dependence: 11 dB round trip at 75 kHz, 120 dB at 250 kHz,
   1920 dB at 1 MHz. (An early version wrongly referenced 1 MHz; correcting it moved
   the operating point by 112 dB.)

4. **Steel sound speed reduced in the FDTD only.** c_steel set to 2500 m/s (real
   density kept) purely to relax the CFL time step; the impedance contrast and hence
   the reflection are preserved, and the oil speed (1440 m/s) is exact so the
   round-trip time-of-flight (~4.17 ms) is physical. Visualization-only trade-off.

5. **Replace the arbitrary noise floor with a modeled one.** Charge-amplifier
   voltage noise times capacitive noise gain (Cin/Cf) dominates, where Cs is the
   (large) thin-film-PZT clamped capacitance; plus current noise and dielectric-loss
   Johnson noise, integrated over the resonant bandwidth f0/Q. This revealed the link
   is receiver-capacitance limited, not channel limited at low frequency.

6. **Diagnose the 1 MHz SNR cliff as physics, not a bug.** Absorption in dB scales
   as f^2 over a fixed path, so amplitude collapses as exp(-k f^2). No noise model
   changes it; it is signal-side.

7. **Tie transmit and receive to one coupling (reciprocity).** The decoupled model
   (power-in transmit, Eq.-17 receive) implied gain > 1 at low frequency (a 140 V
   echo from a 10 V drive). Using the lumped reciprocity relation qdot/V = w00*(Q/Pin)
   and anchoring to the paper's MEASURED 3 nm/V transmit makes the two-way passive
   (V_echo/Vin <= 0.77 everywhere) and area-self-consistent. The reciprocity-derived
   receive sensitivity (27 fC/Pa at 250 kHz) is ~41x smaller than Eq. 17 (1130 fC/Pa);
   the paper's own transmit and receive figures are mutually inconsistent.

8. **Report absolute SNR with an explicit uncertainty caveat.** The paper's
   sensitivity numbers (Eq. 17, Table 3, Fig. 9) disagree by >100x, so absolute SNR
   carries order-of-magnitude uncertainty. The robust, reported results are the
   relative trends, the f^2 cliff, the reciprocity constraint, the
   capacitance-limited regime, and the detection-boundary structure (~135 kHz).

9. **Operating points for the comparison: 100 / 135 / 200 kHz.** Bracket the
   ~135 kHz detection boundary (detectable / knife-edge / buried).
