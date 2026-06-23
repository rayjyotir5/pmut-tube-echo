# Progress log — PMUT pulse-echo in an oil-filled steel pipe

## Status: complete (paper compiled)

## What this is
A system-level study of whether a piezoelectric micromachined ultrasonic
transducer (PMUT) can range a target by pulse-echo through an oil-filled steel
pipe. Starts from the analytical transducer model of Dangi & Pratap (2017),
adds an acoustic-channel model for the pipe waveguide, a first-principles
receiver-noise model, and a 2D FDTD visualization, and ties transmit/receive
together by reciprocity so the absolute link budget is passive and consistent.

## Headline findings
1. **Confinement governs feasibility.** The oil-filled pipe is a waveguide; the
   dominant free-field spreading loss is replaced by a bounded fill factor. The
   pipe case is far more favorable than free-field radiation.
2. **An f^2 absorption cliff sets the operating frequency.** Oil absorption over
   the 6 m round trip is 11 / 120 / 1920 dB at 75 / 250 / 1000 kHz. The
   conventional "use ~1 MHz for resolution" choice is annihilated; the usable band
   is far lower (tens to ~135 kHz).
3. **The receiver is electronics-limited, not acoustics-limited.** A modeled
   charge-amp + thin-film-PZT capacitance noise floor (~16 mV-equivalent) is what
   sets detectability, not the channel; an arbitrary 10 uV floor was ~70 dB
   optimistic.
4. **Reciprocity matters.** Pulling transmit and receive from independent anchors
   gives an impossible gain > 1; tying them to one coupling (anchored to the
   measured 3 nm/V deflection) makes the model passive and shifts the detection
   boundary to ~135 kHz (at 20 dB/m oil).

## Artifacts
- src/pmut_echo_sim.py   — calibrated 1D link budget + noise model + figure
- src/pmut_wave2d.py     — 2D FDTD wave solver + 3-panel synced animation
- src/make_paper_figures.py — schematic + wave montage
- figures/pmut_tube_results.png — 6-panel comparison at 100/135/200 kHz
- figures/fig_schematic.png, figures/fig_wave_montage.png
- animations (in working dir pmut_sim/): 100khz.mp4, 135khz.mp4, 200khz.mp4
- paper/main.tex, paper/appendix.tex, paper/main.pdf

## Reproduce
```
cd src
python3 pmut_echo_sim.py 100 135 200     # -> ../figures/pmut_tube_results.png + reports
python3 pmut_wave2d.py 100               # -> 100khz.mp4 / 100khz.gif (and 135, 200)
python3 make_paper_figures.py            # -> schematic + montage
cd ../paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Honest limitations
- Absolute SNR has order-of-magnitude uncertainty: the source paper's transmit and
  receive sensitivities are mutually inconsistent (>100x spread). Relative trends,
  the cliff, and the reciprocity/capacitance conclusions are robust.
- The 2D FDTD is for visualization (axial slice, reduced steel speed); it is not the
  source of the SNR numbers.
- No experiment; this is a modeling study. Oil attenuation and PZT permittivity are
  representative, not measured for a specific oil/device.
