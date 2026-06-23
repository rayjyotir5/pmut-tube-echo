# Modelling of a PMUT pulse-echo through an oil-filled steel tube

Code, figures, and wave-propagation videos for the paper *"Modelling of a PMUT
Pulse-Echo Through an Oil-Filled Steel Tube"* (Jyotirmoy Ray, Wipro Research).
See [`paper/main.pdf`](paper/main.pdf).

A piezoelectric micromachined ultrasonic transducer (PMUT) is placed at one end of a
3 m oil-filled steel tube (OD 300 mm, 10 mm wall) and used for pulse-echo ranging off
the thick-steel far end. Starting from the transducer model of Dangi & Pratap (2017),
this repo adds a waveguide acoustic-channel link budget, a reciprocity-consistent
transmit/receive coupling, a first-principles receiver-noise model, and a 2D FDTD
wave simulation.

## Key results
- **The operating frequency is set by an f² oil-absorption cliff.** Round-trip loss is
  ~11 / 120 / 1920 dB at 75 / 250 / 1000 kHz. The conventional ~1 MHz imaging choice is
  annihilated; the usable band is near and below ~135 kHz.
- **Confinement, not free-field radiation, governs feasibility.** The tube guides the
  beam, so spherical spreading is replaced by a bounded aperture fill factor.
- **The receiver is the limit, not the channel.** A thin-film-PZT clamped capacitance
  (~479 nF) into a low-noise charge amp sets a noise floor far above a naive 10 µV.
- **Reciprocity matters.** Independent transmit/receive anchors imply impossible voltage
  gain > 1; tying them to one coupling makes the link passive and puts the detection
  boundary at ~135 kHz (20 dB/m oil).

## Repository layout
```
src/   pmut_echo_sim.py     calibrated 1D link budget + receiver-noise model + results figure
       pmut_wave2d.py       2D variable-density FDTD solver + 3-panel synced animation
       make_paper_figures.py  schematic + wave montage
figures/  pmut_tube_results.png, fig_schematic.png, fig_wave_montage.png
media/    100khz.mp4, 135khz.mp4, 200khz.mp4   (wave-propagation movies)
paper/    main.tex, appendix.tex, refs.bib, main.pdf (+ style files)
*.md / prompts/   provenance: decision log, session flow, exact prompts
```

## Requirements
Python 3.11+, with:
```
pip install numpy scipy matplotlib imageio imageio-ffmpeg
```
(`ffmpeg` is needed for MP4 output; `imageio-ffmpeg` bundles it. TeX Live is needed
only to rebuild the paper.)

## Reproduce

**Results figure + link-budget reports** (seconds):
```
cd src
python3 pmut_echo_sim.py 100 135 200
# -> ../figures/pmut_tube_results.png  and a printed link budget per frequency
```

**Wave-propagation movies** (a few minutes each on one CPU core):
```
cd src
python3 pmut_wave2d.py 100     # -> 100khz.mp4 / 100khz.gif
python3 pmut_wave2d.py 135
python3 pmut_wave2d.py 200
```
Pass any centre frequency in kHz, e.g. `python3 pmut_wave2d.py 75`.

**Schematic + wave montage** (reads the rendered MP4s):
```
cd src
python3 make_paper_figures.py
```

**Rebuild the paper**:
```
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

All physical parameters live at the top of `src/pmut_echo_sim.py` (the `Config`
dataclass) and `src/pmut_wave2d.py` (the `P` class), each documented inline.

## Caveat
This is a modelling study; there is no physical experiment. Absolute SNR carries
order-of-magnitude uncertainty because the source paper's own transmit and receive
sensitivities disagree by more than 100×. The robust outputs are the relative trends,
the f² cliff, the passivity (reciprocity) constraint, and the receiver-limited regime.
See the paper's Limitations section and `decisions-log.md`.

## Acknowledgement
The model, code, figures, and paper were produced with Claude (Anthropic) acting as an
autonomous research and writing assistant under human direction.
