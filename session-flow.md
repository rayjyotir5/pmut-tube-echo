# Session flow (human input vs AI effort) — for paper appendix

Provenance record. The human supplied the starting paper, the physical scenario,
a handful of physical parameters, and steering/quality feedback. The AI agent
(Claude, Opus 4.8) did all reading, derivation, coding, debugging, figure
generation, and writing. Each turn below pairs the human steer with the AI work.

## Turn 1 — starting point
- HUMAN: read modelling.pdf (Dangi & Pratap 2017, system-level PMUT model), will specify problem.
- AI: read the 11-page paper; summarized the lumped/analytical transducer model, the kappa^2
  plate-vs-membrane regimes, the transmit/receive/transceiver transfer functions, and validation.

## Turn 2 — the question + first simulator
- HUMAN: does the paper's model show a detectable echo (SNR >= 12 dB) through 10 mm steel + mineral
  oil at 3 m, 1 MHz? If not, build a simulation. (Chose: pulse-echo through wall; sweep noise;
  1D analytical + transfer matrix.)
- AI: established the paper does NOT cover this (it radiates into air, <0.5 m, no SNR). Built a
  free-field link-budget simulator: transfer-matrix steel wall, beam spreading, oil attenuation,
  target reflection, SNR-vs-noise sweep. Found two bugs in own first draft (double-counted spreading;
  unphysical 7 W transmit) and fixed them. Verdict: conditional / marginal in free field.

## Turn 3 — geometry change to the oil-filled tube + FDTD animation
- HUMAN: change to a 3 m oil-filled steel tube (OD 300 mm, 10 mm wall), PMUT at one end, far end a
  thick-steel reflector; re-center frequency for best SNR; produce a wave-propagation animation
  (WaveSimulator2D approach).
- AI: reformulated the channel as a waveguide (spherical spreading -> bounded fill factor; far end
  near-total reflector). Built a 2D variable-density staggered-grid FDTD solver (fixed a sign-error
  instability). Found lower frequency strongly preferred.

## Turn 4 — parameters + richer animation
- HUMAN: oil attenuation 20 dB/m, wall-bounce 0.5 dB/m; add transmitter and receiver time-signal
  panels with a synced cursor and a sim-time readout; include noise at 250 kHz.
- AI: added per-cell oil attenuation to the FDTD; built the 3-panel synced animation.

## Turn 5 — attenuation reference correction
- HUMAN: 20 dB/m is at 250 kHz, not 1 MHz.
- AI: re-referenced attenuation to 250 kHz (f^2 scaling). This is a 112 dB change at the operating
  point; the verdict flips from "easy" to "marginal".

## Turn 6 — contrast clip
- HUMAN: render a second clip near the 50-100 kHz optimum (surviving echo).
- AI: rendered a low-frequency clip; the wave survives the round trip.

## Turn 7 — diagnosis question
- HUMAN: is the 1 MHz SNR cliff correct? why? should we model the noise properly instead of sweeping?
- AI: confirmed the cliff is the f^2 absorption law over the 6 m round trip (signal collapses as
  exp(-k f^2)); implemented a first-principles charge-amplifier + PZT-capacitance noise model,
  replacing the arbitrary floor. Discovered the receiver is voltage-noise/capacitance limited.

## Turn 8 — reciprocity fix
- HUMAN: "yes" (tie transmit and receive to one coupling).
- AI: found the paper's measured transmit (3 nm/V) and Eq. 17 receive sensitivity disagree by ~41x.
  Anchored both ends to one coupling via the lumped reciprocity relation qdot/V = w00 * (Q/Pin),
  removing an impossible gain>1 artifact (140 V echo) and making the model passive and area-consistent.

## Turn 9-15 — calibrated re-runs + figure quality
- HUMAN: regenerate plots/animations with the calibrated model; view above and below the 135 kHz
  detection boundary; single results PNG; operating points 100/135/200 kHz; fix cropped/auto-scaled
  axes; panel 4 should show the real waveform (not a bar chart, not normalized).
- AI: produced the single comparison figure (axes bounded so the exp(-k f^2) tail does not crush the
  plot), and three FDTD clips at 100/135/200 kHz with comparable receiver panels showing the actual
  received waveform on per-frequency scales.

## Turn 16 — documentation + paper (this deliverable)
- HUMAN: treat the whole session as the research; document everything; produce the paper.
- AI: created this spike directory, the logs, and the CAISc-2026-format paper (main + appendix) with
  compiled PDF.
