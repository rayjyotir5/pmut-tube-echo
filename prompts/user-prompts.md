# Verbatim human prompts (this session)

These are the human inputs, in order. They are the sole human-authored steering
of the project; everything else (derivations, code, figures, paper) was produced
by the AI agent (Claude) in response. Light typos preserved as written.

1. "read the @modelling.pdf and understand what it is doing - will specify the actual problem after that"

2. "Has analytical model shown acoustic transmission through 10 mm steel + synthetic mineral oil
   yields detectable echo at 3 m with SNR >= 12 dB at 1 MHz centre frequency? if not i want to set
   up a simulation and visualize the results"
   - Follow-up multiple-choice answers (agent-posed clarifying questions):
     * Geometry: "Pulse-echo through wall"
     * Noise: "Sweep / show sensitivity"
     * Method: "1D analytical + transfer-matrix"

3. "lets now upgrade the sim to a 2d environment - a steel tube is 3m long, 10mm thick and since it
   is the cross section of a circular tube, the outer dia of the tube is 300mm. the pmut is
   positioned at one end and the other end forms the reflecting surface (this end may have increased
   thickness of steel so as to ignore steel-air interface. tube is filled with oil. now run the
   simulation and repeat previous sweeps. feel free to re-center the frequency sweep and not assume
   1MHz is a fixed input so that SNR is optimal. i'd also like to see a video or gif animation of the
   simulation of the wave propagating - use https://github.com/0x23/WaveSimulator2D approach if required."

4. "oil attenuation is 20 db/m and waveguide wall-bounce loss of 0.5db/m is fine. in the mp4, please
   render two more panels below the wave propagation animation. one is the pmut transmitter signal
   over time plot and a vertical line moves along that synced to wave animation. similarly a panel
   for measured pmut receiver signal over time. also have a legend to show physical time of sim
   progress. hope the wave is simulated with noise considered at 250khz"

5. "erm 20db/m was not at 1mhz, it is at 250khz"

6. "render a second clip at the ~50-100 kHz optimum so you can see a clean, surviving echo for contrast"

7. "is the sharp drop at 1MHz in the SNR vs frequency correct? why is this happening? would it help
   to now stop the sweep and model the noise correctly?"

8. "yes"   (proceed to tie transmit and receive to one reciprocity-consistent coupling)

9. "lets now regenerate the results plots and wave propagation animations with the new calibrated
   model and view it above and below 135kHz"

10. "make sure to regenerate pmut_tube_results.png"

11. "i dont want three png files, i want a single one with pmut_tube_results.png"

12. "in pmut_tube_results.png the plots are at 250khz and title says ... pls update to 100, 135, 200"

13. "plot y axes are not scaled correctly - getting cropped"

14. "panel 1 and 4 as well...will log scale help?"

15. "i dont like the bar chart in panel 4...id rather see the waveform, not normalized to noise"

16. "read voila.md and treat this entire session as the research done - accordingly document
    everything and produce the resultant paper along with all assumptions, mathematical models,
    simulations and results"   [the request that produced this spike]
