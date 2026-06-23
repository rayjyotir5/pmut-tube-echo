"""
make_paper_figures.py
Generate the static figures used in the paper that are not produced by
pmut_echo_sim.py (which makes pmut_tube_results.png):
  - fig_schematic.png    : the oil-filled steel-tube pulse-echo geometry
  - fig_wave_montage.png  : FDTD wave-field snapshots, 100 kHz (survives) vs
                            200 kHz (dies), cropped from the rendered animations.
Run from the spike root; reads the animation mp4s from the pmut_sim working dir.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow

OUT = os.path.join(os.path.dirname(__file__), "..", "figures")
PMUT_DIR = "/Users/jyotirmoy/Documents/pmut_sim"


# ----------------------------------------------------------------------------
def schematic():
    fig, ax = plt.subplots(figsize=(9, 3.2))
    L, Rout, wall = 3.0, 0.15, 0.03 * 0.0 + 0.012   # exaggerate wall for clarity
    wall = 0.018
    # oil column
    ax.add_patch(Rectangle((0, -Rout + wall), L, 2 * (Rout - wall),
                           fc="#cfe8ff", ec="none", zorder=1))
    # steel walls (top/bottom)
    for y0 in (-Rout, Rout - wall):
        ax.add_patch(Rectangle((0, y0), L, wall, fc="#9a9a9f", ec="k", lw=0.5, zorder=2))
    # far-end thick steel reflector
    ax.add_patch(Rectangle((L, -Rout), 0.18, 2 * Rout, fc="#6f6f74", ec="k", lw=0.5, zorder=2))
    # near-end cap + PMUT
    ax.add_patch(Rectangle((-0.04, -Rout), 0.04, 2 * Rout, fc="#9a9a9f", ec="k", lw=0.5))
    ax.add_patch(Rectangle((-0.10, -0.05), 0.06, 0.10, fc="#ffb24d", ec="k", lw=0.6))
    ax.text(-0.07, 0.13, "PMUT at\nend cap", ha="center", fontsize=8)
    # outgoing / returning arrows
    ax.add_patch(FancyArrow(0.15, 0.055, 2.4, 0, width=0.004, head_width=0.02,
                            head_length=0.12, fc="C3", ec="C3", length_includes_head=True))
    ax.text(1.3, 0.085, "transmit (volume velocity)", color="C3", fontsize=8)
    ax.add_patch(FancyArrow(2.55, -0.055, -2.4, 0, width=0.004, head_width=0.02,
                            head_length=0.12, fc="C0", ec="C0", length_includes_head=True))
    ax.text(1.3, -0.10, "echo from piston face (incident pressure)", color="C0", fontsize=8)
    # annotations
    ax.annotate("", xy=(0, -Rout - 0.03), xytext=(L, -Rout - 0.03),
                arrowprops=dict(arrowstyle="<->", color="k"))
    ax.text(L / 2, -Rout - 0.075, "oil column: piston distance (up to ~3 m)", ha="center", fontsize=9)
    ax.text(L + 0.09, 0, "piston\nface", ha="center", va="center", fontsize=8, color="white")
    ax.annotate("OD 300 mm", xy=(0.02, Rout), xytext=(0.35, Rout + 0.05),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.text(0.02, -Rout + wall + 0.004, " 10 mm wall", fontsize=7, va="bottom")
    ax.set_xlim(-0.2, 3.35); ax.set_ylim(-0.27, 0.2)
    ax.axis("off")
    fig.tight_layout()
    p = os.path.join(OUT, "fig_schematic.png")
    fig.savefig(p, dpi=150, bbox_inches="tight"); print("saved", p)


# ----------------------------------------------------------------------------
def wave_montage():
    import imageio.v2 as imageio
    sim_dt = 4.70   # ms total
    times = [0.30, 1.50, 4.30]
    rows = [("100khz.mp4", "100 kHz (survives)"), ("200khz.mp4", "200 kHz (dies in transit)")]
    fig, axes = plt.subplots(len(rows), len(times), figsize=(12, 3.4))
    for r, (fn, lab) in enumerate(rows):
        path = os.path.join(PMUT_DIR, fn)
        rd = imageio.get_reader(path); frames = [f for f in rd]; n = len(frames)
        for c, tt in enumerate(times):
            fr = frames[min(n - 1, int(tt / sim_dt * n))]
            h, w = fr.shape[:2]
            # crop the wave-field (top) panel
            crop = fr[int(0.10 * h):int(0.40 * h), int(0.06 * w):int(0.99 * w)]
            ax = axes[r, c]
            ax.imshow(crop); ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(f"t = {tt:.2f} ms", fontsize=9)
            if c == 0:
                ax.set_ylabel(lab, fontsize=9)
    fig.suptitle("FDTD acoustic field in the oil-filled tube (axial slice): low frequency survives the "
                 "3 m round trip, high frequency is absorbed", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUT, "fig_wave_montage.png")
    fig.savefig(p, dpi=140, bbox_inches="tight"); print("saved", p)




# ----------------------------------------------------------------------------
def hero():
    """Opening teaser: the 100 kHz pulse launches at the end cap, travels 3 m down
    the oil column, reflects off the piston face, and returns (FDTD field)."""
    import imageio.v2 as imageio
    sim_dt = 4.70
    stages = [(0.30, "launch"), (1.55, "mid-transit"), (4.32, "echo returns")]
    path = os.path.join(PMUT_DIR, "100khz.mp4")
    rd = imageio.get_reader(path); frames = [f for f in rd]; n = len(frames)
    fig, axes = plt.subplots(len(stages), 1, figsize=(11, 4.3))
    for ax, (tt, lab) in zip(axes, stages):
        fr = frames[min(n - 1, int(tt / sim_dt * n))]
        h, w = fr.shape[:2]
        # crop just the wave-field (drop panel title above and x-axis below)
        crop = fr[int(0.135 * h):int(0.385 * h), int(0.085 * w):int(0.995 * w)]
        ax.imshow(crop); ax.set_xticks([]); ax.set_yticks([])
        ax.set_ylabel(f"{lab}\n(t={tt:.1f} ms)", fontsize=9, rotation=0,
                      ha="right", va="center", labelpad=28)
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle("A 100 kHz pulse launched at the end cap travels the 3 m oil column, "
                 "reflects off the piston face, and returns",
                 fontsize=11, y=0.99)
    fig.tight_layout(rect=[0.02, 0, 1, 0.95])
    p = os.path.join(OUT, "fig_hero.png")
    fig.savefig(p, dpi=150, bbox_inches="tight"); print("saved", p)


if __name__ == "__main__":
    schematic()
    try:
        hero()
    except Exception as e:
        print("hero skipped:", e)
    try:
        wave_montage()
    except Exception as e:
        print("wave_montage skipped:", e)
