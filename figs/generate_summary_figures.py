"""Generates the figures used in SUMMARY.md.

Uses illustrative parameters (see NOTES.md: cost, lysis rate, and environmental
switching rate are all unmeasured for this system) -- chosen to be qualitatively
reasonable, not fit to data. Re-run after changing lib/fitness.py or sim/ to
refresh the figures.
"""
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sweeps"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from run_phase_diagram import run_phase_diagram  # noqa: E402

FIGS_DIR = os.path.dirname(__file__)

# Palette (see dataviz skill reference/palette.md): diverging blue<->red, neutral
# gray midpoint; categorical slot-1 blue for a single-hue binary fill.
BLUE = "#2a78d6"
RED = "#e34948"
NEUTRAL_GRAY = "#f0efec"
MUTED_GRAY = "#c9c8c2"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"

plt.rcParams.update({
    "font.size": 11,
    "text.color": TEXT_PRIMARY,
    "axes.edgecolor": MUTED_GRAY,
    "axes.labelcolor": TEXT_PRIMARY,
    "xtick.color": TEXT_SECONDARY,
    "ytick.color": TEXT_SECONDARY,
    "axes.grid": False,
})

# Illustrative baseline parameters -- see SUMMARY.md / NOTES.md for what's
# established vs. assumed. Nothing here is fit to experimental data.
G_A, G_N = 1.0, 0.4      # aerobic growth rate > anaerobic (TMAO) growth rate
Q_A, Q_N = 0.1, 0.95     # minority "prepared" aerobically; nearly all anaerobically
C = 0.1                  # modest, hard-to-detect cost of needless expression
D = 5.0                  # strong growth arrest if caught unprepared
K_SWITCH = 1.0           # phenotype switches about as fast as cells divide
M, DELTA, BETA, S_STAR = 1.0, 0.001, 3.0, 10.0  # illustrative phage parameters


def make_phase_diagram():
    sigma_values = np.logspace(-3, 1, 25)
    lambda0_values = np.linspace(0.0, 1.2, 25)

    df = run_phase_diagram(
        sigma_values=list(sigma_values), k_switch_values=[K_SWITCH],
        c_values=[C], lambda0_values=list(lambda0_values), alpha_values=[0.0],
        q_A=Q_A, q_N=Q_N, g_A=G_A, g_N=G_N, D=D,
        m=M, delta=DELTA, beta=BETA, S_star=S_STAR,
        target_env_segments=3000, env_seed=1, n_jobs=-1,
    )

    host_grid = df.pivot(index="lambda0", columns="sigma_AN", values="host_gap").values
    phage_grid = df.pivot(index="lambda0", columns="sigma_AN", values="phage_invasion_rate").values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)

    vmax = np.abs(host_grid).max()
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "diverging", [BLUE, NEUTRAL_GRAY, RED])
    im1 = ax1.pcolormesh(sigma_values, lambda0_values, host_grid, cmap=cmap,
                          vmin=-vmax, vmax=vmax, shading="nearest")
    ax1.set_xscale("log")
    ax1.set_xlabel("environmental switching rate\n(stable, host-like $\\rightarrow$ volatile, free-living-like)")
    ax1.set_ylabel("cost of carrying the phage\n(lysis rate)")
    ax1.set_title("Which strategy grows faster?", loc="left", fontsize=12)
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.set_ticks([-vmax * 0.85, vmax * 0.85])
    cbar1.set_ticklabels(["favors\nbet-hedging", "favors\non/off switching"])

    phage_bool = (phage_grid > 0).astype(float)
    cmap2 = matplotlib.colors.ListedColormap([MUTED_GRAY, BLUE])
    ax2.pcolormesh(sigma_values, lambda0_values, phage_bool, cmap=cmap2,
                   vmin=0, vmax=1, shading="nearest")
    ax2.set_xscale("log")
    ax2.set_xlabel("environmental switching rate\n(stable, host-like $\\rightarrow$ volatile, free-living-like)")
    ax2.set_title("Does the phage spread when rare?", loc="left", fontsize=12)
    legend_elems = [
        matplotlib.patches.Patch(facecolor=BLUE, label="phage spreads"),
        matplotlib.patches.Patch(facecolor=MUTED_GRAY, label="phage fails to spread"),
    ]
    ax2.legend(handles=legend_elems, loc="upper right", frameon=False, fontsize=9)

    fig.suptitle("Host and phage don't necessarily agree on the outcome",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    out = os.path.join(FIGS_DIR, "phase_diagram.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def make_crossover_line():
    # Lysis rate is the axis with the cleanest, most robust story (monotonic,
    # crosses zero, and doesn't depend on how well-matched the bet-hedging
    # "prepared fraction" q_A is to the environment's switching rate -- see
    # SUMMARY.md's caveat on the sigma axis for why that one is more subtle).
    lambda0_values = np.linspace(0.0, 0.6, 40)
    sigma_reps = [0.01, 1.0]

    df = run_phase_diagram(
        sigma_values=sigma_reps, k_switch_values=[K_SWITCH],
        c_values=[C], lambda0_values=list(lambda0_values), alpha_values=[0.0],
        q_A=Q_A, q_N=Q_N, g_A=G_A, g_N=G_N, D=D,
        m=M, delta=DELTA, beta=BETA, S_star=S_STAR,
        target_env_segments=3000, env_seed=1, n_jobs=-1,
    )

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    colors = [BLUE, RED]
    for color, sigma in zip(colors, sigma_reps):
        sub = df[df["sigma_AN"] == sigma].sort_values("lambda0")
        label = f"stable environment (sigma={sigma})" if sigma < 1 else f"volatile environment (sigma={sigma})"
        ax.plot(sub["lambda0"], sub["host_gap"], color=color, lw=2.5, label=label)

    ax.axhline(0, color=MUTED_GRAY, lw=1.5, zorder=0)
    ax.set_xlabel("cost of carrying the phage (lysis rate)")
    ax.set_ylabel("growth-rate advantage of\non/off switching over bet-hedging")
    ax.text(0.02, 0.95, "switching wins", transform=ax.transAxes,
            color=TEXT_SECONDARY, fontsize=9, va="top")
    ax.text(0.02, 0.05, "bet-hedging wins", transform=ax.transAxes,
            color=TEXT_SECONDARY, fontsize=9, va="bottom")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title("A cheap phage favors switching; a costly one favors bet-hedging", loc="left", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIGS_DIR, "crossover_line.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    make_phase_diagram()
    make_crossover_line()
