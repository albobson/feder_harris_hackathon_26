"""Generates the figures used in SUMMARY.md.

Uses illustrative parameters (see NOTES.md: the expression cost, the lysis rate and
the environmental switching rate are all unmeasured for this system) -- chosen to be
qualitatively reasonable, not fit to data. Re-run after changing lib/ or sim/.

Efficiency note: lysis is a uniform tax on the lysogen (it hits both phenotypes
equally), so host_gap(sigma, lambda) == growth_advantage(sigma) - lambda_bar
EXACTLY. So the expensive Lyapunov-exponent calculation only has to run once per
sigma, and the lambda axis is then derived analytically rather than swept.
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
from fitness import phage_invasion_growth_rate  # noqa: E402

FIGS_DIR = os.path.dirname(__file__)

# Palette roles (see the dataviz reference palette).
BLUE, RED, AQUA = "#2a78d6", "#e34948", "#1baf7a"
NEUTRAL, MUTED = "#f0efec", "#c9c8c2"
INK, INK_2 = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "font.size": 11, "text.color": INK, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.grid": False,
})

# --- Illustrative baseline parameters (nothing here is fit to data) ---
G_A, G_N = 1.0, 0.4     # aerobic growth rate > anaerobic (TMAO) growth rate
Q_A_S, Q_A_L = 0.1, 0.0  # aerobic prepared fraction: bet-hedger vs lysogen (off)
Q_N = 0.95              # anaerobic target, the same for both genotypes
C = 0.1                 # cost of making torCAD when it isn't needed yet
D = 5.0                 # growth arrest when caught unprepared by anaerobiosis
K_SWITCH = 1.0          # phenotype switches about as fast as cells divide
M, DELTA, BETA, S_STAR = 1.0, 0.01, 4.0, 50.0   # illustrative phage parameters

SIGMAS = np.logspace(-3, 1.5, 130)


def growth_advantage_by_sigma():
    """The lysogen's growth-rate advantage over the bet-hedger, before lysis."""
    df = run_phase_diagram(
        sigma_values=list(SIGMAS), k_switch_values=[K_SWITCH], c_values=[C],
        lambda0_values=[0.0], alpha_values=[0.0],
        q_A_S=Q_A_S, q_A_L=Q_A_L, q_N=Q_N, g_A=G_A, g_N=G_N, D=D,
        m=M, delta=DELTA, beta=BETA, S_star=S_STAR,
        target_env_segments=6000, env_seed=1, n_jobs=-1)
    df = df.sort_values("sigma_AN")
    return df["sigma_AN"].values, df["growth_advantage"].values


def make_crossover_line(sigmas, adv):
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for color, lam, label in [(BLUE, 0.0, "no lysis cost"),
                               (AQUA, 0.002, "lysis rate 0.002"),
                               (RED, 0.01, "lysis rate 0.01")]:
        ax.plot(sigmas, adv - lam, color=color, lw=2.5, label=label)

    ax.axhline(0, color=MUTED, lw=1.5, zorder=0)
    ax.set_xscale("log")
    ax.set_xlabel("how often oxygen availability flips\n"
                  "(stable, host-associated $\\rightarrow$ volatile, free-living)")
    ax.set_ylabel("growth advantage of on/off switching\nover bet-hedging")
    # Anchor the two region labels to the zero line (data y, axes x) so they sit
    # beside it wherever it falls, instead of colliding with the curves.
    blend = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.985, 0.0012, "on/off switching wins", transform=blend, color=INK_2,
            fontsize=9, ha="right", va="bottom")
    ax.text(0.985, -0.0015, "bet-hedging wins", transform=blend, color=INK_2,
            fontsize=9, ha="right", va="top")
    ax.legend(frameon=False, fontsize=9, loc="center left")
    ax.set_title("Losing the hedge only pays in a stable environment,\n"
                 "and only if the phage is almost free to carry",
                 loc="left", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIGS_DIR, "crossover_line.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def make_phase_diagram(sigmas, adv):
    lams = np.logspace(-4, 0, 400)
    SIG, LAM = np.meshgrid(sigmas, lams)
    ADV = np.broadcast_to(adv, SIG.shape)

    host_gap = ADV - LAM                       # exact: lysis is a uniform tax
    phage = np.vectorize(phage_invasion_growth_rate)(ADV, LAM, M, DELTA, S_STAR, BETA)

    host_wins = host_gap > 0
    phage_wins = phage > 0

    # 0 = phage can't establish; 1 = conflict (phage spreads, host worse off);
    # 2 = aligned (both gain).
    region = np.where(~phage_wins, 0, np.where(host_wins, 2, 1))

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    cmap = matplotlib.colors.ListedColormap([MUTED, RED, AQUA])
    ax.pcolormesh(SIG, LAM, region, cmap=cmap, vmin=0, vmax=2, shading="nearest")
    ax.contour(SIG, LAM, host_gap, levels=[0], colors=[INK], linewidths=1.4)
    ax.contour(SIG, LAM, phage, levels=[0], colors=[INK], linewidths=1.4,
               linestyles="--")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("how often oxygen availability flips\n"
                  "(stable, host-associated $\\rightarrow$ volatile, free-living)")
    ax.set_ylabel("cost of carrying the phage (lysis rate)")
    ax.set_title("Almost everywhere the phage can spread,\n"
                 "it spreads at the bacterium's expense", loc="left", fontsize=12)

    ax.text(0.30, 0.72, "PHAGE SPREADS,\nBACTERIUM WORSE OFF", transform=ax.transAxes,
            color="white", fontsize=11, fontweight="bold", ha="center", va="center")
    ax.text(0.12, 0.16, "both gain", transform=ax.transAxes, color=INK,
            fontsize=9.5, ha="center")
    ax.text(0.80, 0.09, "phage cannot establish", transform=ax.transAxes,
            color=INK_2, fontsize=9.5, ha="center")

    fig.tight_layout()
    out = os.path.join(FIGS_DIR, "phase_diagram.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    spreads = phage_wins.sum()
    print(f"  conflict = {100.0 * np.mean(region == 1):.0f}% of the plotted plane; "
          f"of the area where the phage CAN spread, "
          f"{100.0 * (region == 1).sum() / max(spreads, 1):.0f}% has the host worse off")


if __name__ == "__main__":
    sigmas, adv = growth_advantage_by_sigma()
    print(f"growth advantage ranges {adv.min():+.4f} .. {adv.max():+.4f}; "
          f"crosses zero near sigma = "
          f"{sigmas[np.argmin(np.abs(adv))]:.3g}")
    make_crossover_line(sigmas, adv)
    make_phase_diagram(sigmas, adv)
