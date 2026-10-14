from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
MODELS = {"full": "Full", "no_size": "No size"}
DROPOUT_LEVELS = (0.0, 0.2, 0.5)


def write_table(runs, summary):
    means = runs.groupby(["dropout_target", "model"])[
        ["beta_age", "beta_size", "beta_chemo"]
    ].mean()
    chemo = summary.loc[summary["parameter"] == "beta_chemo"].set_index(
        ["dropout_target", "model"]
    )
    rows = []
    for dropout in DROPOUT_LEVELS:
        for model, label in MODELS.items():
            coefficients = means.loc[(dropout, model)]
            performance = chemo.loc[(dropout, model)]
            values = [*coefficients, *performance[["bias", "variance", "mse"]]]
            formatted = ["---" if pd.isna(value) else f"{value:.4f}" for value in values]
            rows.append(
                f"{dropout:.0%}".replace("%", r"\%")
                + " & " + label + " & " + " & ".join(formatted) + r" \\"
            )

    table = r"""\begin{table}[htbp]
\centering
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llrrrrrr}
\toprule
$q$ & Model & $\bar\beta_{\mathrm{age}}$ & $\bar\beta_{\mathrm{size}}$ & $\bar\beta_{\mathrm{chemo}}$ & Bias & Variance & MSE \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\caption{Mean log hazard ratios over successful fits; age is per year and
tumor size per mm. Bias, variance and MSE refer to chemotherapy
($\beta=-0.7$). Here $q=P(C<5)$; the dash denotes an omitted covariate.}
\label{tab:part3_coefficients}
\end{table}
"""
    (ROOT / "tables").mkdir(exist_ok=True)
    (ROOT / "tables" / "coefficients.tex").write_text(table)


def plot_risks(runs):
    times = np.array([1, 3, 5])
    columns = [f"risk_{time}y" for time in times]
    truth = -np.expm1(-0.001 * np.exp(5.25) * times**0.7)
    styles = {"full": ("#0072B2", "o", -0.06), "no_size": ("#D55E00", "s", 0.06)}
    fig, axes = plt.subplots(1, 3, figsize=(7, 2.8), sharex=True, sharey=True)
    fig.set_facecolor("white")

    for ax, dropout in zip(axes, DROPOUT_LEVELS):
        for model, label in MODELS.items():
            risk = runs.loc[
                (runs["dropout_target"] == dropout) & (runs["model"] == model), columns
            ]
            color, marker, offset = styles[model]
            ax.errorbar(
                times + offset, risk.mean(), yerr=risk.std(ddof=0),
                color=color, marker=marker, markersize=4, linewidth=1.2,
                capsize=2, label=f"{label} (mean ± SD)",
            )
        ax.plot(times, truth, color="0.25", linestyle="--", linewidth=1.2, label="True risk")
        ax.set_title(f"$q={dropout:.0%}$".replace("%", r"\%"), fontsize=10)
        ax.set_xticks(times)
        ax.set_xlim(0.6, 5.4)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, fontsize=8)
    fig.supxlabel("Years since baseline", fontsize=9)
    fig.supylabel("Recurrence risk", fontsize=9)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.19, top=0.77, wspace=0.12)
    (ROOT / "figures").mkdir(exist_ok=True)
    for extension in ("pdf", "png"):
        fig.savefig(ROOT / "figures" / f"recurrence_risk.{extension}", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    runs = pd.read_csv(ROOT / "data" / "simulation_runs.csv")
    summary = pd.read_csv(ROOT / "data" / "simulation_summary.csv")
    write_table(runs, summary)
    plot_risks(runs)
