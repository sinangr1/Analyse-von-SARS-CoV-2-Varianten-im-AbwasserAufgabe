#!/usr/bin/env python3
"""Create three labeled SARS-CoV-2 wastewater plots from LolliPop output.

Place this script at Abgabe/scripts/plot_with_percentages.py.
Expected input: Abgabe/results/deconvoluted.tsv
Alternative: Abgabe/results/deconvoluted.tsv.zst
Output: Abgabe/plots/*_labeled.png

Install requirements: pandas, numpy, matplotlib (and zstandard for .zst input).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
PLOTS = ROOT / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

VARIANTS = [
    "BA.1", "BA.2", "BA.2.12.1", "BA.2.75",
    "BA.4", "BA.5", "BQ.1", "undetermined",
]
COLORS = dict(zip(VARIANTS, plt.get_cmap("tab10").colors[:8]))


def load_results():
    candidates = [RESULTS / "deconvoluted.tsv", RESULTS / "deconvoluted.tsv.zst"]
    source = next((path for path in candidates if path.is_file()), None)
    if source is None:
        raise FileNotFoundError(
            "Missing input: place deconvoluted.tsv or deconvoluted.tsv.zst "
            "in the Abgabe/results folder."
        )

    df = pd.read_csv(source, sep="\t")
    required = {"location", "variant", "date", "proportion"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")
    df["date"] = pd.to_datetime(df["date"])
    df["proportion"] = pd.to_numeric(df["proportion"], errors="raise")
    print(f"INPUT: {source}")
    return df


def plot_location(df, location, filename):
    sub = df.loc[df["location"] == location]
    if sub.empty:
        raise ValueError(f"No data for location: {location}")

    table = sub.pivot(index="date", columns="variant", values="proportion")
    table = (table.reindex(columns=VARIANTS, fill_value=0).fillna(0) * 100)
    table = table.sort_index()
    dates = table.index.strftime("%d.%m.%Y")
    x = np.arange(len(table))

    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 2.4], hspace=0.36)
    ax = fig.add_subplot(gs[0])
    ax_table = fig.add_subplot(gs[1])
    ax_table.axis("off")

    bottom = np.zeros(len(table))
    for variant in VARIANTS:
        values = table[variant].to_numpy()
        ax.bar(x, values, bottom=bottom, width=0.66, color=COLORS[variant], label=variant)
        for i, value in enumerate(values):
            if value >= 4.0:
                ax.text(
                    x[i], bottom[i] + value / 2, f"{value:.1f} %",
                    ha="center", va="center", fontsize=10, fontweight="bold",
                    color="white" if variant in {"BA.1", "BA.2.12.1", "BA.2.75", "BA.4", "BA.5"} else "black",
                )
        bottom += values

    ax.set_title(f"{location}: SARS-CoV-2-Varianten im Abwasser", fontsize=16, pad=15)
    ax.set_ylabel("Geschätzter Variantenanteil (%)")
    ax.set_xlabel("Probenahmedatum")
    ax.set_xticks(x)
    ax.set_xticklabels(dates)
    ax.set_ylim(0, 100)
    ax.legend(title="Variante", loc="upper left", bbox_to_anchor=(1.02, 1))

    cell_text = [[f"{v:.2f}" for v in table[variant]] for variant in VARIANTS]
    tab = ax_table.table(
        cellText=cell_text, rowLabels=VARIANTS, colLabels=dates,
        cellLoc="center", rowLoc="center", loc="center", bbox=[0, 0.05, 1, 0.9],
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(10)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if r == 0:
            cell.set_facecolor("#eeeeee")
            cell.set_text_props(weight="bold")

    fig.text(
        0.12, 0.025,
        "Alle Tabellenwerte in %. LolliPop-Schätzungen; frühe BQ.1-Zuordnungen sind unsicher.",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.14, right=0.80, top=0.93, bottom=0.07)
    output = PLOTS / filename
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"CREATED: {output}")


def plot_ba5_comparison(df):
    sub = df.loc[df["variant"] == "BA.5"].copy()
    sub["percent"] = sub["proportion"] * 100
    fig, ax = plt.subplots(figsize=(12, 7))

    for location, color, offset in [
        ("Davos", "tab:blue", 12),
        ("S-chanf", "tab:orange", -18),
    ]:
        part = sub[sub["location"] == location].sort_values("date")
        ax.plot(
            part["date"], part["percent"], marker="o", markersize=9,
            linewidth=2.5, color=color, label=location,
        )
        for _, row in part.iterrows():
            ax.annotate(
                f'{row["percent"]:.2f} %',
                (row["date"], row["percent"]),
                xytext=(0, offset), textcoords="offset points", ha="center",
                fontsize=10, fontweight="bold", color=color,
            )

    ax.set_title("BA.5: Vergleich Davos und S-chanf", fontsize=16)
    ax.set_xlabel("Probenahmedatum")
    ax.set_ylabel("Geschätzter BA.5-Anteil (%)")
    ax.set_ylim(0, 105)
    dates = sorted(sub["date"].unique())
    ax.set_xticks(dates)
    ax.set_xticklabels(pd.to_datetime(dates).strftime("%d.%m.%Y"))
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    output = PLOTS / "ba5_comparison_labeled.png"
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"CREATED: {output}")


def main():
    df = load_results()
    plot_location(df, "Davos", "davos_variants_labeled.png")
    plot_location(df, "S-chanf", "schanf_variants_labeled.png")
    plot_ba5_comparison(df)
    print("ALL THREE LABELED CHARTS COMPLETED")


if __name__ == "__main__":
    main()
