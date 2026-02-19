"""
plot_rot_bsse_vs_R.py

Purpose
-------
Select the top N rotational orientations (alpha, beta, gamma) with the
largest number of distinct R values, then plot a chosen quantity
(bsse, no_vmfc, or vmfc) vs R for each orientation on the same figure.

Supports SCF, MP2, and CCSD(T) methods.

Why this matters
----------------
If BSSE is weakly dependent on orientation, the BSSE(R) curves for
different (alpha, beta, gamma) should be similar and differences should
shrink at larger R.

Usage
-----
Terminal:
  # SCF plots
  python3 plot_rot_bsse_vs_R.py --method SCF --top 5 --minR 5
  
  # MP2 plots
  python3 plot_rot_bsse_vs_R.py --method MP2 --top 5 --minR 5 --abs
  
  # CCSD(T) plots
  python3 plot_rot_bsse_vs_R.py --method CCSD_T --top 7 --minR 5 --abs --save
  
  # Compare all methods (creates 3 separate plots)
  python3 plot_rot_bsse_vs_R.py --method all --top 5 --minR 5

IPython:
  run plot_rot_bsse_vs_R.py --method MP2 --top 5 --minR 5 --abs
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Use TeX Live for all text rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{amsmath} \usepackage{amssymb}",
})

# Mapping from CLI value key to LaTeX y-axis label
VALUE_LABELS = {
    "bsse":    r"BSSE (Hartree)",
    "no_vmfc": r"$E_{\mathrm{int}}$ no VMFC (Hartree)",
    "vmfc":    r"$E_{\mathrm{int}}$ VMFC (Hartree)",
}


def plot_method_data(
    csv_path: str,
    method: str,
    top: int,
    minR: int,
    value_col: str,
    use_abs: bool,
    save: bool
) -> None:
    """
    Plot BSSE vs R for top orientations for a single method.

    Parameters
    ----------
    csv_path : str
        Path to the CSV file
    method : str
        Method name (for labels)
    top : int
        Number of top orientations to plot
    minR : int
        Minimum number of R values required
    value_col : str
        Column to plot ('bsse', 'no_vmfc', 'vmfc')
    use_abs : bool
        Whether to plot absolute values
    save : bool
        Whether to save the figure
    """

    if not Path(csv_path).exists():
        print(f"[SKIP] {csv_path} not found. Skipping {method}.")
        return

    df = pd.read_csv(csv_path)

    # Count distinct R per orientation
    coverage = (
        df.groupby(["alpha", "beta", "gamma"])["R"]
          .nunique()
          .sort_values(ascending=False)
    )

    eligible = coverage[coverage >= minR]
    if eligible.empty:
        print(f"[INFO] {method}: No orientations found with >= {minR} distinct R values.")
        return

    picked = eligible.head(top)
    print(f"\n[INFO] {method}: Plotting {len(picked)} orientations")

    fig, ax = plt.subplots(figsize=(10, 6))

    base_ylabel = VALUE_LABELS.get(value_col, f"{value_col} (Hartree)")
    ylabel = rf"$|${base_ylabel}$|$" if use_abs else base_ylabel

    # Use a color cycle for better distinction
    colors = plt.cm.tab10(np.linspace(0, 1, len(picked)))

    for idx, ((a, b, g), nR) in enumerate(picked.items()):
        sel = (
            df[(df["alpha"] == a) & (df["beta"] == b) & (df["gamma"] == g)]
            .sort_values("R")
        )

        xvals = sel["R"].to_numpy()
        yvals = sel[value_col].to_numpy()

        if use_abs:
            yvals = np.abs(yvals)

        label = (
            rf"$\alpha={a:g}^\circ,\;\beta={b:g}^\circ,\;\gamma={g:g}^\circ$"
            rf" ($n_R={int(nR)}$)"
        )
        ax.plot(xvals, yvals, marker="o", label=label,
                color=colors[idx], linewidth=2, markersize=6)

    ax.set_xlabel(r"{\large $R$} (\AA)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(
        rf"{method}: {base_ylabel} vs $R$ — top {len(picked)} orientations",
        fontsize=13
    )
    ax.legend(fontsize=9, loc="best")
    #ax.grid(True, alpha=0.3)  # Uncomment to show grid lines
    plt.tight_layout()

    if save:
        suffix = "abs" if use_abs else "signed"
        basis = Path(csv_path).stem.split("_")[-1]
        out = f"rot_{method}_{basis}_top{len(picked)}_{value_col}_vs_R_{suffix}.png"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[SAVE] {out}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot BSSE vs R for rotational orientations"
    )
    parser.add_argument(
        "--method",
        choices=["SCF", "MP2", "CCSD_T", "all"],
        default="SCF",
        help="Correlation method to plot (or 'all' for all methods)"
    )
    parser.add_argument(
        "--basis",
        default="aug-cc-pvdz",
        help="Basis set name"
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="Number of top orientations to plot"
    )
    parser.add_argument(
        "--minR",
        type=int,
        default=5,
        help="Minimum number of R values required for an orientation"
    )
    parser.add_argument(
        "--value",
        choices=["bsse", "no_vmfc", "vmfc"],
        default="bsse",
        help="Quantity to plot"
    )
    parser.add_argument(
        "--abs",
        action="store_true",
        help="Plot absolute value of the selected quantity"
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save figures as PNG files"
    )
    args = parser.parse_args()

    methods_to_plot = ["SCF", "MP2", "CCSD_T"] if args.method == "all" else [args.method]

    for method in methods_to_plot:
        csv_name = f"rot_{method}_{args.basis}.csv"
        plot_method_data(
            csv_path=csv_name,
            method=method,
            top=args.top,
            minR=args.minR,
            value_col=args.value,
            use_abs=args.abs,
            save=args.save
        )


if __name__ == "__main__":
    main()
