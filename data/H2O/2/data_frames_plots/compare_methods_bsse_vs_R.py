"""
compare_methods_bsse_vs_R.py

Purpose
-------
Compare BSSE vs R curves across different methods (SCF, MP2, CCSD_T)
for the same orientations. This helps visualize how electron correlation
affects BSSE at different distances and orientations.

Usage
-----
Terminal:
  # All eligible orientations, default 3-column grid
  python3 compare_methods_bsse_vs_R.py

  # Plot orientations 1 through 9 in a 3x3 grid
  python3 compare_methods_bsse_vs_R.py --select 1:9 --ncols 3

  # Plot specific orientations by index
  python3 compare_methods_bsse_vs_R.py --select 1 2 30 32

  # Plot orientations 10 to 20, absolute values, save
  python3 compare_methods_bsse_vs_R.py --select 10:20 --abs --save

  # All 32 orientations in a 6x6 grid
  python3 compare_methods_bsse_vs_R.py --select 1:32 --ncols 6 --save

IPython:
  run compare_methods_bsse_vs_R.py --select 1:9 --ncols 3 --abs
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple

# Use TeX Live for all text rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{amsmath} \usepackage{amssymb}",
})

def parse_selection(select_args, n_total):
    """
    Parse --select argument into a list of 0-based indices.

    Accepts:
      - A single range string  "10:20"  -> indices 10..19 (1-based -> 0-based)
      - One or more integers   "1 2 30" -> those specific 1-based positions
      - None                           -> all indices 0..n_total-1
    """
    if not select_args:
        return list(range(n_total))

    # Single range token  e.g. "10:20"
    if len(select_args) == 1 and ":" in select_args[0]:
        parts = select_args[0].split(":")
        start = int(parts[0]) - 1          # convert 1-based to 0-based
        stop  = int(parts[1])              # stop is exclusive, keep as-is
        return list(range(start, stop))

    # Explicit list of 1-based integers
    return [int(x) - 1 for x in select_args]


# Mapping from CLI value key to LaTeX y-axis label
VALUE_LABELS = {
    "bsse":    r"BSSE (Hartree)",
    "no_vmfc": r"$E_{\mathrm{int}}$ no VMFC (Hartree)",
    "vmfc":    r"$E_{\mathrm{int}}$ VMFC (Hartree)",
}


def load_method_data(
    method: str,
    basis: str,
    alpha: float,
    beta: float,
    gamma: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load R and BSSE data for a specific orientation and method.

    Returns
    -------
    R_values, BSSE_values (both as numpy arrays)
    """

    csv_path = f"rot_{method}_{basis}.csv"

    if not Path(csv_path).exists():
        return np.array([]), np.array([])

    df = pd.read_csv(csv_path)

    sel = df[
        (df["alpha"] == alpha) &
        (df["beta"] == beta) &
        (df["gamma"] == gamma)
    ].sort_values("R")

    return sel["R"].to_numpy(), sel["bsse"].to_numpy()


def find_common_orientations(
    basis: str,
    methods: List[str],
    minR: int,
    top: int
) -> List[Tuple[float, float, float]]:
    """
    Find orientations that have sufficient R coverage in ALL methods.

    Returns
    -------
    List of (alpha, beta, gamma) tuples
    """

    method_orientations = {}

    for method in methods:
        csv_path = f"rot_{method}_{basis}.csv"

        if not Path(csv_path).exists():
            continue

        df = pd.read_csv(csv_path)

        coverage = (
            df.groupby(["alpha", "beta", "gamma"])["R"]
              .nunique()
        )

        eligible = coverage[coverage >= minR]
        method_orientations[method] = set(eligible.index)

    if not method_orientations:
        print("[ERROR] No method data found!")
        return []

    common = set.intersection(*method_orientations.values())

    if not common:
        print(f"[WARN] No orientations found with >= {minR} R values in ALL methods.")
        print(f"       Trying with data from available methods only...")
        common = set.union(*method_orientations.values())

    # Rank by total R coverage across all methods
    coverage_scores = {}
    for orient in common:
        total_coverage = 0
        for method in methods:
            csv_path = f"rot_{method}_{basis}.csv"
            if Path(csv_path).exists():
                df = pd.read_csv(csv_path)
                count = len(df[
                    (df["alpha"] == orient[0]) &
                    (df["beta"] == orient[1]) &
                    (df["gamma"] == orient[2])
                ])
                total_coverage += count
        coverage_scores[orient] = total_coverage

    sorted_orients = sorted(
        coverage_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )

    return [orient for orient, _ in sorted_orients]


def plot_comparison(
    orientations: List[Tuple[float, float, float]],
    basis: str,
    methods: List[str],
    use_abs: bool,
    save: bool,
    ncols: int = 3,
) -> None:
    """
    Create comparison plots for each orientation across methods.
    Orientations are arranged in a grid with ncols columns per row.
    """

    n_orients = len(orientations)

    if n_orients == 0:
        print("[ERROR] No orientations to plot!")
        return

    ncols = min(ncols, n_orients)
    nrows = int(np.ceil(n_orients / ncols))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(6 * ncols, 5 * nrows),
        squeeze=False
    )
    # hspace scales with number of rows to prevent title/label overlap
    hspace = 0.35 + 0.08 * nrows
    fig.subplots_adjust(hspace=hspace, wspace=0.35)

    # Flatten to 1D for easy indexing; hide unused axes
    axes_flat = axes.flatten()
    for idx in range(n_orients, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    colors       = {"SCF": "C0", "MP2": "C1", "CCSD_T": "C2"}
    markers      = {"SCF": "o",  "MP2": "s",  "CCSD_T": "^"}
    method_labels = {"SCF": "SCF", "MP2": "MP2", "CCSD_T": r"CCSD(T)"}

    base_ylabel = VALUE_LABELS.get("bsse", "BSSE (Hartree)")
    ylabel = rf"$|${base_ylabel}$|$" if use_abs else base_ylabel

    for idx, (alpha, beta, gamma) in enumerate(orientations):
        ax = axes_flat[idx]

        for method in methods:
            R_vals, bsse_vals = load_method_data(method, basis, alpha, beta, gamma)

            if len(R_vals) == 0:
                continue

            if use_abs:
                bsse_vals = np.abs(bsse_vals)

            ax.plot(
                R_vals, bsse_vals,
                marker=markers.get(method, "o"),
                color=colors.get(method, "black"),
                label=method_labels.get(method, method),
                linewidth=2,
                markersize=7,
                alpha=0.8
            )

        ax.set_xlabel(r"{\large $R$} (\AA)", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(
            rf"$\alpha={alpha:g}^\circ,\;\beta={beta:g}^\circ,\;\gamma={gamma:g}^\circ$",
            fontsize=11
        )
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    if save:
        suffix = "abs" if use_abs else "signed"
        out = f"compare_methods_{basis}_top{n_orients}_bsse_vs_R_{suffix}.png"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        print(f"[SAVE] {out}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Compare BSSE vs R across methods"
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=["SCF", "MP2", "CCSD_T"],
        default=["SCF", "MP2", "CCSD_T"],
        help="Methods to compare"
    )
    parser.add_argument(
        "--basis",
        default="aug-cc-pvdz",
        help="Basis set name"
    )
    parser.add_argument(
        "--top",
        type=int,
        default=32,
        help="Maximum pool of top-ranked orientations to consider (default: 32)"
    )
    parser.add_argument(
        "--select",
        nargs="+",
        default=None,
        metavar="N",
        help=(
            "Select which orientations to plot (1-based). "
            "Use a range '10:20' or explicit indices '1 2 30 32'. "
            "Defaults to all eligible orientations up to --top."
        )
    )
    parser.add_argument(
        "--minR",
        type=int,
        default=5,
        help="Minimum number of R values required"
    )
    parser.add_argument(
        "--abs",
        action="store_true",
        help="Plot absolute values"
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save figure as PNG"
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=3,
        help="Number of subplot columns per row (default: 3)"
    )
    args = parser.parse_args()

    print(f"\n[INFO] Finding common orientations across {args.methods}...")

    all_orientations = find_common_orientations(
        args.basis,
        args.methods,
        args.minR,
        args.top,
    )

    if not all_orientations:
        print("[ERROR] No suitable orientations found for comparison.")
        print(f"        Try lowering --minR or generating more data.")
        return

    print(f"[INFO] {len(all_orientations)} eligible orientation(s) available")

    # Apply --select filter
    indices = parse_selection(args.select, len(all_orientations))
    # Clamp to valid range
    indices = [i for i in indices if 0 <= i < len(all_orientations)]

    if not indices:
        print("[ERROR] --select produced no valid indices.")
        print(f"        Valid range: 1 to {len(all_orientations)}")
        return

    orientations = [all_orientations[i] for i in indices]
    print(f"[INFO] Plotting {len(orientations)} orientation(s) (indices {[i+1 for i in indices]})")

    plot_comparison(
        orientations,
        args.basis,
        args.methods,
        args.abs,
        args.save,
        args.ncols,
    )


if __name__ == "__main__":
    main()
