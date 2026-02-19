"""
plot_trnl_contours.py

Purpose
-------
Make x-z contour plots at fixed y for translational BSSE scans.
Grid completeness analysis has been moved to analyze_trnl_grid_coverage.py.

Usage
-----
  python3 plot_trnl_contours.py --method SCF --value bsse --save
  python3 plot_trnl_contours.py --method MP2 --value vmfc
  python3 plot_trnl_contours.py --method CCSD_T --value bsse --save

  # Check grid completeness first:
  python3 analyze_trnl_grid_coverage.py
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Use TeX Live for all text rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{amsmath} \usepackage{amssymb}",
})

# Mapping from CLI value key to LaTeX label for colorbar
VALUE_LABELS = {
    "bsse":    r"BSSE (Hartree)",
    "no_vmfc": r"$E_{\mathrm{int}}$ no VMFC (Hartree)",
    "vmfc":    r"$E_{\mathrm{int}}$ VMFC (Hartree)",
}


def contour_fixed_y(
    df: pd.DataFrame,
    y_values=(2.0, 2.25, 2.5, 2.75),
    value_col="bsse",
    levels=30,
    save=False,
    prefix="",
) -> None:
    """Make x-z contour plots at fixed y. Assumes these y slices are full grids."""

    cbar_label = VALUE_LABELS.get(value_col, f"{value_col} (Hartree)")

    for y0 in y_values:
        dfy = df[df["y"] == y0].copy()

        # pivot to rectangular grid: rows = x, cols = z
        grid = dfy.pivot(index="x", columns="z", values=value_col)

        missing = int(grid.isna().sum().sum())
        if missing != 0:
            print(f"[WARN] y={y0}: grid has {missing} missing points; skipping.")
            continue

        X, Z = np.meshgrid(grid.columns.to_numpy(), grid.index.to_numpy())
        V = grid.to_numpy()

        plt.figure()
        plt.contourf(X, Z, V, levels=levels)
        plt.colorbar(label=cbar_label)
        plt.xlabel(r"{\Large{$z$}} (\AA)")
        plt.ylabel(r"{\Large{$x$}} (\AA)")
        plt.title(
            rf"{cbar_label.split('(')[0].strip()} contour at "
            rf"$y = {y0}\,\mathrm{{\AA}}$ ({prefix})"
        )
        plt.tight_layout()

        if save:
            fname = f"contour_xz_y{y0:g}_{value_col}_{prefix}.png"
            plt.savefig(fname, dpi=200)
            print(f"[SAVE] {fname}")

        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot translational BSSE contours at fixed y slices."
    )
    parser.add_argument(
        "--method",
        choices=["SCF", "MP2", "CCSD_T"],
        default="SCF",
    )
    parser.add_argument(
        "--basis",
        default="aug-cc-pvdz",
    )
    parser.add_argument(
        "--value",
        choices=["bsse", "no_vmfc", "vmfc"],
        default="bsse",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save PNGs to disk.",
    )
    args = parser.parse_args()

    csv_name = f"trnl_{args.method}_{args.basis}.csv"
    df = pd.read_csv(csv_name)

    print(f"[INFO] Loaded {csv_name}  ({len(df)} rows)")
    print(f"[INFO] To check grid completeness run: python3 analyze_trnl_grid_coverage.py")

    prefix = f"{args.method}_{args.basis}"
    contour_fixed_y(
        df,
        y_values=(2.0, 2.25, 2.5, 2.75),
        value_col=args.value,
        levels=30,
        save=args.save,
        prefix=prefix,
    )


if __name__ == "__main__":
    main()
