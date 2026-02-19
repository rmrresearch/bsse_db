"""
analyze_trnl_grid_coverage.py

Purpose
-------
Analyze the completeness of the translational (x, y, z) grid without
producing any plots. Separates the analysis concern from plot_trnl_contours.py.

Since all methods (SCF, MP2, CCSD_T) share the same (x,y,z) grid, the
grid analysis is run once on whichever CSV is available first. A warning
is printed if a requested method's file is missing.

Usage
-----
  # Analyse all available methods (uses first found for grid)
  python3 analyze_trnl_grid_coverage.py

  # Filter to a single method
  python3 analyze_trnl_grid_coverage.py --method MP2

  # Different basis set
  python3 analyze_trnl_grid_coverage.py --basis cc-pvtz
"""

import argparse
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

METHODS = ["SCF", "MP2", "CCSD_T"]


def load_first_available(methods: list[str], basis: str) -> tuple[pd.DataFrame, str]:
    """Return (df, method_name) for the first CSV that exists."""
    for m in methods:
        p = Path(f"trnl_{m}_{basis}.csv")
        if p.exists():
            return pd.read_csv(p), m
    raise FileNotFoundError(
        f"No translational CSV found for methods {methods} with basis '{basis}'.\n"
        f"Expected files like: trnl_SCF_{basis}.csv"
    )


def check_missing_files(methods: list[str], basis: str) -> None:
    """Warn about any method CSVs that are absent."""
    missing = [m for m in methods if not Path(f"trnl_{m}_{basis}.csv").exists()]
    if missing:
        for m in missing:
            print(f"[WARN] trnl_{m}_{basis}.csv not found — skipping.")


# ---------------------------------------------------------------------------
# Grid summary tables (extracted from plot_trnl_contours.py)
# ---------------------------------------------------------------------------

def _summary_fixed_x(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("x")[["y", "z"]]
          .apply(lambda g: pd.Series({
              "n_points":  len(g),
              "n_y":       g["y"].nunique(),
              "n_z":       g["z"].nunique(),
              "expected":  g["y"].nunique() * g["z"].nunique(),
              "missing":   g["y"].nunique() * g["z"].nunique() - len(g),
          }))
          .sort_values(["missing", "expected"], ascending=[True, False])
    )


def _summary_fixed_z(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("z")[["x", "y"]]
          .apply(lambda g: pd.Series({
              "n_points":  len(g),
              "n_x":       g["x"].nunique(),
              "n_y":       g["y"].nunique(),
              "expected":  g["x"].nunique() * g["y"].nunique(),
              "missing":   g["x"].nunique() * g["y"].nunique() - len(g),
          }))
          .sort_values(["missing", "expected"], ascending=[True, False])
    )


def _summary_fixed_y(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("y")[["x", "z"]]
          .apply(lambda g: pd.Series({
              "n_points":  len(g),
              "n_x":       g["x"].nunique(),
              "n_z":       g["z"].nunique(),
              "expected":  g["x"].nunique() * g["z"].nunique(),
              "missing":   g["x"].nunique() * g["z"].nunique() - len(g),
          }))
          .sort_values(["missing", "expected"], ascending=[True, False])
    )


# ---------------------------------------------------------------------------
# Per-axis recommendations
# ---------------------------------------------------------------------------

def _recommend(axis: str, fixed_val: float, missing: int, expected: int) -> str:
    """Return a one-line status string for a single slice."""
    if missing == 0:
        return f"  ✓  {axis}={fixed_val:g}  — complete ({expected} pts)  →  ready to plot"
    pct = 100 * missing / expected
    return (
        f"  ✗  {axis}={fixed_val:g}  — {missing}/{expected} pts missing ({pct:.0f}%)  "
        f"→  incomplete"
    )


def print_recommendations(
    sx: pd.DataFrame,
    sz: pd.DataFrame,
    sy: pd.DataFrame,
) -> None:
    """Print per-axis completeness recommendations."""

    width = 78
    print("\n" + "=" * width)
    print("RECOMMENDATIONS BY AXIS")
    print("=" * width)

    for label, summary, ax in [
        ("Fixed x  →  contour over (y, z)", sx, "x"),
        ("Fixed z  →  contour over (x, y)", sz, "z"),
        ("Fixed y  →  contour over (x, z)", sy, "y"),
    ]:
        complete   = summary[summary["missing"] == 0]
        incomplete = summary[summary["missing"] >  0]

        print(f"\n── {label} ──")
        print(f"   {len(complete)} complete / {len(summary)} total slices")

        if not complete.empty:
            print("   Fully plottable:")
            for val, row in complete.iterrows():
                print(_recommend(ax, val, 0, int(row["expected"])))

        if not incomplete.empty:
            print("   Partially complete (need more data to plot):")
            for val, row in incomplete.iterrows():
                print(_recommend(ax, val, int(row["missing"]), int(row["expected"])))

    print("\n" + "=" * width)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyse translational grid completeness (no plotting)."
    )
    parser.add_argument(
        "--method",
        choices=METHODS + ["all"],
        default="all",
        help="Method to report (default: all). Grid analysis runs on first "
             "available CSV since all methods share the same (x,y,z) grid.",
    )
    parser.add_argument(
        "--basis",
        default="aug-cc-pvdz",
        help="Basis set suffix used in CSV filenames (default: aug-cc-pvdz).",
    )
    args = parser.parse_args()

    methods = METHODS if args.method == "all" else [args.method]

    # Warn about absent files
    check_missing_files(methods, args.basis)

    # Load grid from first available CSV
    df, source_method = load_first_available(methods, args.basis)
    print(f"\n[INFO] Grid loaded from trnl_{source_method}_{args.basis}.csv")
    print(f"       Total data points : {len(df)}")
    print(f"       Unique x values   : {sorted(df['x'].unique())}")
    print(f"       Unique y values   : {sorted(df['y'].unique())}")
    print(f"       Unique z values   : {sorted(df['z'].unique())}")

    # Build the three summary tables
    sx = _summary_fixed_x(df)
    sz = _summary_fixed_z(df)
    sy = _summary_fixed_y(df)

    # Print the same detailed tables as plot_trnl_contours.py did
    print("\n=== Fixed x  →  (y, z) grid completeness ===")
    print(sx.to_string())

    print("\n=== Fixed z  →  (x, y) grid completeness ===")
    print(sz.to_string())

    print("\n=== Fixed y  →  (x, z) grid completeness ===")
    print(sy.to_string())

    # Concise recommendations
    print_recommendations(sx, sz, sy)


if __name__ == "__main__":
    main()
