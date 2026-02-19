"""
analyze_rot_coverage.py

Purpose
-------
Analyze rotational data coverage across methods (SCF, MP2, CCSD_T).
Shows which orientations have the most R values and identifies datasets
complete enough for plotting.

Absorbs the functionality of:
  - verification_of_rot_orientations_vs_R.py  (R-count per orientation)
  - print_top_orientations_bsse_vs_R.py        (BSSE(R) tables per orientation)

Usage
-----
  # Summary table for all methods
  python3 analyze_rot_coverage.py

  # Top-10 orientations per method
  python3 analyze_rot_coverage.py --detailed

  # BSSE(R) value tables for top 5 orientations (SCF only)
  python3 analyze_rot_coverage.py --method SCF --print-tables

  # BSSE(R) tables for top 7, looser R requirement
  python3 analyze_rot_coverage.py --method MP2 --print-tables --top 7 --minR 3

  # All options together
  python3 analyze_rot_coverage.py --method SCF --detailed --print-tables --top 5
"""

import argparse
import pandas as pd
from pathlib import Path


METHODS = ["SCF", "MP2", "CCSD_T"]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_method(method: str, basis: str) -> pd.DataFrame | None:
    p = Path(f"rot_{method}_{basis}.csv")
    if not p.exists():
        return None
    return pd.read_csv(p)


# ---------------------------------------------------------------------------
# Coverage analysis
# ---------------------------------------------------------------------------

def analyze_coverage(df: pd.DataFrame, minR: int) -> dict:
    coverage = (
        df.groupby(["alpha", "beta", "gamma"])["R"]
          .nunique()
          .sort_values(ascending=False)
    )
    eligible = coverage[coverage >= minR]
    return {
        "total_points":   len(df),
        "n_orientations": len(coverage),
        "n_eligible":     len(eligible),
        "max_R":          int(coverage.max()),
        "min_R":          int(coverage.min()),
        "mean_R":         round(coverage.mean(), 1),
        "coverage":       coverage,   # full Series, used by other functions
        "eligible":       eligible,
    }


# ---------------------------------------------------------------------------
# Printers
# ---------------------------------------------------------------------------

def print_summary(results: dict, minR: int) -> None:
    W = 72
    print("\n" + "=" * W)
    print(f"ROTATIONAL DATA COVERAGE SUMMARY  (minR = {minR})")
    print("=" * W)
    print(f"{'Method':<10} {'Exists':<8} {'Points':<10} {'Orient.':<10}"
          f"{'Eligible':<10} {'Max R':<8} {'Mean R':<8}")
    print("-" * W)

    for method, info in results.items():
        if info is None:
            print(f"{method:<10} {'No':<8} {'---':<10} {'---':<10}"
                  f"{'---':<10} {'---':<8} {'---':<8}")
        else:
            s = info["stats"]
            print(f"{method:<10} {'Yes':<8} {s['total_points']:<10}"
                  f"{s['n_orientations']:<10} {s['n_eligible']:<10}"
                  f"{s['max_R']:<8} {s['mean_R']:<8}")

    print("=" * W)
    print(f"\nLegend:")
    print(f"  Points:   Total (R, α, β, γ) data points")
    print(f"  Orient.:  Unique (α, β, γ) combinations")
    print(f"  Eligible: Orientations with >= {minR} distinct R values")
    print(f"  Max R:    Most R values held by any single orientation")
    print(f"  Mean R:   Average R values per orientation")


def print_detailed(results: dict) -> None:
    for method, info in results.items():
        if info is None:
            continue
        top = info["stats"]["coverage"].head(10)
        print(f"\n{'=' * 72}")
        print(f"{method}  —  TOP 10 ORIENTATIONS BY R COVERAGE")
        print(f"{'=' * 72}")
        print(f"{'Alpha':>8} {'Beta':>8} {'Gamma':>8}  {'# R values':>10}")
        print("-" * 40)
        for (a, b, g), nR in top.items():
            print(f"{a:>8.1f} {b:>8.1f} {g:>8.1f}  {int(nR):>10}")


def print_bsse_tables(results: dict, top: int, value_col: str) -> None:
    """Print BSSE(R) value tables for the top N eligible orientations."""
    for method, info in results.items():
        if info is None:
            continue

        df       = info["df"]
        eligible = info["stats"]["eligible"]

        if eligible.empty:
            print(f"\n[{method}] No eligible orientations to tabulate.")
            continue

        picked = eligible.head(top)

        print(f"\n{'=' * 72}")
        print(f"{method}  —  {value_col.upper()}(R) TABLES  "
              f"(top {len(picked)} orientations by R coverage)")
        print(f"{'=' * 72}")

        # Header ranking block (mirrors print_top_orientations_bsse_vs_R.py)
        print("\n[Top orientations by distinct R count]")
        print(picked.to_string())
        print()

        for (a, b, g), nR in picked.items():
            sel = (
                df[(df["alpha"] == a) & (df["beta"] == b) & (df["gamma"] == g)]
                .sort_values("R")
            )
            print("=" * 60)
            print(f"orientation (α,β,γ) = ({a:g}, {b:g}, {g:g})"
                  f"  |  distinct R = {int(nR)}")
            print(
                sel[["alpha", "beta", "gamma", "R", value_col]]
                .to_string(index=False)
            )
            print()


def print_recommendations(results: dict, minR: int) -> None:
    print(f"\n{'RECOMMENDATIONS':─^72}")
    for method, info in results.items():
        if info is None:
            print(f"  {method}: CSV not found — generate data first.")
            continue
        s = info["stats"]
        if s["n_eligible"] == 0:
            print(f"  {method}: No orientations with >= {minR} R values. "
                  f"Lower --minR or generate more R points.")
        elif s["n_eligible"] < 5:
            print(f"  {method}: Only {s['n_eligible']} eligible orientations — "
                  f"consider generating more R points.")
        else:
            print(f"  {method}: Good coverage — {s['n_eligible']} orientations "
                  f"with >= {minR} R values, ready to plot.")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyse rotational data coverage across methods."
    )
    parser.add_argument(
        "--method",
        choices=METHODS + ["all"],
        default="all",
        help="Method(s) to analyse (default: all).",
    )
    parser.add_argument(
        "--basis",
        default="aug-cc-pvdz",
        help="Basis set suffix used in CSV filenames.",
    )
    parser.add_argument(
        "--minR",
        type=int,
        default=5,
        help="Minimum distinct R values for an orientation to be 'eligible'.",
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Show top-10 orientations by R coverage for each method.",
    )
    parser.add_argument(
        "--print-tables",
        action="store_true",
        help="Print BSSE(R) value tables for the top --top orientations. "
             "Best used with --method SCF|MP2|CCSD_T.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="Number of orientations shown in --print-tables (default: 5).",
    )
    parser.add_argument(
        "--value",
        choices=["bsse", "no_vmfc", "vmfc"],
        default="bsse",
        help="Value column shown in --print-tables (default: bsse).",
    )
    args = parser.parse_args()

    methods = METHODS if args.method == "all" else [args.method]

    # Load data and compute coverage for each requested method
    results = {}
    for method in methods:
        df = load_method(method, args.basis)
        if df is None:
            print(f"[WARN] rot_{method}_{args.basis}.csv not found — skipping.")
            results[method] = None
        else:
            results[method] = {
                "df":    df,
                "stats": analyze_coverage(df, args.minR),
            }

    # Always: summary table
    print_summary(results, args.minR)

    # Optional: top-10 orientation lists
    if args.detailed:
        print_detailed(results)

    # Optional: BSSE(R) value tables
    if args.print_tables:
        print_bsse_tables(results, args.top, args.value)

    # Always: recommendations
    print_recommendations(results, args.minR)


if __name__ == "__main__":
    main()
