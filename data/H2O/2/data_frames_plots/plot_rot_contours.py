# plot_rot_contours.py

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def summaries_rot(df: pd.DataFrame) -> None:
    """Print grid-completeness summaries for several 2D slice families."""

    # Fix (R, gamma) -> vary (alpha, beta)
    summary_Rg = (
        df.groupby(["R", "gamma"])[["alpha", "beta"]]
          .apply(lambda g: pd.Series({
              "n_points": len(g),
              "n_alpha": g["alpha"].nunique(),
              "n_beta": g["beta"].nunique(),
              "expected": g["alpha"].nunique() * g["beta"].nunique(),
              "missing": g["alpha"].nunique() * g["beta"].nunique() - len(g),
          }))
          .sort_values(["missing", "expected"], ascending=[True, False])
    )

    # Fix (alpha, beta) -> vary (R, gamma)
    summary_ab = (
        df.groupby(["alpha", "beta"])[["R", "gamma"]]
          .apply(lambda g: pd.Series({
              "n_points": len(g),
              "n_R": g["R"].nunique(),
              "n_gamma": g["gamma"].nunique(),
              "expected": g["R"].nunique() * g["gamma"].nunique(),
              "missing": g["R"].nunique() * g["gamma"].nunique() - len(g),
          }))
          .sort_values(["missing", "expected"], ascending=[True, False])
    )

    # Fix (alpha, gamma) -> vary (beta, R)
    summary_ag = (
        df.groupby(["alpha", "gamma"])[["beta", "R"]]
          .apply(lambda g: pd.Series({
              "n_points": len(g),
              "n_beta": g["beta"].nunique(),
              "n_R": g["R"].nunique(),
              "expected": g["beta"].nunique() * g["R"].nunique(),
              "missing": g["beta"].nunique() * g["R"].nunique() - len(g),
          }))
          .sort_values(["missing", "expected"], ascending=[True, False])
    )

    print("\n=== Fix (R, gamma) -> alpha–beta grid completeness ===")
    print(summary_Rg)

    print("\n=== Fix (alpha, beta) -> R–gamma grid completeness ===")
    print(summary_ab)

    print("\n=== Fix (alpha, gamma) -> beta–R grid completeness ===")
    print(summary_ag)


def plot_alpha_beta_2x4(
    df: pd.DataFrame,
    value_col: str = "bsse",
    levels: int = 20,
    save: bool = False,
    prefix: str = "",
) -> None:
    """Plot only slices where fixing (R,gamma) yields a complete 2x4 alpha-beta grid."""

    # Build the (R,gamma) summary (same as in summaries)
    summary_Rg = (
        df.groupby(["R", "gamma"])[["alpha", "beta"]]
          .apply(lambda g: pd.Series({
              "n_points": len(g),
              "n_alpha": g["alpha"].nunique(),
              "n_beta": g["beta"].nunique(),
              "expected": g["alpha"].nunique() * g["beta"].nunique(),
              "missing": g["alpha"].nunique() * g["beta"].nunique() - len(g),
          }))
    )

    # Keep only the "true 2D" 2x4 grids
    good = summary_Rg.query("missing == 0 and n_alpha == 2 and n_beta == 4").sort_index()

    if good.empty:
        print("[INFO] No complete 2x4 (alpha,beta) grids found for any (R,gamma).")
        return

    print("\n=== Will plot these (R, gamma) slices (2x4 alpha–beta grids) ===")
    print(good)

    for (R0, g0), stats in good.iterrows():
        sl = df[(df["R"] == R0) & (df["gamma"] == g0)].copy()

        # Pivot to rectangular grid: rows=alpha, cols=beta
        grid = sl.pivot(index="alpha", columns="beta", values=value_col)

        # Extra safety check
        missing = int(grid.isna().sum().sum())
        if missing != 0:
            print(f"[WARN] Skipping R={R0}, gamma={g0}: pivot has {missing} NaNs.")
            continue

        A, B = np.meshgrid(grid.columns.to_numpy(), grid.index.to_numpy())
        V = grid.to_numpy()

        plt.figure()
        plt.contourf(A, B, V, levels=levels)
        plt.colorbar(label=f"{value_col} (Hartree)")
        plt.xlabel("beta (deg)")
        plt.ylabel("alpha (deg)")
        plt.title(f"{value_col} at R={R0}, gamma={g0} ({prefix})")
        plt.tight_layout()

        if save:
            fname = f"rot_ab_R{R0:g}_g{g0:g}_{value_col}_{prefix}.png"
            plt.savefig(fname, dpi=200)
            print("saved:", fname)

        plt.show()


def load_method_df(method: str, basis: str = "aug-cc-pvdz") -> pd.DataFrame:
    csv_name = f"rot_{method}_{basis}.csv"
    return pd.read_csv(csv_name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=["SCF", "MP2", "CCSD_T"], default="SCF")
    parser.add_argument("--basis", default="aug-cc-pvdz")
    parser.add_argument("--value", choices=["bsse", "no_vmfc", "vmfc"], default="bsse")
    parser.add_argument("--save", action="store_true")
    parser.add_argument("--levels", type=int, default=20)
    args = parser.parse_args()

    df = load_method_df(args.method, args.basis)
    prefix = f"{args.method}_{args.basis}"

    # 1) Print summaries (for explaining to PI)
    summaries_rot(df)

    # 2) Plot only the valid 2D 2x4 slices
    plot_alpha_beta_2x4(
        df,
        value_col=args.value,
        levels=args.levels,
        save=args.save,
        prefix=prefix,
    )


if __name__ == "__main__":
    main()
