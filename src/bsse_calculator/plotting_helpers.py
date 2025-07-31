import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata, UnivariateSpline
from scipy.optimize import curve_fit
from itertools import combinations
from collections import defaultdict
from typing import Optional
import json
import re


def make_bsse_contour_plots(bsse_data_nwchem, monomer_name, title=""):
    num_plots = bsse_data_nwchem["y"].nunique()
    if num_plots > 8:
        num_plots = 8
    fig, axes = plt.subplots(nrows=num_plots // 2, ncols=2, figsize=(20, 20))
    axes = axes.flatten()
    cmap = cm.get_cmap("viridis")
    vmin = bsse_data_nwchem["BSSE_A"].min()
    vmax = bsse_data_nwchem["BSSE_A"].max()
    normalizer = Normalize(vmin=vmin, vmax=vmax)

    for ax, fixed_y_value in zip(
        axes, (np.arange(0.25, bsse_data_nwchem["y"].max(), 0.25))
    ):
        data = bsse_data_nwchem[bsse_data_nwchem["y"] == fixed_y_value][
            ["x", "z", "BSSE_A"]
        ]
        xs = data["x"].values
        zs = data["z"].values
        bsse = data["BSSE_A"].values

        xi = np.linspace(min(xs), max(xs), 50)
        zi = np.linspace(min(zs), max(zs), 50)
        Xi, Zi = np.meshgrid(xi, zi)
        bssei = griddata((xs, zs), bsse, (Xi, Zi), method="cubic")

        contour = ax.contourf(
            Xi, Zi, bssei, levels=8, cmap="viridis", vmin=vmin, vmax=vmax
        )
        ax.set_title(f"Contour Plot for y = {fixed_y_value} {title}")
        ax.set_xlabel(f"X Coordinate of Second {monomer_name} (Å)")
        ax.set_ylabel(f"Z Coordinate of Second {monomer_name} (Å)")
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=cmap), ax=axes.ravel().tolist()
    )
    cbar.set_label("BSSE (kcal/mol)", labelpad=20)
    plt.show()

    return fig


def fit_gaussian(bsse_data_nwchem, num_gaussians=2):
    def gauss(xyz, *p):
        output = 0
        for i in range(num_gaussians):
            r = (
                (xyz[:, 0] - p[3 * i]) ** 2
                + (xyz[:, 1] - p[3 * i + 1]) ** 2
                + (xyz[:, 2] - p[3 * i + 2]) ** 2
            )
            if p[3 * i + 4] < 0:
                width = 0
            else:
                width = p[3 * i + 4]
            output += p[3 * i + 3] * np.exp(-r * width)
        return output

    p0 = []
    for _ in range(num_gaussians):
        p0.append(0)
        p0.append(0)
        p0.append(0)
        p0.append(1)
        p0.append(1)
    coords = bsse_data_nwchem[["x", "y", "z"]].values
    bsses = bsse_data_nwchem["BSSE_A"].values
    coeff, var_matrix = curve_fit(gauss, coords, bsses, p0=p0, maxfev=10000)

    def fit_gauss(xyz):
        return gauss(xyz, *coeff)

    return fit_gauss


def make_contour_plots_of_fitted_gaussian(bsse_data_nwchem, fit_gauss):
    num_plots = bsse_data_nwchem["y"].nunique()
    fig, axes = plt.subplots(nrows=num_plots // 2, ncols=2, figsize=(25, 25))
    cmap = cm.get_cmap("viridis")
    vmin = bsse_data_nwchem["BSSE_A"].min()
    vmax = bsse_data_nwchem["BSSE_A"].max()
    normalizer = Normalize(vmin=vmin, vmax=vmax)

    for i, fixed_y_value in enumerate(
        np.arange(0.25, bsse_data_nwchem["y"].max(), 0.25)
    ):
        ax = axes[i // 2, i % 2]
        x = np.linspace(0, 2, 200)
        z = np.linspace(0, 2, 200)
        X, Z = np.meshgrid(x, z)
        grid_points = np.column_stack(
            [X.ravel(), fixed_y_value * np.ones_like(X.ravel()), Z.ravel()]
        )
        Z_values = fit_gauss(grid_points)
        Z_reshaped = Z_values.reshape(X.shape)
        mask = X**2 + fixed_y_value**2 + Z**2 < 1.5
        Z_masked = np.where(mask, np.nan, Z_reshaped)
        contour = ax.contourf(
            X, Z, Z_masked, levels=50, cmap="viridis", vmin=vmin, vmax=vmax
        )
        ax.set_title(f"Estimated BSSE for y = {fixed_y_value}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")

    cbar = fig.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=cmap), ax=axes.ravel().tolist()
    )
    cbar.set_label("Estimated BSSE (kcal/mol)", rotation=270, labelpad=20)
    plt.show()


def predicted_v_actual_plot(bsse_data_nwchem, fit_gauss, title):
    predicted_bsses = fit_gauss(bsse_data_nwchem[["x", "y", "z"]].values)
    actual_bsse = bsse_data_nwchem["BSSE_A"].values
    bsse_data_nwchem["predicted_bsse"] = predicted_bsses
    plt.figure(figsize=(5, 5))
    plt.scatter(predicted_bsses, actual_bsse, alpha=0.6, s=60)

    # Add y=x line (perfect prediction line)
    min_val = min(predicted_bsses.min(), actual_bsse.min())
    max_val = max(predicted_bsses.max(), actual_bsse.max())
    plt.plot(
        [min_val, max_val],
        [min_val, max_val],
        "r--",
        linewidth=2,
        label="Perfect Prediction (y=x)",
    )

    plt.xlabel("Predicted BSSE (kcal/mol)")
    plt.ylabel("Actual BSSE (kcal/mol)")
    plt.title(f"Predicted vs Actual BSSE Values {title}")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # Make the plot square for better comparison
    plt.axis("equal")
    plt.tight_layout()
    plt.show()
    #   print("RMSE:", np.sqrt(np.mean((predicted_bsses - actual_bsse) ** 2)))
    #   print("MAE:", np.mean(np.abs(predicted_bsses - actual_bsse)))
    #   print(r"% Error:", np.mean(np.abs((predicted_bsses - actual_bsse) / actual_bsse)) * 100)
    return {
        "RMSE": np.sqrt(np.mean((predicted_bsses - actual_bsse) ** 2)),
        "MAE": np.mean(np.abs(predicted_bsses - actual_bsse)),
        "% Error": np.mean(np.abs((predicted_bsses - actual_bsse) / actual_bsse)) * 100,
    }


def make_error_contour_plots(bsse_data_nwchem, fit_gauss):
    predicted_bsses = fit_gauss(bsse_data_nwchem[["x", "y", "z"]].values)
    actual_bsse = bsse_data_nwchem["BSSE_A"].values
    num_plots = bsse_data_nwchem["y"].nunique()
    fig, axes = plt.subplots(nrows=num_plots // 2, ncols=2, figsize=(25, 25))
    error = predicted_bsses - actual_bsse
    cmap = cm.get_cmap("viridis")
    vmin = error.min()
    vmax = error.max()
    normalizer = Normalize(vmin=vmin, vmax=vmax)
    for i, fixed_y_value in enumerate(
        np.arange(0.25, bsse_data_nwchem["y"].max(), 0.25)
    ):
        ax = axes[i // 2, i % 2]

        mask = bsse_data_nwchem["y"] == fixed_y_value
        x_data = bsse_data_nwchem[mask]["x"].values
        z_data = bsse_data_nwchem[mask]["z"].values
        error_data = error[mask]

        x_grid = np.linspace(x_data.min(), x_data.max(), 100)
        z_grid = np.linspace(z_data.min(), z_data.max(), 100)
        X_grid, Z_grid = np.meshgrid(x_grid, z_grid)

        error_grid = griddata(
            (x_data, z_data),
            error_data,
            (X_grid, Z_grid),
            method="cubic",
            fill_value=np.nan,
        )

        contour = ax.contourf(
            X_grid, Z_grid, error_grid, levels=10, cmap=cmap, vmin=vmin, vmax=vmax
        )

        ax.scatter(
            x_data,
            z_data,
            c=error_data,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            s=30,
            edgecolors="white",
            linewidth=0.5,
            alpha=0.8,
        )

        ax.set_title(f"Absolute Error for y = {fixed_y_value}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")

    cbar = fig.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=cmap), ax=axes.ravel().tolist()
    )
    cbar.set_label("Predicted BSSE - Actual BSSE", rotation=270, labelpad=20)
    plt.show()


def prepare_bsse_data(file):
    bsse_data_nwchem_kcal = pd.read_csv(file)
    for column in bsse_data_nwchem_kcal.columns:
        if "energy" in column:
            bsse_data_nwchem_kcal[column] = bsse_data_nwchem_kcal[column].apply(
                lambda x: float(x.strip("()").split(",")[0])
            )
            bsse_data_nwchem_kcal[column] = bsse_data_nwchem_kcal[column].apply(
                lambda x: x * 627.5096
            )
    return bsse_data_nwchem_kcal


def plot_energy_types_by_distance(df, monomer, fit_spline=False):
    energy_types = [
        "Total Energy",
        "SCF Energy",
        "MP2 Correlation",
        "CCSD Correlation",
        "(T) Correlation",
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    delta_color = "blue"
    bsse_color = "red"

    all_lines = []
    all_labels = []

    for i, label in enumerate(energy_types):
        delta_col = [
            col
            for col in df.columns
            if label.title() in col.replace("_", " ").title() and "Delta" in col
        ][0]
        bsse_col = [
            col
            for col in df.columns
            if label.title() in col.replace("_", " ").title() and "BSSE" in col
        ][0]

        if i < len(axes):
            ax1 = axes[i]
            x = df["distance"].values
            y_delta = df[delta_col].values
            y_bsse = df[bsse_col].values

            line1 = ax1.scatter(x, y_delta, color=delta_color, alpha=0.7, s=50)
            ax1.set_xlabel(f"{monomer} Distance", fontsize=10)
            ax1.set_ylabel("Delta_EAB_EAB", fontsize=10)

            if fit_spline:
                delta_spline = UnivariateSpline(x, y_delta, s=0)
                ax1.plot(
                    np.sort(x),
                    delta_spline(np.sort(x)),
                    color=delta_color,
                    linestyle="--",
                )

            ax2 = ax1.twinx()
            line2 = ax2.scatter(x, y_bsse, color=bsse_color, alpha=0.7, s=50)
            ax2.set_ylabel("BSSE_A", fontsize=10)

            if fit_spline:
                bsse_spline = UnivariateSpline(x, y_bsse, s=0)
                ax2.plot(
                    np.sort(x),
                    bsse_spline(np.sort(x)),
                    color=bsse_color,
                    linestyle="--",
                )

            ax1.set_title(f"{label}", fontsize=11)
            ax1.grid(True, alpha=0.3)

            # Store one line from each for the shared legend
            if i == 0:
                all_lines.extend([line1, line2])
                all_labels.extend(["Delta_EAB_EAB", "BSSE_A"])

    # Turn off unused subplot if energy_types < axes
    if len(energy_types) < len(axes):
        for j in range(len(energy_types), len(axes)):
            axes[j].axis("off")
    fig.legend(
        all_lines,
        all_labels,
        loc="upper center",
        ncol=2,
        fontsize=12,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(
        f"Delta_EAB_EAB vs BSSE_A by Energy Component for {monomer}",
        fontsize=14,
        y=1.08,
    )
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)


def create_df(energies):
    energy_types = next(iter(next(iter(energies.values())).values())).keys()
    regrouped_energies = defaultdict(list)
    for side_length in energies.keys():
        for energy_type in energy_types:
            row = {"side_length": side_length}
            for calculation in energies[side_length].keys():
                calculation_new = calculation.replace("_", "").replace(")", ")_")
                row[calculation_new] = (
                    energies[side_length][calculation][energy_type] * 627.5096
                )
            regrouped_energies[energy_type].append(row)
    return {
        energy_type: pd.DataFrame(regrouped_energies[energy_type])
        for energy_type in energy_types
    }


def create_BSSEs_and_interaction_energies(total_energies):
    # get interaction for each dimer
    for monomers in combinations(["A", "B", "C"], 2):
        dimer_basis_name = (
            f"Delta_E({monomers[0]}{monomers[1]})_{monomers[0]}{monomers[1]}"
        )
        dimer = f"{monomers[0]}{monomers[1]}"
        total_energies[dimer_basis_name] = (
            total_energies[f"E({dimer})_{dimer}"]
            - total_energies[f"E({dimer})_{monomers[0]}"]
            - total_energies[f"E({dimer})_{monomers[1]}"]
        )
        # BSSE of adding third basis
        full_basis_name = f"Delta_E(ABC)_{monomers[0]}{monomers[1]}"
        total_energies[full_basis_name] = (
            total_energies[f"E(ABC)_{dimer}"]
            - total_energies[f"E(ABC)_{monomers[0]}"]
            - total_energies[f"E(ABC)_{monomers[1]}"]
        )
        total_energies[f"BSSE_({dimer})"] = (
            total_energies[dimer_basis_name] - total_energies[full_basis_name]
        )
    working_interaction_energy = total_energies["E(ABC)_ABC"].copy()
    for monomer in ["A", "B", "C"]:
        working_interaction_energy -= total_energies[f"E(ABC)_{monomer}"]
    for monomers in combinations(["A", "B", "C"], 2):
        working_interaction_energy += (
            -total_energies[f"Delta_E(ABC)_{monomers[0]}{monomers[1]}"]
            + total_energies[
                f"Delta_E({monomers[0]}{monomers[1]})_{monomers[0]}{monomers[1]}"
            ]
        )
    total_energies["Total Interaction Energy (BSSE free)"] = working_interaction_energy
    working_total_interaction_energy_wo_bsse_correction = total_energies[
        "E(ABC)_ABC"
    ].copy()
    for monomer in ["A", "B", "C"]:
        working_total_interaction_energy_wo_bsse_correction -= total_energies[
            f"E({monomer})_{monomer}"
        ]
    # E(A_B_C)_A_B_C - E(A)_A - E(B)_B - E(C)_C
    total_energies["Total Interaction Energy (w/o BSSE correction)"] = (
        working_total_interaction_energy_wo_bsse_correction
    )
    total_energies["Total BSSE"] = (
        total_energies["Total Interaction Energy (w/o BSSE correction)"]
        - total_energies["Total Interaction Energy (BSSE free)"]
    )
    return total_energies


def plot_interaction_and_bsse(total_energies):

    fig, ax1 = plt.subplots(figsize=(8, 6))

    sorted_indices = np.argsort(total_energies["side_length"].values.astype(float))
    x_sorted = np.array(total_energies["side_length"])[sorted_indices].astype(float)
    # x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), 300)

    y_delta_e = np.array(total_energies["Delta_E(ABC)_ABC"])[sorted_indices]
    # spline_delta_e = UnivariateSpline(x_sorted, y_delta_e, s=0)
    ax1.plot(x_sorted, y_delta_e, color="black", label="ΔE(ABC)", zorder=3, marker="o")
    # ax1.plot(x_smooth, spline_delta_e(x_smooth), color='black', alpha=0.7, linewidth=2)
    ax1.set_ylabel("ΔE(ABC)", color="black")
    ax1.tick_params(axis="y", labelcolor="black")

    ax2 = ax1.twinx()

    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    color_idx = 0

    for monomers in combinations(["B", "C"], 2):
        y_bsse = np.array(total_energies[f"BSSE_({monomers[0]}_{monomers[1]})"])[
            sorted_indices
        ]
        # spline_bsse = UnivariateSpline(x_sorted, y_bsse, s=0)
        ax2.plot(
            x_sorted,
            y_bsse,
            label=f"BSSE {monomers[0]}{monomers[1]}",
            color=colors[color_idx],
            marker="o",
        )
        # ax2.plot(x_smooth, spline_bsse(x_smooth), color=colors[color_idx], alpha=0.7, linewidth=2)
        color_idx += 1

    for monomer in ["C"]:
        y_bsse_mono = np.array(total_energies[f"BSSE_({monomer})_ABC"])[sorted_indices]
        # spline_bsse_mono = UnivariateSpline(x_sorted, y_bsse_mono, s=0)
        ax2.plot(
            x_sorted,
            y_bsse_mono,
            label=f"BSSE {monomer}",
            color=colors[color_idx],
            marker="o",
        )
        # ax2.plot(x_smooth, spline_bsse_mono(x_smooth), color=colors[color_idx], alpha=0.7, linewidth=2)
        color_idx += 1

    ax2.set_ylabel("BSSE terms")
    ax1.set_xlabel("Side Length")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(
        lines1 + lines2, labels1 + labels2, bbox_to_anchor=(1.15, 1), loc="upper left"
    )

    fig.tight_layout()
    plt.show()


def plot_similarity(ax, list1, list2, title="", x_label="", y_label=""):
    x = np.array(list1)
    y = np.array(list2)
    min_len = min(len(x), len(y))
    x = x[:min_len]
    y = y[:min_len]

    percent_error = np.mean(np.abs(x - y) / x) * 100
    min_val = min(np.min(x), np.min(y))
    max_val = max(np.max(x), np.max(y))
    range_val = max_val - min_val
    padding = range_val * 0.05

    ax.scatter(x, y, alpha=0.7, s=50, color="blue", label="Data points")
    ax.plot(
        [min_val, max_val],
        [min_val, max_val],
        "r--",
        linewidth=2,
        label="y = x (perfect agreement)",
    )
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(f"{title}\nPercent Error: {percent_error:.2f}%", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(min_val - padding, max_val + padding)
    ax.set_ylim(min_val - padding, max_val + padding)
    ax.set_aspect("equal")
    ax.legend()

    return percent_error


def plot_interaction_with_extra_basis(total_energies):
    # Create 3 subplots (for 3 combinations of monomers A, B, C)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax, monomers in zip(axes, combinations(["A", "B", "C"], 2)):
        plot_similarity(
            ax,
            total_energies[
                f"Delta_E({monomers[0]}{monomers[1]})_{monomers[0]}{monomers[1]}"
            ],
            total_energies[f"Delta_E({monomers[0]}{monomers[1]})_ABC"],
            f"ΔE({monomers[0]}{monomers[1]})_{monomers[0]}{monomers[1]} vs ΔE({monomers[0]}{monomers[1]})_ABC",
            f"ΔE({monomers[0]}{monomers[1]})_AB",
            f"ΔE({monomers[0]}{monomers[1]})_ABC",
        )

    plt.tight_layout()
    plt.show()


def calculate_total_interaction_energy_with_bsse(df):
    working_total_interaction_energy_with_bsse = df["E(ABCD)_ABCD"].copy()
    for monomer in ["A", "B", "C", "D"]:
        working_total_interaction_energy_with_bsse -= df[f"E({monomer})_{monomer}"]
    df["Total Interaction Energy (w/ BSSE)"] = (
        working_total_interaction_energy_with_bsse
    )


def calculate_dimer_interaction_energies(df):
    for monomer1, monomer2 in combinations(["A", "B", "C", "D"], 2):
        df[f"Delta_E({monomer1}{monomer2})_{monomer1}{monomer2}"] = (
            df[f"E({monomer1}{monomer2})_{monomer1}{monomer2}"]
            - df[f"E({monomer1}{monomer2})_{monomer1}"]
            - df[f"E({monomer1}{monomer2})_{monomer2}"]
        )
    for monomer1, monomer2, monomer3 in combinations(["A", "B", "C", "D"], 3):
        for dimer in combinations([monomer1, monomer2, monomer3], 2):
            df[f"Delta_E({monomer1}{monomer2}{monomer3})_{dimer[0]}{dimer[1]}"] = (
                df[f"E({monomer1}{monomer2}{monomer3})_{dimer[0]}{dimer[1]}"]
                - df[f"E({monomer1}{monomer2}{monomer3})_{dimer[0]}"]
                - df[f"E({monomer1}{monomer2}{monomer3})_{dimer[1]}"]
            )
    for dimer in combinations(["A", "B", "C", "D"], 2):
        df[f"Delta_E(ABCD)_{dimer[0]}{dimer[1]}"] = (
            df[f"E(ABCD)_{dimer[0]}{dimer[1]}"]
            - df[f"E(ABCD)_{dimer[0]}"]
            - df[f"E(ABCD)_{dimer[1]}"]
        )


def calculate_trimer_interaction_energies(df):
    for monomer1, monomer2, monomer3 in combinations(["A", "B", "C", "D"], 3):
        working_interaction_energy = df[
            f"E({monomer1}{monomer2}{monomer3})_{monomer1}{monomer2}{monomer3}"
        ].copy()
        for monomer in [monomer1, monomer2, monomer3]:
            working_interaction_energy -= df[
                f"E({monomer1}{monomer2}{monomer3})_{monomer}"
            ]
        for dimer in combinations([monomer1, monomer2, monomer3], 2):
            working_interaction_energy -= df[
                f"Delta_E({monomer1}{monomer2}{monomer3})_{dimer[0]}{dimer[1]}"
            ]
        df[
            f"Delta_E({monomer1}{monomer2}{monomer3})_{monomer1}{monomer2}{monomer3}"
        ] = working_interaction_energy
    for monomer1, monomer2, monomer3 in combinations(["A", "B", "C", "D"], 3):
        working_interaction_energy_tetramer_basis = df[
            f"E(ABCD)_{monomer1}{monomer2}{monomer3}"
        ].copy()
        for monomer in [monomer1, monomer2, monomer3]:
            working_interaction_energy_tetramer_basis -= df[f"E(ABCD)_{monomer}"]
        for dimer in combinations([monomer1, monomer2, monomer3], 2):
            working_interaction_energy_tetramer_basis -= df[
                f"Delta_E(ABCD)_{dimer[0]}{dimer[1]}"
            ]
        df[f"Delta_E(ABCD)_{monomer1}{monomer2}{monomer3}"] = (
            working_interaction_energy_tetramer_basis
        )


def calculate_total_interaction_energy_with_correction(df):
    working_total_interaction_energy_with_correction = df["E(ABCD)_ABCD"].copy()
    for monomer in ["A", "B", "C", "D"]:
        working_total_interaction_energy_with_correction -= df[f"E(ABCD)_{monomer}"]
    for dimer in combinations(["A", "B", "C", "D"], 2):
        working_total_interaction_energy_with_correction += (
            -df[f"Delta_E(ABCD)_{dimer[0]}{dimer[1]}"]
            + df[f"Delta_E({dimer[0]}{dimer[1]})_{dimer[0]}{dimer[1]}"]
        )
    for trimer in combinations(["A", "B", "C", "D"], 3):
        working_total_interaction_energy_with_correction += (
            -df[f"Delta_E(ABCD)_{trimer[0]}{trimer[1]}{trimer[2]}"]
            + df[
                f"Delta_E({trimer[0]}{trimer[1]}{trimer[2]})_{trimer[0]}{trimer[1]}{trimer[2]}"
            ]
        )
    df["Total Interaction Energy (BSSE free)"] = (
        working_total_interaction_energy_with_correction
    )


def calculate_total_bsse(df):
    df["Total BSSE"] = (
        df["Total Interaction Energy (w/ BSSE)"]
        - df["Total Interaction Energy (BSSE free)"]
    )


def get_dimer_bsse(df):
    for dimer in combinations(["A", "B", "C", "D"], 2):
        df[f"BSSE_Dimer({dimer[0]}{dimer[1]})"] = (
            df[f"Delta_E({dimer[0]}{dimer[1]})_{dimer[0]}{dimer[1]}"]
            - df[f"Delta_E(ABCD)_{dimer[0]}{dimer[1]}"]
        )


def get_trimer_bsse(df):
    for trimer in combinations(["A", "B", "C", "D"], 3):
        df[f"BSSE_Trimer({trimer[0]}{trimer[1]}{trimer[2]})"] = (
            df[
                f"Delta_E({trimer[0]}{trimer[1]}{trimer[2]})_{trimer[0]}{trimer[1]}{trimer[2]}"
            ]
            - df[f"Delta_E(ABCD)_{trimer[0]}{trimer[1]}{trimer[2]}"]
        )


def get_monomer_bsse(df):
    for monomer in ["A", "B", "C", "D"]:
        df[f"BSSE_Monomer({monomer})"] = (
            df[f"E({monomer})_{monomer}"] - df[f"E(ABCD)_{monomer}"]
        )


def group_bsse(df):
    total_dimer_bsse = pd.Series(0, index=df.index)
    for column in df.columns:
        if "BSSE_Dimer" in column:
            total_dimer_bsse += df[column]
    df["BSSE Dimer Sum"] = total_dimer_bsse

    total_trimer_bsse = pd.Series(0, index=df.index)
    for column in df.columns:
        if "BSSE_Trimer" in column:
            total_trimer_bsse += df[column]
    df["BSSE Trimer Sum"] = total_trimer_bsse

    total_Monomer_bsse = pd.Series(0, index=df.index)
    for column in df.columns:
        if "BSSE_Monomer" in column:
            total_Monomer_bsse += df[column]
    df["BSSE Monomer Sum"] = total_Monomer_bsse


def run_all_calculations(df):
    calculate_total_interaction_energy_with_bsse(df)
    calculate_dimer_interaction_energies(df)
    calculate_trimer_interaction_energies(df)
    calculate_total_interaction_energy_with_correction(df)
    calculate_total_bsse(df)
    get_dimer_bsse(df)
    get_trimer_bsse(df)
    get_monomer_bsse(df)
    group_bsse(df)
    return df


import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional


def plot_dual_axis(
    series1: pd.Series,
    series2: pd.Series,
    x_series: pd.Series,
    type: str = "line",
    *,
    ax1: Optional[plt.Axes] = None,
    title: str = "Dual Axis Plot",
    x_label: str = "X-axis",
    y1_label: str = "Left Y-axis",
    y2_label: str = "Right Y-axis",
    label1: str = "Series 1",
    label2: str = "Series 2",
    show_legend: bool = True,
):
    fig, ax1 = plt.subplots(figsize=(10, 6)) if ax1 is None else (None, ax1)

    if type == "line":
        ax1.plot(
            x_series.values, series1.values, color="tab:blue", label=label1, marker="o"
        )
        ax2 = ax1.twinx()
        ax2.plot(
            x_series.values, series2.values, color="tab:red", label=label2, marker="o"
        )
    elif type == "scatter":
        ax1.scatter(x_series.values, series1.values, color="tab:blue", label=label1)
        ax2 = ax1.twinx()
        ax2.scatter(x_series.values, series2.values, color="tab:red", label=label2)
    else:  # bar plot
        width = 0.4
        x = x_series.values
        x1 = [i - width / 2 for i in range(len(x))]
        x2 = [i + width / 2 for i in range(len(x))]

        ax1.bar(
            x1,
            series1.values,
            width=width,
            color="tab:blue",
            label=label1,
            align="center",
        )
        ax2 = ax1.twinx()
        ax2.bar(
            x2,
            series2.values,
            width=width,
            color="tab:red",
            label=label2,
            align="center",
        )

        ax1.set_xticks(range(len(x)))
        ax1.set_xticklabels(x)

    ax1.set_ylabel(y1_label)
    ax2.set_ylabel(y2_label)
    ax1.set_xlabel(x_label)
    ax1.set_title(title)
    ax1.tick_params(axis="y")
    ax2.tick_params(axis="y")

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    if show_legend:
        ax1.legend(lines1 + lines2, labels1 + labels2)
    return lines1 + lines2, labels1 + labels2


def create_df_euler_angles(filename):
    energies = json.loads(open(filename, "r").read())
    energies_df = []
    for phi in energies:
        for theta in energies[phi]:
            for psi in energies[phi][theta]:
                for translation in energies[phi][theta][psi]:
                    record = {
                        "phi": phi,
                        "theta": theta,
                        "psi": psi,
                        "translation": translation,
                        **{
                            k: energy * 627.5096
                            for k, energy in energies[phi][theta][psi][
                                translation
                            ].items()
                        },
                    }
                    energies_df.append(record)
    energies_df = pd.DataFrame(energies_df)
    for angle in ["phi", "theta", "psi"]:
        energies_df[angle] = energies_df[angle].astype(float) * 180 / np.pi
    grouped_df = energies_df.groupby("translation").agg(["mean", "min", "max"])
    return grouped_df


def plot_energy_vs_bsse_euler_angles(grouped_df):
    energy_types = [
        column.replace("_BSSE_A", "")
        for column in grouped_df.columns.get_level_values(0).unique()
        if "BSSE_A" in column
    ]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    # To collect legend entries
    handles_labels = []
    for ax, energy_type in zip(axes, energy_types):
        formated_title = energy_type.replace("_", " ").title()

        feature1 = grouped_df[f"{energy_type}_Delta_E_AB_AB"]
        feature2 = grouped_df[f"{energy_type}_BSSE_A"]

        yerr1 = [feature1["mean"] - feature1["min"], feature1["max"] - feature1["mean"]]
        yerr2 = [feature2["mean"] - feature2["min"], feature2["max"] - feature2["mean"]]

        h1 = ax.errorbar(
            x=feature1.index,
            y=feature1["mean"],
            yerr=yerr1,
            fmt="o",
            capsize=5,
            color="black",
            ecolor="gray",
            label="Total_energy_BSSE_A",
        )
        ax.set_ylabel("Total Interaction Energy (kcal/mol)")
        ax.set_title(formated_title)

        ax2 = ax.twinx()
        h2 = ax2.errorbar(
            x=feature2.index,
            y=feature2["mean"],
            yerr=yerr2,
            fmt="s",
            capsize=5,
            color="blue",
            ecolor="lightblue",
            label="BSSE_A",
        )
        ax2.set_ylabel("BSSE (kcal/mol)")

    plt.tight_layout()
    bbox = axes[-1].get_position()
    fig.delaxes(axes[-1])  # Remove empty subplot if needed

    # Add a joint legend to the figure
    fig.legend(
        handles=[h1, h2],
        labels=["Delta_E_AB_AB", "BSSE_A"],
        loc="center",
        bbox_to_anchor=(bbox.x0 + bbox.width / 2, bbox.y0 + bbox.height / 2),
        frameon=False,
    )


def plot_bsse_components(filename, monomer_name):
    energies = json.loads(open(filename, "r").read())
    energy_dfs = create_df(energies)
    energy_types_df = pd.DataFrame(
        columns=energy_dfs["Total_energy"].columns.tolist() + ["type"]
    )
    for energy_type, df in energy_dfs.items():
        df = df.copy()
        df["type"] = energy_type
        energy_types_df = pd.concat([energy_types_df, df], axis=0)

    energy_types_df.reset_index(drop=True, inplace=True)
    energy_dfs = run_all_calculations(energy_types_df)
    tmp_df = energy_types_df[
        ["type", "BSSE Monomer Sum", "BSSE Dimer Sum", "BSSE Trimer Sum"]
    ]
    tmp_df.plot(x="type", kind="barh")
    plt.title(f"BSSE components for Tetramer ({monomer_name})")
    plt.tight_layout()


def plot_interaction_v_bsse(filename, monomer_name):
    energies = json.loads(open(filename, "r").read())
    energy_dfs = create_df(energies)
    energy_types_df = pd.DataFrame(
        columns=energy_dfs["Total_energy"].columns.tolist() + ["type"]
    )

    for energy_type, df in energy_dfs.items():
        df = df.copy()
        df["type"] = energy_type
        energy_types_df = pd.concat([energy_types_df, df], axis=0)
    energy_types_df.reset_index(drop=True, inplace=True)
    energy_dfs = run_all_calculations(energy_types_df)
    energy_dfs = energy_dfs.set_index("type")
    energy_dfs = energy_dfs.loc[
        [
            "Total_energy",
            "SCF_energy",
            "MP2_correlation_energy",
            "CCSD_correlation_energy",
            "(T)_correlation_energy",
        ]
    ]
    lines, labels = plot_dual_axis(
        energy_dfs["Total Interaction Energy (BSSE free)"],
        energy_dfs["Total BSSE"],
        energy_dfs.index,
        title=f"Interaction Energy vs BSSE for ({monomer_name})",
        x_label="Side Length (Å)",
        y1_label="Interaction Energy (kcal/mol)",
        y2_label="BSSE (kcal/mol)",
        label1="Interaction Energy",
        label2="BSSE",
    )
    plt.tight_layout()
    plt.legend(lines, labels, loc="upper left")
    return energy_dfs


def plot_bsse_composition_trimer(filename, monomer_name):
    energies = json.loads(open(filename, "r").read())
    new_energies = defaultdict(dict)
    labels = {
        "BSSE monomers (one body approximation)": r"$\Sigma_{i} \theta(ABC)_{i}$",
        "BSSE monomers (two body approximation)": r"$\Sigma_{i}(\theta(ABC)_{i} - \theta(ij)_{i})$",
        "BSSE Dimers (two body approximation)": r"$\Sigma_{i<j} \theta(ABC)_{ij}$",
    }
    for side_length, calculations in energies.items():
        for calculation_type, energies_types in calculations.items():
            match = re.search(r"\((.*?)\)", calculation_type)
            if match:
                result = match.group(1)  # 'afd'
            else:
                raise Exception
            monomers = result.replace("_", "")
            basis = (
                calculation_type.replace(f"E({result})", "")
                .replace("_", "")
                .replace("E(", "")
                .replace(")", "")
            )
            new_energies[side_length][f"E({basis})_{monomers}"] = energies_types
    energy_dfs = create_df(new_energies)
    energy_dfs = {
        k: create_BSSEs_and_interaction_energies(df) for k, df in energy_dfs.items()
    }
    energy_dfs = {
        k: TrimerBsseCalculator.calculate_all_bsse(df) for k, df in energy_dfs.items()
    }
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    for ax, (k, df) in zip(axes, energy_dfs.items()):
        formated_title = k.replace("_", " ")
        tmp_df = df[["side_length"] + list(labels.keys())]
        tmp_df.plot(x="side_length", ax=ax, kind="bar")
        ax.set_title(f"BSSE for {formated_title} ({monomer_name})")
        ax.legend(
            labels=[labels[col] for col in tmp_df.columns[1:]],
        )
    plt.suptitle(f"BSSE components for Trimer ({monomer_name})", fontsize=20)
    plt.tight_layout()
    fig.delaxes(axes[-1])


def plot_bsse_composition_tetramer(filename, monomer_name):
    energy_dfs = plot_interaction_v_bsse(filename, monomer_name)
    energy_dfs = TetramerBSSECalculator.calculate_all_bsse(energy_dfs)
    monomer_name = "HF"
    labels = {
        "BSSE monomers (one body approximation)": r"$\Sigma_{i} \theta_{i}(ABCD)$",
        "BSSE monomers (two body approximation)": r"$\Sigma_{i}(\theta_{i}(ABCD) - \theta_{i}(ij))$",
        "BSSE monomers (three body approximation)": r"$\Sigma_{i}(\theta_{i}(ABCD) - \theta_{i}(ijk))$",
        "BSSE Dimers (two body approximation)": r"$\Sigma_{i<j} \theta_{ij}(ABCD)$",
        "BSSE Dimers (three body approximation)": r"$\Sigma_{i<j}(\theta_{ij}(ABCD) - \theta_{ij}(ijk))$",
        "BSSE Trimers (three body approximation)": r"$\Sigma_{i<j<k} \theta_{ijk}(ABCD)$",
    }
    tmp_df = energy_dfs[list(labels.keys())]
    tmp_df.index = tmp_df.index.map(
        lambda x: x.replace("_", " ").replace("energy", "").replace("correlation", "")
    )

    ax = tmp_df.plot(kind="barh", figsize=(7, 5))
    # Title and legend
    plt.title(f"BSSE components for Tetramer ({monomer_name})")
    ax.legend(
        labels=[labels[col] for col in tmp_df.columns],
    )
    ax.set_ylabel("Energy Component (kcal/mol)")
    plt.tight_layout()
    plt.show()


class TrimerBsseCalculator:
    @staticmethod
    def switch_convention(df):
        df_copy = df.copy()
        for column in df.columns:
            match = re.search(r"\((.*?)\)", column)
            if match:
                result = match.group(1)  # 'afd'
            else:
                continue
            monomers = result.replace("_", "")
            basis = (
                column.replace(f"E({result})", "")
                .replace("_", "")
                .replace("E(", "")
                .replace(")", "")
            )
            df_copy[f"E({basis})_{monomers}"] = df[column]
        return df_copy

    @staticmethod
    def calculate_monomer_bsse_one_approx(df):
        total_dimer_bsse = pd.Series(0, index=df.index)
        for monomer in ["A", "B", "C"]:
            total_dimer_bsse += df[f"E({monomer})_{monomer}"] - df[f"E(ABC)_{monomer}"]
        df["BSSE monomers (one body approximation)"] = total_dimer_bsse
        return df

    @staticmethod
    def calculate_monomer_bsse_two_approx(df):
        total_dimer_bsse = pd.Series(0, index=df.index)
        for monomer in ["A", "B", "C"]:
            two_body_approx = pd.Series(0, index=df.index)
            for dimer in combinations(["A", "B", "C"], 2):
                if monomer in dimer:
                    two_body_approx += df[f"E({dimer[0]}{dimer[1]})_{monomer}"] / 2
            total_dimer_bsse += two_body_approx - df[f"E(ABC)_{monomer}"]
        df["BSSE monomers (two body approximation)"] = total_dimer_bsse
        return df

    @staticmethod
    def calculate_dimer_two_approx(df):
        total_dimer_bsse = pd.Series(0, index=df.index)
        for dimer in combinations(["A", "B", "C"], 2):
            total_dimer_bsse += (
                df[f"Delta_E({dimer[0]}{dimer[1]})_{dimer[0]}{dimer[1]}"]
                - df[f"Delta_E(ABC)_{dimer[0]}{dimer[1]}"]
            )
        df["BSSE Dimers (two body approximation)"] = total_dimer_bsse
        return df

    @staticmethod
    def calculate_all_bsse(df, switch_convention=False):
        if switch_convention:
            df = TrimerBsseCalculator.switch_convention(df)
        df = TrimerBsseCalculator.calculate_monomer_bsse_one_approx(df)
        df = TrimerBsseCalculator.calculate_monomer_bsse_two_approx(df)
        df = TrimerBsseCalculator.calculate_dimer_two_approx(df)
        return df


class TetramerBSSECalculator:
    @staticmethod
    def calculate_monomer_bsse_one_approx(df):
        bsse = pd.Series(0, index=df.index)
        for monomer in ["A", "B", "C", "D"]:
            bsse += df[f"E({monomer})_{monomer}"] - df[f"E(ABCD)_{monomer}"]
        df["BSSE monomers (one body approximation)"] = bsse
        return df

    @staticmethod
    def calculate_monomer_bsse_two_approx(df):
        bsse = pd.Series(0, index=df.index)
        for monomer in ["A", "B", "C", "D"]:
            two_body_approx = pd.Series(0, index=df.index)
            for dimer in combinations(["A", "B", "C", "D"], 2):
                if monomer in dimer:
                    two_body_approx += df[f"E({dimer[0]}{dimer[1]})_{monomer}"] / 3
            bsse += two_body_approx - df[f"E(ABCD)_{monomer}"]
        df["BSSE monomers (two body approximation)"] = bsse
        return df

    @staticmethod
    def calculate_monomer_bsse_three_approx(df):
        bsse = pd.Series(0, index=df.index)
        for monomer in ["A", "B", "C", "D"]:
            three_body_approx = pd.Series(0, index=df.index)
            for trimer in combinations(["A", "B", "C", "D"], 3):
                if monomer in trimer:
                    three_body_approx += (
                        df[f"E({trimer[0]}{trimer[1]}{trimer[2]})_{monomer}"] / 3
                    )
            bsse += three_body_approx - df[f"E(ABCD)_{monomer}"]
        df["BSSE monomers (three body approximation)"] = bsse
        return df

    @staticmethod
    def calculate_dimer_two_approx(df):
        bsse = pd.Series(0, index=df.index)
        for dimer in combinations(["A", "B", "C", "D"], 2):
            bsse += (
                df[f"Delta_E({dimer[0]}{dimer[1]})_{dimer[0]}{dimer[1]}"]
                - df[f"Delta_E(ABCD)_{dimer[0]}{dimer[1]}"]
            )
        df["BSSE Dimers (two body approximation)"] = bsse
        return df

    @staticmethod
    def calculate_dimer_three_approx(df):
        bsse = pd.Series(0, index=df.index)
        for dimer in combinations(["A", "B", "C", "D"], 2):
            three_body_approx = pd.Series(0, index=df.index)
            for trimer in combinations(["A", "B", "C", "D"], 3):
                if set(dimer).issubset(set(trimer)):
                    three_body_approx += (
                        df[
                            f"Delta_E({trimer[0]}{trimer[1]}{trimer[2]})_{dimer[0]}{dimer[1]}"
                        ]
                        / 2
                    )
            bsse += three_body_approx - df[f"Delta_E(ABCD)_{dimer[0]}{dimer[1]}"]
        df["BSSE Dimers (three body approximation)"] = bsse
        return df

    @staticmethod
    def calculate_trimer_three_approx(df):
        bsse = pd.Series(0, index=df.index)
        for trimer in combinations(["A", "B", "C", "D"], 3):
            bsse += (
                df[
                    f"Delta_E({trimer[0]}{trimer[1]}{trimer[2]})_{trimer[0]}{trimer[1]}{trimer[2]}"
                ]
                - df[f"Delta_E(ABCD)_{trimer[0]}{trimer[1]}{trimer[2]}"]
            )
        df["BSSE Trimers (three body approximation)"] = bsse
        return df

    @staticmethod
    def calculate_all_bsse(df, switch_convention=False):
        if switch_convention:
            df = TetramerBSSECalculator.switch_convention(df)
        df = TetramerBSSECalculator.calculate_monomer_bsse_one_approx(df)
        df = TetramerBSSECalculator.calculate_monomer_bsse_two_approx(df)
        df = TetramerBSSECalculator.calculate_monomer_bsse_three_approx(df)
        df = TetramerBSSECalculator.calculate_dimer_two_approx(df)
        df = TetramerBSSECalculator.calculate_dimer_three_approx(df)
        df = TetramerBSSECalculator.calculate_trimer_three_approx(df)
        return df
