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


def make_bsse_contour_plots(bsse_data_nwchem, monomer_name):
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
        data = bsse_data_nwchem[bsse_data_nwchem["y"] == fixed_y_value][
            ["x", "z", "BSSE_A"]
        ]
        xs = data["x"].values
        zs = data["z"].values
        bsse = data["BSSE_A"].values

        # Grid and interpolation
        xi = np.linspace(min(xs), max(xs), 50)
        zi = np.linspace(min(zs), max(zs), 50)
        Xi, Zi = np.meshgrid(xi, zi)
        bssei = griddata((xs, zs), bsse, (Xi, Zi), method="cubic")

        # Filled contour plot with fixed vmin/vmax
        contour = ax.contourf(
            Xi, Zi, bssei, levels=8, cmap="viridis", vmin=vmin, vmax=vmax
        )
        ax.set_title(f"Contour Plot for y = {fixed_y_value}")
        ax.set_xlabel(f"X Coordinate of Second {monomer_name} (Å)")
        ax.set_ylabel(f"Z Coordinate of Second {monomer_name} (Å)")

    # Create colorbar using the normalizer and colormap directly
    cbar = fig.colorbar(
        cm.ScalarMappable(norm=normalizer, cmap=cmap), ax=axes.ravel().tolist()
    )
    cbar.set_label("BSSE (kcal/mol)", labelpad=20)

    plt.show()


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


def predicted_v_actual_plot(bsse_data_nwchem, fit_gauss):
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
    plt.title("Predicted vs Actual BSSE Values")
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
                row[calculation] = (
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
        dimer = f"{monomers[0]}_{monomers[1]}"
        total_energies[dimer_basis_name] = (
            total_energies[f"E({dimer})_{dimer}"]
            - total_energies[f"E({monomers[0]})_{dimer}"]
            - total_energies[f"E({monomers[1]})_{dimer}"]
        )
        # BSSE of adding third basis
        full_basis_name = f"Delta_E({monomers[0]}{monomers[1]})_ABC"
        total_energies[full_basis_name] = (
            total_energies[f"E({dimer})_A_B_C"]
            - total_energies[f"E({monomers[0]})_A_B_C"]
            - total_energies[f"E({monomers[1]})_A_B_C"]
        )
        total_energies[f"BSSE_({dimer})"] = (
            total_energies[dimer_basis_name] - total_energies[full_basis_name]
        )
    for monomer in ["A", "B", "C"]:
        column_name = f"BSSE_({monomer})_ABC"
        total_energies[column_name] = (
            total_energies[f"E({monomer})_{monomer}"]
            - total_energies[f"E({monomer})_A_B_C"]
        )
    working_interaction_energy = total_energies["E(A_B_C)_A_B_C"].copy()
    for monomer in ["A", "B", "C"]:
        working_interaction_energy -= total_energies[f"E({monomer})_A_B_C"]
    for monomers in combinations(["A", "B", "C"], 2):
        working_interaction_energy += (
            -total_energies[f"Delta_E({monomers[0]}{monomers[1]})_ABC"]
            + total_energies[
                f"Delta_E({monomers[0]}{monomers[1]})_{monomers[0]}{monomers[1]}"
            ]
        )
    total_energies["Total Interaction Energy (BSSE free)"] = working_interaction_energy
    working_total_interaction_energy_wo_bsse_correction = total_energies[
        "E(A_B_C)_A_B_C"
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

    y_delta_e = np.array(total_energies["Delta_E(A_B_C)_A_B_C"])[sorted_indices]
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


def plot_dual_axis(
    series1: pd.Series,
    series2: pd.Series,
    x_series: pd.Series,
    *,
    ax1: Optional[plt.Axes] = None,
    title: str = "Dual Axis Plot",
    x_label: str = "X-axis",
    y1_label: str = "Left Y-axis",
    y2_label: str = "Right Y-axis",
    label1: str = "Series 1",
    label2: str = "Series 2",
):
    fig, ax1 = plt.subplots(figsize=(10, 6)) if ax1 is None else (None, ax1)

    ax1.plot(
        x_series.values, series1.values, color="tab:blue", label=label1, marker="o"
    )
    ax1.set_ylabel(y1_label)
    ax1.tick_params(axis="y")

    ax2 = ax1.twinx()
    ax2.plot(x_series.values, series2.values, color="tab:red", label=label2, marker="o")
    ax2.set_ylabel(y2_label)
    ax2.tick_params(axis="y")

    ax1.set_xlabel(x_label)
    ax1.set_title(title)

    # Optional: combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    return lines1 + lines2, labels1 + labels2
    # ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")
