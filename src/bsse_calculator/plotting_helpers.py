import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy.optimize import curve_fit


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
