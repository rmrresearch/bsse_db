# data_frames_plots

Scripts and data for translational and rotational BSSE scans of the H₂O dimer
at the SCF, MP2, and CCSD(T) levels with the aug-cc-pVDZ basis set.

---

## Directory structure

```
data_frames_plots/
│
├── build_rot_dataframe.py          # data generation
├── build_trnl_dataframe.py
│
├── analyze_rot_coverage.py         # analysis
├── analyze_trnl_grid_coverage.py
│
├── plot_rot_bsse_vs_R.py           # plotting
├── plot_rot_contours.py
├── plot_trnl_contours.py
│
├── compare_methods_bsse_vs_R.py    # comparison
│
├── rot_SCF_aug-cc-pvdz.csv         # data
├── rot_MP2_aug-cc-pvdz.csv
├── rot_CCSD_T_aug-cc-pvdz.csv
├── trnl_SCF_aug-cc-pvdz.csv
├── trnl_MP2_aug-cc-pvdz.csv
└── trnl_CCSD_T_aug-cc-pvdz.csv
```

---

## CSV file formats

### Rotational — `rot_{METHOD}_{BASIS}.csv`

| Column    | Type  | Description                                     |
|-----------|-------|-------------------------------------------------|
| `R`       | float | Intermolecular distance (Ångströms)              |
| `alpha`   | float | Euler angle α of monomer 2 (degrees)            |
| `beta`    | float | Euler angle β of monomer 2 (degrees)            |
| `gamma`   | float | Euler angle γ of monomer 2 (degrees)            |
| `bsse`    | float | BSSE correction energy (Hartree)                |
| `no_vmfc` | float | Interaction energy without VMFC (Hartree)       |
| `vmfc`    | float | VMFC-corrected interaction energy (Hartree)     |

Monomer 1 is fixed; (α, β, γ) define the orientation of monomer 2.

### Translational — `trnl_{METHOD}_{BASIS}.csv`

| Column    | Type  | Description                                     |
|-----------|-------|-------------------------------------------------|
| `x`       | float | x-displacement of monomer 2 (Ångströms)          |
| `y`       | float | y-displacement of monomer 2 (Ångströms)          |
| `z`       | float | z-displacement of monomer 2 (Ångströms)          |
| `bsse`    | float | BSSE correction energy (Hartree)                |
| `no_vmfc` | float | Interaction energy without VMFC (Hartree)       |
| `vmfc`    | float | VMFC-corrected interaction energy (Hartree)     |

---

## Scripts

### Data generation

#### `build_rot_dataframe.py`
Reads raw rotational calculation output and writes `rot_{METHOD}_{BASIS}.csv`.

#### `build_trnl_dataframe.py`
Reads raw translational calculation output and writes `trnl_{METHOD}_{BASIS}.csv`.

---

### Analysis

#### `analyze_rot_coverage.py`
Reports how many distinct R values each (α, β, γ) orientation has across
all methods, identifies which orientations are eligible for plotting, and
optionally prints BSSE(R) value tables.

```bash
# Summary table for all methods
python3 analyze_rot_coverage.py

# Top-10 orientations per method
python3 analyze_rot_coverage.py --detailed

# BSSE(R) tables for top 5 orientations, SCF only
python3 analyze_rot_coverage.py --method SCF --print-tables

# Tables for a different value column, top 7, looser R threshold
python3 analyze_rot_coverage.py --method MP2 --print-tables --top 7 --minR 3 --value vmfc
```

| Flag | Default | Description |
|------|---------|-------------|
| `--method` | `all` | `SCF`, `MP2`, `CCSD_T`, or `all` |
| `--basis` | `aug-cc-pvdz` | Basis set suffix in CSV filename |
| `--minR` | `5` | Minimum distinct R values to be considered eligible |
| `--detailed` | off | Show top-10 orientations per method |
| `--print-tables` | off | Print BSSE(R) value tables for top orientations |
| `--top` | `5` | Number of orientations shown in `--print-tables` |
| `--value` | `bsse` | Column shown in tables: `bsse`, `no_vmfc`, `vmfc` |

---

#### `analyze_trnl_grid_coverage.py`
Reports completeness of the (x, y, z) translational grid for each fixed-axis
family (fixed x, fixed y, fixed z). Identifies which slices have complete grids
and are therefore ready to plot. Grid analysis runs once since all methods share
the same geometry.

```bash
# Coverage for all methods
python3 analyze_trnl_grid_coverage.py

# Single method
python3 analyze_trnl_grid_coverage.py --method MP2
```

| Flag | Default | Description |
|------|---------|-------------|
| `--method` | `all` | `SCF`, `MP2`, `CCSD_T`, or `all` |
| `--basis` | `aug-cc-pvdz` | Basis set suffix in CSV filename |

**Known status:** only fixed-y slices (y = 2.0, 2.25, 2.5, 2.75 Å) have
complete grids. Fixed-x and fixed-z slices have missing points and cannot be
contour-plotted yet.

---

### Plotting

#### `plot_rot_bsse_vs_R.py`
Plots the chosen quantity vs R for the top N orientations with the most R
coverage. One figure per method. The R axis is labeled in **Ångströms** and
the y axis in **Hartree**.

```bash
# Top 5 orientations, SCF
python3 plot_rot_bsse_vs_R.py --method SCF --top 5

# All methods in one call, save figures
python3 plot_rot_bsse_vs_R.py --method all --top 5 --save

# Absolute values, MP2
python3 plot_rot_bsse_vs_R.py --method MP2 --abs --save

# VMFC-corrected values, CCSD(T)
python3 plot_rot_bsse_vs_R.py --method CCSD_T --value vmfc --save
```

| Flag | Default | Description |
|------|---------|-------------|
| `--method` | `SCF` | `SCF`, `MP2`, `CCSD_T`, or `all` |
| `--basis` | `aug-cc-pvdz` | Basis set suffix in CSV filename |
| `--top` | `5` | Number of orientations to plot |
| `--minR` | `5` | Minimum R values required for an orientation |
| `--value` | `bsse` | Column to plot: `bsse`, `no_vmfc`, `vmfc` |
| `--abs` | off | Plot absolute values |
| `--save` | off | Save figure as PNG |
| `--ncols` | `3` | Number of subplot columns per row |
| `--select` | all | Orientations to plot: range `10:20` or list `1 2 30 32` (1-based) |
| `--top` | `32` | Maximum pool of ranked orientations to draw from |

---

#### `plot_rot_contours.py`
Rotational contour plots.

```bash
python3 plot_rot_contours.py
```

---

#### `plot_trnl_contours.py`
x–z contour plots at fixed y. Only runs for y slices with complete grids
(y = 2.0, 2.25, 2.5, 2.75 Å). Run `analyze_trnl_grid_coverage.py` first
to confirm which slices are available. Both spatial axes (x, z) are labeled
in **Ångströms** and the colorbar in **Hartree**.

```bash
# SCF, BSSE
python3 plot_trnl_contours.py --method SCF --value bsse

# CCSD(T), VMFC-corrected, save figures
python3 plot_trnl_contours.py --method CCSD_T --value vmfc --save
```

| Flag | Default | Description |
|------|---------|-------------|
| `--method` | `SCF` | `SCF`, `MP2`, or `CCSD_T` |
| `--basis` | `aug-cc-pvdz` | Basis set suffix in CSV filename |
| `--value` | `bsse` | Column to plot: `bsse`, `no_vmfc`, `vmfc` |
| `--save` | off | Save figures as PNG |

---

### Comparison

#### `compare_methods_bsse_vs_R.py`
Plots SCF, MP2, and CCSD(T) on the same axes for the same orientations,
one subplot per orientation. Useful for seeing the effect of electron
correlation on BSSE directly. The R axis is labeled in **Ångströms** and
the y axis in **Hartree**.

```bash
# All eligible orientations, default 3-column grid
python3 compare_methods_bsse_vs_R.py

# Plot orientations 1 through 9 in a 3x3 grid
python3 compare_methods_bsse_vs_R.py --select 1:9 --ncols 3

# Plot specific orientations by index
python3 compare_methods_bsse_vs_R.py --select 1 2 30 32

# Orientations 10 to 20, absolute values, save
python3 compare_methods_bsse_vs_R.py --select 10:20 --abs --save

# All 32 orientations in a 6-column grid, save
python3 compare_methods_bsse_vs_R.py --select 1:32 --ncols 6 --save

# SCF and MP2 only, 2 columns per row
python3 compare_methods_bsse_vs_R.py --methods SCF MP2 --select 1:6 --ncols 2 --abs --save
```

| Flag | Default | Description |
|------|---------|-------------|
| `--methods` | `SCF MP2 CCSD_T` | Methods to overlay |
| `--basis` | `aug-cc-pvdz` | Basis set suffix in CSV filename |
| `--top` | `3` | Number of orientations (subplots) |
| `--minR` | `5` | Minimum R values required for an orientation |
| `--abs` | off | Plot absolute values |
| `--save` | off | Save figure as PNG |
| `--ncols` | `3` | Number of subplot columns per row |
| `--select` | all | Orientations to plot: range `10:20` or list `1 2 30 32` (1-based) |
| `--top` | `32` | Maximum pool of ranked orientations to draw from |

---

## Typical workflow

```bash
# 1. Check rotational coverage
python3 analyze_rot_coverage.py --detailed

# 3. Check translational grid completeness
python3 analyze_trnl_grid_coverage.py

# 4. Plot rotational BSSE vs R per method
python3 plot_rot_bsse_vs_R.py --method SCF --top 5 --save
python3 plot_rot_bsse_vs_R.py --method MP2 --top 5 --save
python3 plot_rot_bsse_vs_R.py --method CCSD_T --top 5 --save

# 5. Compare methods directly
python3 compare_methods_bsse_vs_R.py --select 1:9 --ncols 3 --save
python3 compare_methods_bsse_vs_R.py --select 1:32 --ncols 6 --save

# 6. Translational contours
python3 plot_trnl_contours.py --method SCF --save
python3 plot_trnl_contours.py --method MP2 --save
python3 plot_trnl_contours.py --method CCSD_T --save
```
