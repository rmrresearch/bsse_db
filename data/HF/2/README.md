# HF Dimer Dataset (`HF/2/`)

This directory contains the reorganized NWChem input/output files for
the HF (hydrogen fluoride) dimer, along with the scripts used to
perform the reorganization. The data is used to compute the
basis set superposition error (BSSE) using the Variational
Many-body Function Correction (VMFC) method.

---

## Directory Structure

```
HF/2/
├── original_data/
│   ├── HF_dimers.tar           ← translational dataset backup (aug-cc-pVDZ)
│   ├── HF_pvtz_dimers.tar      ← radial dataset backup (aug-cc-pVTZ)
│   └── HF_files.tar            ← rotational dataset backup (aug-cc-pVDZ, partial)
├── translational/              ← 268 configs, aug-cc-pVDZ
├── radial/                     ← 3 configs, aug-cc-pVTZ
├── rotational/                 ← 152 configs, aug-cc-pVDZ
├── pymol_xyz/                  ← 268 .xyz files (translational, n_x_y_z.xyz)
├── pymol_radial_xyz/           ← 3 .xyz files (radial, n_0.0_0.0_R.xyz)
├── pymol_rotational_xyz/       ← 152 .xyz files (rotational, n_alpha_beta_gamma_R.xyz)
├── load_dimers.pml             ← PyMOL loader for translational configs
├── load_radial.pml             ← PyMOL loader for radial configs
├── load_rotational.pml         ← PyMOL loader for rotational configs
├── data_frames_plots/          ← pipeline output (CSVs, plots)
├── HF_dimer_energies_euler.json ← pre-computed rotational BSSE (304 configs)
├── HF_dimer_energies_pvtz.json  ← pre-computed pvtz BSSE energies
├── README.md                   ← this file
└── scripts/
    ├── copy_rename_files_HF.py
    ├── copy_rename_files_HF_pvtz.py
    ├── copy_rename_files_HF_rot.py
    ├── make_pymol_xyz_HF.py
    ├── make_pymol_radial_xyz_HF.py
    ├── make_pymol_rotational_xyz_HF.py
    └── legacy/                 ← original config_scripts/ and misc scripts
```

---

## Monomer Geometry

All HF dimer calculations use rigid monomers. Monomer A is fixed at
the origin for all configurations:

```
H   0.000   0.000  -0.924   (Angstroms)
F   0.000   0.000   0.000
```

Monomer B has the same bond geometry as monomer A, displaced and/or
rotated in space according to the configuration parameters.

**Note on B_B = A_A:** Since both HF monomers have identical bond
geometries, no separate `B_B` calculation was performed. The BSSE
scripts use `A_A` for both monomers when computing the correction.

---

## File Naming Convention

In each innermost directory (`CCSD_T/aug-cc-pvdz/` or
`CCSD_T/aug-cc-pvtz/`) the following files exist:

| File | Description |
|------|-------------|
| `A_A.nw` / `A_A.out` | Monomer A in its own basis (no ghost functions) |
| `AB_AB.nw` / `AB_AB.out` | Full dimer AB in the dimer basis |
| `A_AB.nw` / `A_AB.out` | Monomer A in the dimer basis (ghost functions at B) |
| `B_AB.nw` / `B_AB.out` | Monomer B in the dimer basis (ghost functions at A) |

File naming follows the `<real>_<basis>.out` convention from Issue #46:
- `A_A` = monomer A (real=A, basis=A)
- `AB_AB` = dimer AB (real=AB, basis=AB)
- `A_AB` = monomer A computed in dimer basis (real=A, basis=AB)
- `B_AB` = monomer B computed in dimer basis (real=B, basis=AB)

---

## Configuration Naming Convention

Configurations follow the PI's `<n>_<desc>` convention (Issue #46):
- `n` is the 0-based configuration index
- `desc` encodes the geometric parameters

This makes the mapping self-documenting and eliminates the need for
a separate parameter lookup table (though the `parameters` dictionary
built by the parser also stores this information).

---

## Datasets

### 1. Translational (aug-cc-pVDZ, 268 configurations)

**Source:** `original_data/HF_dimers.tar`

**Description:** Monomer A is fixed at the origin. Monomer B is
translated in 3D space with the F atom positioned at coordinates
(x, y, z) relative to the F atom of monomer A.

**Configuration label:** `n_x_y_z`

**Parameters:**
- x, y, z ∈ {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75} Å
- Grid is incomplete — not all 7×7×7 = 343 combinations are present
- Configurations are in lexicographic order by (x, y, z) tuple

**Path structure:**
```
translational/n_x_y_z/CCSD_T/aug-cc-pvdz/
```

**Example:**
```
translational/0_0.25_0.25_1.5/CCSD_T/aug-cc-pvdz/A_A.out
```

**File counts:** 268 configs × 4 files × 2 extensions = 2144 files

---

### 2. Radial (aug-cc-pVTZ, 3 configurations)

**Source:** `original_data/HF_pvtz_dimers.tar`

**Description:** Monomer A is fixed at the origin. Monomer B is
displaced along the z-axis at distance R. This dataset uses the
larger aug-cc-pVTZ basis set to study BSSE as a function of
separation distance.

**Configuration label:** `n_0.0_0.0_R`

**Parameters:**
- R ∈ {2.5, 3.0, 3.5} Å along the z-axis
- Configurations are ordered by ascending R value
- **Note:** The HF radial dataset is incomplete — only 3 of the 5 R
  values present in the Ne radial dataset (1.5–3.5 Å) are available.
  R = 1.5 Å and R = 2.0 Å are missing from the original data.

**Path structure:**
```
radial/n_0.0_0.0_R/CCSD_T/aug-cc-pvtz/
```

**Example:**
```
radial/0_0.0_0.0_2.5/CCSD_T/aug-cc-pvtz/A_A.out
```

**File counts:** 3 configs × 4 files × 2 extensions = 24 files

---

### 3. Rotational (aug-cc-pVDZ, 152 configurations)

**Source:** `original_data/HF_files.tar` (partially corrupted)

**Description:** Monomer A is fixed at the origin. Monomer B is
rotated using three Euler angles (alpha, beta, gamma) and placed at
separation distance R. This dataset is used to assess the orientation
dependence of BSSE.

**Configuration label:** `n_alpha_beta_gamma_R`
- Angles are in **degrees**
- R is in Angstroms

**Parameters:**
- alpha, beta, gamma ∈ {0.0°, 90.0°, 180.0°, 270.0°}
- R ∈ {1.75, 2.25, 2.75, 3.25, 3.75} Å
  (except alpha=0.0°: R ∈ {2.75, 3.25, 3.75} Å only)
- Configurations are in lexicographic order by (alpha, beta, gamma, R)

**Path structure:**
```
rotational/n_alpha_beta_gamma_R/CCSD_T/aug-cc-pvdz/
```

**Example:**
```
rotational/0_0.0_0.0_0.0_2.75/CCSD_T/aug-cc-pvdz/A_A.out
```

**Dataset completeness note:**
The full rotational dataset contains 304 configurations. Only 152 raw
output files were recovered from `HF_files.tar` due to file corruption
(alpha = 180.0° and alpha = 270.0° are missing). The complete
pre-computed BSSE energies for all 304 configurations are available in
`HF_dimer_energies_euler.json`. Per PI guidance, 152 configurations
are sufficient to demonstrate orientation independence of BSSE.

**File counts:** 152 configs × 4 files × 2 extensions = 1216 files

---

## Pre-computed Energy Files

### `HF_dimer_energies_euler.json`
Contains pre-computed BSSE results for all **304** rotational
configurations. Structure:

```
d[alpha_rad][beta_rad][gamma_rad][R] = {
    "Total_energy_Delta_E_AB_AB": ...,
    "Total_energy_BSSE_A": ...,
    "Total_energy_BSSE_B": ...,
    "SCF_energy_Delta_E_AB_AB": ...,
    ...
}
```

Angles are stored as radian strings. This file covers the 152 configs
missing from the raw output files.

### `HF_dimer_energies_pvtz.json`
Contains pre-computed BSSE results for the pvtz dataset.

---

## Scripts

### `scripts/copy_rename_files_HF.py`
Reorganizes the translational dataset from the original `HF_dimers/`
flat folder structure into the standardized `translational/n_x_y_z/`
layout.

**Run from:** `HF/2/`
```bash
python3 scripts/copy_rename_files_HF.py
```
**Result:** 268 configurations processed ✅

---

### `scripts/copy_rename_files_HF_pvtz.py`
Extracts the radial dataset directly from `original_data/HF_pvtz_dimers.tar`
into the standardized `radial/n_0.0_0.0_R/` layout.

**Run from:** `HF/2/`
```bash
python3 scripts/copy_rename_files_HF_pvtz.py
```
**Result:** 3 configurations processed ✅

---

### `scripts/copy_rename_files_HF_rot.py`
Reorganizes the rotational dataset from the `HF_files_extracted/` folder
into the standardized `rotational/n_alpha_beta_gamma_R/` layout.
Converts Euler angles from radians to degrees.

**Run from:** `HF/2/`
```bash
python3 scripts/copy_rename_files_HF_rot.py
```
**Result:** 152 configurations processed ✅

---

### `scripts/make_pymol_xyz_HF.py`
Generates `pymol_xyz/n_x_y_z.xyz` files and `load_dimers.pml` from
the translational dataset.

**Run from:** `HF/2/`
```bash
python3 scripts/make_pymol_xyz_HF.py
```
**Result:** 268 XYZ files written ✅

---

### `scripts/make_pymol_radial_xyz_HF.py`
Generates `pymol_radial_xyz/n_0.0_0.0_R.xyz` files and `load_radial.pml`
from the radial dataset.

**Run from:** `HF/2/`
```bash
python3 scripts/make_pymol_radial_xyz_HF.py
```
**Result:** 3 XYZ files written ✅

---

### `scripts/make_pymol_rotational_xyz_HF.py`
Generates `pymol_rotational_xyz/n_alpha_beta_gamma_R.xyz` files and
`load_rotational.pml` from the rotational dataset.

**Run from:** `HF/2/`
```bash
python3 scripts/make_pymol_rotational_xyz_HF.py
```
**Result:** 152 XYZ files written ✅

---

### `scripts/legacy/`
Contains the original configuration setup scripts used to generate the
HF dimer calculations:
- `batch_euler_angles_config_folder.py`
- `batch_setup_config_folder.py`
- `rename_cp_rotational_radial.py`
- `rename_cp_translational.py`
- `setup_config_folder.py`
- `setup_euler_config_folder.py`
- `load_hf_dimers.py`
- `make_contour_plots.ipynb`

These scripts are preserved for reference only and are superseded by
the reorganization scripts above.

---

## Level of Theory

All calculations use NWChem 7.2.3 at the CCSD(T) level of theory.

| Dataset | Basis Set |
|---------|-----------|
| Translational | aug-cc-pVDZ |
| Radial | aug-cc-pVTZ |
| Rotational | aug-cc-pVDZ |

---

## PyMOL Visualization

XYZ files are provided for all three datasets. All `.pml` loader scripts
must be run from `HF/2/`.

### Translational (268 configurations)
```bash
pymol load_dimers.pml
```
XYZ files in `pymol_xyz/` named `n_x_y_z.xyz`. Generated by
`scripts/make_pymol_xyz_HF.py`.

### Radial (3 configurations)
```bash
pymol load_radial.pml
```
XYZ files in `pymol_radial_xyz/` named `n_0.0_0.0_R.xyz`. Generated by
`scripts/make_pymol_radial_xyz_HF.py`.

### Rotational (152 configurations)
```bash
pymol load_rotational.pml
```
XYZ files in `pymol_rotational_xyz/` named `n_alpha_beta_gamma_R.xyz`.
Generated by `scripts/make_pymol_rotational_xyz_HF.py`.

All loaders use `util.cbaw` color scheme (element-based, visible on black
background) with ball-and-stick representation.