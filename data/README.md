## Data Directory Overview

This `data/` directory contains the raw NWChem quantum-chemistry results, parsing utilities, and serialized datasets used for **Basis Set Superposition Error (BSSE)** analysis.

The directory is organized so that raw calculations are preserved, while parsed and validated data can be reused efficiently in later analysis steps.

---

### Molecule Folders

Each molecule (or atom) has its own directory inside `data/`, for example:

* `H2O/`
* `HF/`
* `Ne/`

These folders contain all raw NWChem calculations associated with that system. Typically, each molecule directory includes one or more compressed archives (`*.tar.gz`), where the leading number indicates the cluster size:

* `2.tar.gz` → dimers
* `3.tar.gz` → trimers (if present)

#### Dimer Archives (`2.tar.gz`)

Dimer archives are organized by **motion type** and **configuration**:

* `translational/`
  Cartesian translations between monomers.

* `rotational_radial/`
  Rotational scans defined by Euler angles and an intermolecular separation.

Configuration folders encode the geometric parameters directly in their names:

* Translational example: `0_0.25_0.25_1.5` representing `(config_id, x, y, z)`.
* Rotational–radial example: `3_0.0_0.0_90.0_2.76` representing `(config_id, α, β, γ, R)`.

#### Trimer Archives (`3.tar.gz`)

Trimer archives have no motion type level. Configurations are organized directly by O-O distance:

```
3/<config_id>_<distance>/<method>/<basis>/
```

For example:
```
3/0_2.75/CCSD_T/aug-cc-pvdz/
```

Configuration folders encode the O-O distance directly in their names, e.g. `0_2.75` representing `(config_id, O-O distance in Angstroms)`.

Within each configuration directory, calculations are organized as:

```
<configuration>/<method>/<basis>/
```

containing the corresponding NWChem output (`.out`) files for the BSSE components (e.g., `A_A`, `A_AB`, `ABC_ABC`).

Only configurations for which **all required methods converge** are retained for BSSE analysis.

---

### Per-Molecule Documentation

Each molecule directory contains one or more README files that provide **detailed scientific context** for the calculations performed.

For example, in the `H2O/` directory:

* `README_H2O.md` describes the configuration space, geometry definitions, and how translational and rotational datasets were generated for the dimers.
* `3/README.md` summarizes the trimer directory layout, configuration mapping, and provides utilities or guidance for visualization.

These files are the **authoritative reference** for interpreting configuration labels and geometric parameters used by the parsers and stored in the pickled data.

---

### Python Scripts

The `data/` directory contains Python scripts used to parse, validate, and store the results of the quantum-chemistry calculations. Scripts are named with a trailing `_<cluster_size>` suffix for cluster-specific parsers, making it clear which cluster size they target.

* **`tar_output_discovery.py`**
  Dynamically locates and streams `.out` files inside compressed `.tar.gz` archives for dimers (cluster size 2) without extracting them to disk.

* **`tar_output_discovery_3.py`**
  Same as above but for trimers (cluster size 3). Uses non-recursive `glob` instead of `rglob` to avoid picking up nested backup archives inside the molecule folders.

* **`build_raw_data.py`**
  Parses total SCF, MP2, and CCSD(T) energies from dimer output files, enforces convergence requirements, and constructs the nested `raw_data` dictionary used for BSSE analysis. Non-converged configurations are tracked separately for diagnostic purposes.

* **`build_raw_data_3.py`**
  Same as above but for trimers (cluster size 3). Adapted to handle the trimer directory structure (no motion type level) and extracts O-O distance as the configuration parameter.

* **`pickle_data.py`**
  Serializes the converged `raw_data` and `parameters` dictionaries into molecule-specific pickle files for fast reuse in later analysis steps. Supports an optional `cluster` argument to distinguish pickle files by cluster size.

These scripts are designed to be reusable across different molecules. Cluster-specific scripts (e.g., `build_raw_data_3.py`) can be adapted for HF and Ne trimers by changing the molecule argument.

---

### Pickle Files

After successful parsing, the converged data are saved as pickle files in the `data/` directory:

* **`raw_data_<MOLECULE>.pickle`**
  Example: `raw_data_H2O.pickle` — dimer data for H2O.

* **`raw_data_<MOLECULE>_<CLUSTER>.pickle`**
  Example: `raw_data_H2O_3.pickle` — trimer data for H2O.

Each pickle file contains a bundle with two keys:
* `'raw_data'` — the nested energy dictionary
* `'parameters'` — configuration metadata (distances, angles, etc.)

The `raw_data` dictionary is keyed by:
```
raw_data[molecule][cluster][config_key][level_of_theory][basis][real_basis]
```

Stored energy information includes:

* SCF total energies
* MP2 correlation energies (total MP2 minus total SCF)
* CCSD(T) correlation components (total CCSD(T) minus total MP2)

The `parameters` dictionary is keyed by:
```
parameters[molecule][cluster][config_key] = {
    'config_label': <full folder name>,
    'params': <dict of geometric parameters>
}
```

Only **fully converged configurations** are included. Any configuration with missing or non-converged calculations is excluded to ensure that all BSSE computations are consistent and physically meaningful.

---

This organization allows the parsing step to be performed once, after which all BSSE evaluation, plotting, and comparative analysis can be carried out efficiently using the pickled datasets.
