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

Each archive contains calculations organized by **motion type** and **configuration**.

#### Motion Types

* `translational/`
  Cartesian translations between monomers.

* `rotational_radial/`
  Rotational scans defined by Euler angles and an intermolecular separation.

#### Configuration Directories

Configuration folders encode the geometric parameters directly in their names:

* Translational example:

  ```
  0_0.25_0.25_1.5
  ```

  representing `(config_id, x, y, z)`.

* Rotational–radial example:

  ```
  3_0.0_0.0_90.0_2.76
  ```

  representing `(config_id, α, β, γ, R)`.

Within each configuration directory, calculations are further organized as:

```
<configuration>/<method>/<basis>/
```

containing the corresponding NWChem input (`.nw`) and output (`.out`) files for the BSSE components (e.g., `A_A`, `A_AB`, `AB_AB`).

Only configurations for which **all required methods converge** are retained for BSSE analysis.

---

### Per-Molecule Documentation

Each molecule directory contains one or more README files (`README_n.md`) that provide **detailed scientific context** for the calculations performed.

For example, in the `H2O/` directory:

* `README_1.md` describes the configuration space, geometry definitions, and how translational and rotational datasets were generated.
* `README_2.md` summarizes the directory layout and provides utilities or guidance for visualization.

These files are the **authoritative reference** for interpreting configuration labels and geometric parameters used by the parsers and stored in the pickled data.

---

### Python Scripts

The `data/` directory also contains Python scripts used to parse, validate, and store the results of the quantum-chemistry calculations:

* **`tar_output_discovery.py`**
  Dynamically locates and streams `.out` files inside compressed `.tar.gz` archives without extracting them to disk.

* **`build_raw_data.py`**
  Parses total SCF, MP2, and CCSD(T) energies from output files, enforces convergence requirements, and constructs the nested `raw_data` dictionary used for BSSE analysis. Non-converged configurations are tracked separately for diagnostic purposes.

* **`pickle_data.py`**
  Serializes the converged `raw_data` dictionary into molecule-specific pickle files for fast reuse in later analysis steps.

These scripts are designed to be reusable across different molecules and cluster sizes.

---

### Pickle Files

After successful parsing, the converged data are saved as pickle files in the `data/` directory:

* **`raw_data_<MOLECULE>.pickle`**
  Example: `raw_data_H2O.pickle`

These pickle files contain **validated energy and structural information** required to compute **Basis Set Superposition Error (BSSE)** for specific molecular configurations (both translational and rotational).

Stored information includes:

* SCF total energies
* MP2 correlation energies
* CCSD(T) correlation components (relative to MP2)
* Configuration metadata (translation vectors or Euler angles and intermolecular distance)

Only **fully converged configurations** are included. Any configuration with missing or non-converged calculations is excluded to ensure that all BSSE computations are consistent and physically meaningful.

---

This organization allows the parsing step to be performed once, after which all BSSE evaluation, plotting, and comparative analysis can be carried out efficiently using the pickled datasets.
