# H₂O Cluster BSSE Dataset

This directory contains quantum-chemistry calculations for **water (H₂O) clusters**
of increasing size, with a focus on analyzing **Basis Set Superposition Error (BSSE)**
as a function of intermolecular geometry. All calculations are performed at the
CCSD(T)/aug-cc-pVDZ level of theory using NWChem. All monomers are treated as
**rigid** throughout the dataset.

The dataset is organized by cluster size. Each subdirectory `n/` contains the
complete set of configurations, NWChem input/output files, parsed data frames, and
visualization utilities for that cluster. A detailed description of each cluster's
geometry generation, naming conventions, and folder layout is given in the
corresponding `n/README.md`.

Energy parsing across all clusters is handled by `../build_raw_data.py`, which reads
the NWChem output tarballs directly and extracts SCF, MP2, and CCSD(T) energies into
a nested dictionary keyed by molecule, cluster, motion type, configuration, method,
and basis set. See `../README.md` for details on the parser.

---

## Cluster 2 — H₂O Dimer (`2/`)

The dimer dataset systematically explores BSSE as a function of intermolecular
geometry by sampling two classes of relative motion between the two rigid H₂O
monomers.

### Geometry generation

All geometries are defined by the oxygen–oxygen displacement vector

$$\mathbf{r} = \mathbf{r}_2 - \mathbf{r}_1$$

where $\mathbf{r}_1$ and $\mathbf{r}_2$ are the Cartesian coordinates of the oxygen
atoms in monomer 1 and monomer 2 respectively. Coordinates are given in **Ångströms**.

**Translational configurations** — the relative orientation of the two monomers is
held fixed and the O–O displacement vector is varied in increments of **0.25 Å**.
Not all grid points were retained; configurations where the translational sampling
did not provide sufficient geometric coverage were excluded from the final dataset.

**Rotational–radial configurations** — both the relative orientation and the
intermolecular separation are varied. Euler angles are sampled randomly and the
O–O separation distance R is varied independently. Configurations with an
insufficient number of distinct R values per orientation were excluded from the
analysis.

A detailed account of which configurations were excluded and why for both datasets
is provided in `2/data_frames_plots/README.md`.

### Configuration labels

Translational configurations are labeled `n_x_y_z`, where `n` is the configuration
index and `(x, y, z)` are the Cartesian components of **r** in Ångströms.

Rotational–radial configurations are labeled `n_α_β_γ_R`, where `(α, β, γ)` are
Euler angles in degrees and `R` is the O–O separation in Ångströms.

### Folder contents

See `2/README.md` for the complete description of the internal folder structure,
NWChem file organization, parsed data frames, and visualization utilities.

---

## Cluster 3 — H₂O Trimer (`H2O_trimers/`)

*In progress.* This section will describe the geometry generation strategy,
configuration labeling convention, and folder layout for the trimer dataset once
the sampling and calculations are complete.

---

## Current Status

| Cluster | Directory      | Status      |
|---------|----------------|-------------|
| Dimer   | `2/`           | Complete    |
| Trimer  | `H2O_trimers/` | In progress |
