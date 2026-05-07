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

---

## File Naming Convention

All NWChem files follow the `<real>_<basis>` convention defined in GitHub Issue #46:

- `<real>` — the monomer(s) with real atoms
- `<basis>` — the monomer(s) contributing basis functions
- `.nw` — NWChem input file
- `.out` — NWChem output file

**Example:** monomer A computed in the full trimer basis → `A_ABC.nw` / `A_ABC.out`

---

## Energy Parsing Convention

NWChem energies are stored as additive components across all cluster sizes:

| Key | Quantity stored |
|---|---|
| `SCF` | Total SCF energy |
| `MP2` | MP2 correlation component: E_MP2 − E_SCF |
| `CCSD_T` | Triples component: E_CCSD(T) − E_MP2 |

This decomposition is defined in Step 3 of the PI task list (GitHub Issue #46).
The total BSSE at CCSD(T) level is the sum of contributions from all three
components; they are not independent estimates at different levels of theory.

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
Only fixed-y slices at y = 2.0, 2.25, 2.5, 2.75 Å have full 7×7 grid coverage
and are used for contour analysis. All other slices have incomplete grids and are
excluded.

**Rotational–radial configurations** — both the relative orientation and the
intermolecular separation are varied. Euler angles are sampled randomly and the
O–O separation R is varied independently. Of 64 unique orientations sampled,
32 have at least 5 distinct R values and are retained for analysis.

### Configuration labels

Translational: `n_x_y_z` — index `n`, O–O displacement `(x, y, z)` in Ångströms.

Rotational–radial: `n_α_β_γ_R` — index `n`, Euler angles `(α, β, γ)` in degrees,
O–O separation `R` in Ångströms.

### Key finding

BSSE does not scale purely with interaction energy in the dimer. The AC pair
(trimer labeling) shows orientation-dependent behavior including a
repulsive→attractive sign flip at R = 3.0 Å.

See `2/README.md` for the complete folder structure, scripts, and plot descriptions.

---

## Cluster 3 — H₂O Trimer (`3/`)

The trimer dataset samples 5 O–O distances for an equilateral triangle geometry,
covering the range R = 2.75–3.75 Å in steps of 0.25 Å. Each configuration
requires 19 unique NWChem calculations for the 3-body VMFC correction.

### Configuration labels

Configurations are labeled `n_distance`, where `n` is the configuration index
(0–4) and `distance` is the O–O separation in Ångströms.

| Configuration | O–O Distance (Å) |
|---|---|
| 0 | 2.75 |
| 1 | 3.00 |
| 2 | 3.25 |
| 3 | 3.50 |
| 4 | 3.75 |

### Key finding

BSSE does not scale purely with interaction energy in the trimer. The AC pair
shows an orientation-dependent repulsive→attractive sign flip at R = 3.0 Å,
consistent with the dimer result. The 2-body BSSE dominates; 3-body BSSE is
approximately one order of magnitude smaller.

See `3/README.md` for the complete folder structure, scripts, and plot descriptions.

---

## Cluster 4 — H₂O Tetramer (`4/`)

The tetramer dataset uses a single geometry: the optimal water tetramer obtained
from a geometry optimization. The focus is on decomposing the total BSSE into
its 2-body, 3-body, and 4-body contributions using two complementary formulations
(FORM(I) and FORM(II)). This requires 65 unique NWChem calculations.

### BSSE formulations

**FORM(I)** — body-ordered grouping: BSSE = Σ 2b_bsse + Σ 3b_bsse + 4b_bsse.
Compact and efficient; does not expose the internal sign structure.

**FORM(II)** — explicit basis-set-difference expansion: sums 50 individual
`[E_X(Y) − E_X(X)]` terms organized into 6 sign groups by (|X|, |Y|). The two
formulations are algebraically equivalent; their numerical equality was verified
for all three methods and serves as the primary consistency check.

### Key findings

- 2-body BSSE dominates: `|Σ 2b_bsse| ≈ 6.4 × 10⁻³ Ha`, approximately 80× larger
  than `|Σ 3b_bsse|`.
- The AD pair has positive 2b_bsse (repulsive pair geometry), mirroring the
  orientation-dependent behavior seen at smaller cluster sizes.
- The body-order hierarchy `|2b_bsse| >> |3b_bsse| ≈ |4b_bsse|` is explained
  by the alternating sign structure of FORM(II): the 2b level is the only level
  with no internal cancellation.
- The three method components (SCF, MP2, CCSD(T)) are additive parts of the total
  BSSE, not independent estimates. Correlation components account for ~90% of the
  total 2-body BSSE.

See `4/README.md` for the complete folder structure, calculation count, numerical
results, and plot descriptions.

---

## Auxiliary Scripts

| Script | Description |
|---|---|
| `check_convergence_4.py` | Verifies SCF, MP2, CCSD(T) convergence in all tetramer output files |
| `nwchem_to_xyz.py` | Converts NWChem input files to clean XYZ format for PyMOL |

The downstream pipeline scripts (`build_raw_data_4.py`, `many_body_interactions_4.py`,
`four_body_bsse.py`) live in `bsse_db/data/` alongside their trimer counterparts.

---

## Current Status

| Cluster | Directory | Geometries | Status |
|---|---|---|---|
| Dimer | `2/` | Translational + rotational–radial scan | ✅ Complete |
| Trimer | `3/` | 5 O–O distances (equilateral triangle) | ✅ Complete |
| Tetramer | `4/` | Single optimal geometry | ✅ Complete |

---

## References

- Valiron & Mayer (1997) *Chem. Phys. Lett.* **275**, 46–55 — VMFC method
- Richard, Bakr & Sherrill (2018) *J. Chem. Theory Comput.* **14**, 2386 — PI's unified VMFC/MBE framework
- GitHub Issue #46 — PI's revised task list and file naming conventions