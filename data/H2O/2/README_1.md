Going through the data and looking for how to represent it in a way that:
1. Follows the instructed file structure
2. Can be plotted for comparison and analysis

This is what I have collected in the Readme file (This is to avoid further confusion from my part)

# H₂O Dimer: Geometry and Directory Structure

## 1. Configuration and Geometry Parameters

### **Configuration (general idea)**

A molecular configuration is defined by the position of the center of mass (COM) of a molecule by  
its $(x, y, z)$ coordinates, and its relative orientation expressed through the Euler angles  
$(\alpha, \beta, \gamma)$. Thus a configuration vector is:

$\mathbf{s} = (x, y, z, \alpha, \beta, \gamma)$


# Configuration Used in the Original H₂O Dimer Dataset

In the original work (performed by a previous researcher), the H₂O dimer geometry was studied using two different configuration parameter sets:

1. **Translational configurations:**  
  $r = (x, y, z)$ with fixed relative Euler angles.

2. **Rotational–radial configurations:**  
   $s = (\alpha, \beta, \gamma, R)$,  
   where $R$ is the oxygen–oxygen distance between the monomers.

These two datasets explore different physical aspects of the dimer interaction:  
Dataset (1) varies **relative translation**, whereas Dataset (2) varies **relative orientation and intermolecular distance**.

---

## 1. Translational Configuration Dataset

The **relative position** of monomer 2 with respect to monomer 1 is represented by:

$\mathbf{r} = (x, y, z)$

with

$\mathbf{r} = \mathbf{r}_2 - \mathbf{r}_1 = (x_2 - x_1, y_2 - y_1, z_2 - z_1).$

and $r_1, r_2$ corresponds to the coordinates of the oxygen atom in monomer-1 and monomer-2 respectively

The coordinates sampled in the original dataset were:

- **x-coordinate values:**  
  {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75}

- **y-coordinate values:**  
  {1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75}

- **z-coordinate values:**  
  {0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75}

The number of (y,z) points sampled for each fixed x-value was:

| x-value | # of points | Meaning |
|:-------:|:-----------:|---------|
| 0.25    |     38      | 38 (y,z) grid points |
| 0.50    |     38      | 38 (y,z) grid points |
| 0.75    |     40      | 40 (y,z) grid points |
| 1.00    |     47      | 47 (y,z) grid points |
| 1.25    |     53      | 53 (y,z) grid points |
| 1.50    |     60      | 60 (y,z) grid points |
| 1.75    |     69      | 69 (y,z) grid points |


Then, each dataset corresponds to sampling the rectangle

$[1.00, 2.75] \times [0.25, 1.75] \subset \mathbb{R}^2$

in the yz-plane, with different densities for each x.

---

## 2. Rotational–Tanslational Configuration Dataset

The second dataset uses the configuration vector:

$$
s = (\alpha, \beta, \gamma, R)
$$

where the Euler angles describe the orientation of monomer 2 relative to monomer 1, and $R$ is the O–O separation distance.

### Euler angles
$$
\alpha, \beta, \gamma \in \{0^\circ, 90^\circ, 180^\circ, 270^\circ\}
$$

### Intermolecular distance

$R \in [1.61,3.91]$


Values of $R$ that produced atomic overlap were excluded to avoid unphysical configurations.

---

## Proposed Analysis

For fixed $\alpha$ and $\beta$, we analyze how the BSSE varies as:

- $\gamma$ changes,  
- $R$ changes.

This yields contour plots of BSSE as a function of $(\gamma, R)$ at each $(\alpha, \beta)$, enabling:

- comparison across orientations,
- identification of trends,
- characterization of anisotropy in the BSSE landscape.



---

# 3. Goal: Organizing the Data Systematically

To comply with the PI's organizational requirements, we map the original dataset into a structured
directory layout, based on the two analyzed configurations:
1. $r = (x, y, z)$ (Translational configuration data set)
2. $s = (\alpha, \beta, \gamma, R)$ (Rotational–Radial configuration data set)

To account for these configurations and organize the dimer data in a systematic way the following folder structure has been proposed:

`chemicalspecies/cluster/configuration_type/configuration/theorylevel/basis`

Where:
 - **chemicalspecies**: Molecule or Atom to be considered (e.g., Ne, H2O, HF)
 - **cluster**: cluster size, e.g., 2, 3, or 4
 - **configuration_type**: one of the two above (`tranlational` or `rotation_radial`)
 - **configuration**: number of the configuration $1,2,\dots,n$
 - **theory level**: Highest level of theory used (in our case CCSD(T), represented as `CCSD_T`)
 - **basis**: `aug-cc-pvdz`, `aug-cc-pvtz`, and possibly `aug-cc-pvqz`

Each configuration is assigned a unique index `n`, representing its position within the sequence of
all sampled configurations. The folder naming scheme is:

`<config-number>_<config>`

Where:
- `<config-number>` is an integer starting from **0**, uniquely identifying the configuration.
- `<config>` is either the position vector $r = (x, y, z)$ or the rotation–distance vector
  $s = (\alpha, \beta, \gamma, R)$.

---

## 3.1 Translational Configuration: $r = (x, y, z)$

For dimers defined by the translational configuration $r = (x, y, z)$, the folder containing all
associated data (input files, output files, and BSSE calculations) follows the pattern:

`n_x_y_z`

### **Example**

If:
- $r = (0.25, 0.25, 1.50)$  
- The configuration index is $n = 0$

Then the folder name is:

`0_0.25_0.25_1.5`


And the full path to this folder, when using the aug-cc-pvdz basis set, is:

`H2O/2/translational/0_0.25_0.25_1.5/CCSD_T/aug-cc-pvdz/`

---

## 3.2 Rotational–Radial Configuration: $s = (\alpha, \beta, \gamma, R)$

For dimers defined by rotational and radial parameters  
$s = (\alpha, \beta, \gamma, R)$, the folder containing all associated data follows:

`n_alpha_beta_gamma_R`

### **Example**

If:
- $s = (0,0,30,2.25)$  
- Configuration index is $n = 20$

Then the folder name becomes:

`20_0_0_30_2.25`

And the corresponding directory path is:

`H2O/2/rotational_radial/20_0_0_30_2.25/CCSD_T/aug-cc-pvdz/`



---

This naming convention ensures:

- Full traceability between configuration parameters and folder names  
- Compatibility with automated parsing and batch BSSE analysis  
- Clear separation between translational and rotational–radial datasets  


 







