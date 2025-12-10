# README_2.md  
## H2O/2 — Folder Contents and Description

This directory contains all computational and structural data for the **H₂O–H₂O dimer study**, including translational scans, rotational–radial scans, and the corresponding 3D coordinate files used for visualization and analysis.

Below is a description of each subfolder and file contained in `H2O/2`.

---

## 📁 **Folders**

### **1. `translational/`**  
Contains computations in which **monomer 2** (specifically the second oxygen atom) is translated with respect to **monomer 1**.  
These configurations represent **center-of-mass translations** without changing orientation.

---

### **2. `rotational_radial/`**  
Contains computations where **monomer 2 is rotated and placed at a radial distance R from monomer 1**.  
The relative orientation is defined using **Euler angles**:

$\alpha,\ \beta,\, \gamma$


Angles are stored and interpreted in **degree units**.  
These configurations form the dataset used in the rotational–radial study of the H₂O dimer.

---

### **3. `H2O_dimers_xyz/`**  
Contains all **XYZ coordinate files** used for the translational study.  
These files can be loaded into PyMOL for visualization of each translational geometry.

---

### **4. `H2O_dimers_euler_angles_xyz/`**  
Contains all coordinate files describing configurations characterized by Euler angles:

$(\alpha, \beta, \gamma, R)$

These are the input geometries for the rotational–radial study and can also be visualized in PyMOL.

---

## 📄 **Files**

### **`README_1`**  
A detailed description of the entire folder structure and the full methodology behind the H₂O dimer study.

---

### **`README_2`**  
*This file* — a concise summary of the contents of the `H2O/2` directory.

---

### **`load_h2o_dimers.py`**  
A Python script used to load all **translational XYZ coordinate files** into PyMOL.  
Inside PyMOL, run:

`run load_h2o_dimers.py`

to automatically load and visualize all geometries.

---

### **`load_h2o_dimers_euler_angles.py`**  
A Python script used to load all **Euler-angle XYZ coordinate files** into PyMOL.  
Inside PyMOL, run:

`run load_h2o_dimers_euler_angles.py`

to visualize all rotational–radial configurations.

---