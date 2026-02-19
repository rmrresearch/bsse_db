#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 18:43:43 2025

@author: leothan
"""


from pathlib import Path
import json
import random
import numpy as np
#--
from RotationMatrix import RotationMatrix


class MolecularRotation:

    def __init__(self, molecule):
        self._atoms = molecule['atoms']        
        # --- load your atomic data once ---
        DATA_PATH = Path(__file__).parent / "data" / "element_data.json"
        elements_data = json.loads(DATA_PATH.read_text())
        self._E = elements_data["elements"]  # shorthand
        #--- Dictionary of atoms and its coordiantes -------------------
        atoms_dict = {}
        new_atoms_dict = {}
        rot_atoms_dict = {}
        i = 1
        for a in self._atoms:
            element = f"{a[0]}"
            coord = a[1]
            atoms_dict[f'atom_{i}']={"atom":element,"coordinate":coord}
            new_atoms_dict[f'atom_{i}']={"atom":element,"coordinate":coord}
            rot_atoms_dict[f'atom_{i}']={"atom":element,"coordinate":coord}
            i+=1
        
        self.atoms_dict = atoms_dict
        self.new_atoms_dict = new_atoms_dict
        self.rot_atoms_dict = rot_atoms_dict
    # --- help functions ----------------------------------------------
    def mass(self,sym):
        """Atomic mass in amu (same numeric value as g/mol)."""
        return self._E[sym]["mass"]["value"]
    #--
    def atom_coord(self,atom):
        """Atom coordinates relative to the origin of the xyz-axis."""
        return self.atoms["coordinate"]
    
    # def vdw(sym: str) -> float:
    #     return float(E[sym]["vdw_radii"]["value"])    
    def cm_coord(self):
        """
        Computation of the center of mass (cm) of the molecule      
        m = numpy array of masses of the atoms in the molecule
        pos = numpy array of positions of the atoms in the moleucle
        """
        atoms = self.atoms_dict.copy()
        #--
        m = np.array([self.mass(a['atom']) for a in atoms.values()])
        pos = np.array([a['coordinate'] for a in atoms.values()])
        #--
        cm_num = (m[:,None]*pos).sum(axis=0)
        cm_den = m.sum()
        cm = cm_num/cm_den
        #--       
        return cm
    #--
    def atom_coord_relative_to_cm(self,atom_coord):
        """
        Determines the atoms coordinates of the molecule
        relative to the center of mass
        """
        center_of_mass = self.cm_coord()
        #--
        new_coord = atom_coord - center_of_mass
        return new_coord        
    #--
    def new_atoms_coord(self):
        new_atoms_coord = self.new_atoms_dict
        rot_atoms_dict = self.rot_atoms_dict
        for a in new_atoms_coord.values():
            a['coordinate'] = self.atom_coord_relative_to_cm(a['coordinate'])
        return new_atoms_coord, rot_atoms_dict
    #--
    
        
        
        
        
        
        
    
