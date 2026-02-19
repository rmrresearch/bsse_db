#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 18:43:43 2025

@author: leothan
"""

from pathlib import Path
import json
import numpy as np


class AtomicInfo:
    
    def __init__(self,molecule):
        self.molecule = molecule
    
    def get_atoms(self):
        atoms = self.molecule['atoms']
        return atoms
    
    def get_atomic_mass_coord(self):
        # --- load your atomic data once ---
        DATA_PATH = Path(__file__).parent / "data" / "element_data.json"
        elements_data = json.loads(DATA_PATH.read_text())
        E = elements_data["elements"]  # shorthand
        #--- Atomic Masses ----------------------------------------------
        atoms = self.get_atoms()
        for a in atoms:
            mass_of_a = a[0]

# # --- helpers you can call anywhere ---
# def mass(sym: str) -> float:
#     """Atomic mass in amu (same numeric value as g/mol)."""
#     return float(E[sym]["mass"]["value"])

# def vdw(sym: str) -> float:
#     """van der Waals radius in Å from your JSON."""
#     return float(E[sym]["vdw_radii"]["value"])