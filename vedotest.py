#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 09:27:21 2025

@author: leothan
"""
import numpy as np
import 

class PeriodicTable:
    """
    A Vedo-compatible class for accessing periodic table data, wrapping vtkPeriodicTable.
 
    This class provides access to element properties such as atomic numbers, names,
    symbols, and covalent radii, using VTK's built-in periodic table database.
 
    Attributes:
        periodic_table (vtkPeriodicTable): The underlying VTK periodic table object.
    """
 
    def __init__(self):
        """
        Initialize the PeriodicTable with VTK's built-in periodic table data.
        """
        self.periodic_table = vtkPeriodicTable()
    
    def get_atomic_number(self, symbol):
        """
        Get the atomic number of the element with the given symbol.
    
        Arguments:
            symbol : (str)
                The symbol of the element.
    
        Returns:
            The atomic number of the element.
        """
        return self.periodic_table.GetAtomicNumber(symbol)

        
    def get_vdw_radius(self, atomic_number):
        """
        Get the Van der Waals radius of the element with the given atomic number.

        Arguments:
            atomic_number: (int)
                The atomic number of the element.
        
        Returns:
            The Van der Waals radius of the element.
        """
        return self.periodic_table.GetVDWRadius(atomic_number)
        
        
        
        
        
        
        