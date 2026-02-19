#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 12:08:59 2025

@author: leothan
"""
import numpy as np
import random

class RotationMatrix:
    """
    Right-handed active rotations. By default uses intrinsic Y-X-Z order:
    Rotation_Matrix = Ry * Rx * Rz
    """
    
    def __init__(self, theta_x=None, theta_y=None, theta_z=None):
        if theta_x is not None:
            self.theta_x = float(theta_x) 
        else:
            random.uniform(0, 2*np.pi)
        if theta_y is not None:
            self.theta_y = float(theta_y) 
        else:
            random.uniform(0, 2*np.pi)
        if theta_z is not None:
            self.theta_z = float(theta_z) 
        else:
            random.uniform(0, 2*np.pi)

    #--  
    def trig_funct_eval(self,theta):
        c, s = np.cos(theta), np.sin(theta)
        return c, s
    #--
    def x_rotation_matrix(self):
        theta_x = self.theta_x
        c, s = self.trig_funct_eval(theta_x)
        rot_x = np.array([
            [1, 0, 0],
            [0, c,-s],
            [0, s, c],
            ])
        return rot_x
    #----------------------------------------------------------------------
    def y_rotation_matrix(self):
        theta_y = self.theta_y
        c, s = self.trig_funct_eval(theta_y)
        rot_y = np.array([
            [c, 0, s],
            [0, 1, 0],
            [-s,0, c],
            ])
        return rot_y
    #----------------------------------------------------------------------
    def z_rotation_matrix(self):
        theta_z = self.theta_z
        c, s = self.trig_funct_eval(theta_z)
        rot_z = np.array([
            [c,-s, 0],
            [s, c, 0],
            [0, 0, 1],
            ])
        return rot_z
    #----
    def rotation_matrix(self):
        rot_z = self.z_rotation_matrix()
        rot_x = self.x_rotation_matrix()
        rot_y = self.y_rotation_matrix()
        Rotation_Matrix = rot_y @ rot_x @ rot_z
        return Rotation_Matrix