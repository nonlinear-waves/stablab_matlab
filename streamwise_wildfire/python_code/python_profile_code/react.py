#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:07:29 2024

@author: oliviachandrasekhar
"""
import numpy as np

from python_profile_code.parameters import vstar, l

#reaction term


#reaction term
def react(u):
    out=np.zeros(u.shape)
    for j in range(0, len(u)):
        if u[j] > 0:
            out[j]=np.exp(-1/u[j])
    return out

#derivative of the reaction term 
def react_prime(u):
    if (u <= 0):
        return -l
    else:
        return (vstar/u**2)*np.exp(-1/u)-l
    
