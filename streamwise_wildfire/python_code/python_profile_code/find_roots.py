#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 14:14:19 2024

@author: oliviachandrasekhar
"""

import numpy as np
import scipy as scipy

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import CubicSpline

from python_profile_code.ode_theta import ODE_theta

from python_profile_code.react import react_prime

from python_profile_code.parameters import alpha, vstar, l

#from find_roots import find_roots



from python_profile_code.react import react

#spacing algoritm to cluster around either end--cluster towards 0 if spacing >1 
def my_lin(lb, ub, steps, spacing=5):
    span = (ub-lb)
    dx = 1.0 / (steps-1)
    return [lb + (i*dx)**spacing*span for i in range(steps)]



def find_roots(vstar, l):
    
    #fixed points are given by roots of this function
    def func(u):
        return vstar*react(u)-l*u
    

    lin = my_lin(0, 10, 11)
    umax=37  #upper bound based on largest zero of func(u) when v_val=1 
    urange=my_lin(0,umax, 10)

    zero_init=0 #initial guess
    all_zeros=[] #make sure list contains that initial guess-specific to this problem because I know 0 is always a zero
    tol=10 #decimal place tolerance--could probably do something where I vary this 

    #find all zeros--number of zeros will vary on number of entries in urange (fineness of search)
    for x0 in urange:
        zero = scipy.optimize.fsolve(func, x0) 
        zero_num=zero[0]
        all_zeros.append(zero_num)
    
        #print("I found all of these zeros:  " + str(all_zeros))

    distinct_zeros=[]

    for i in range(len(all_zeros)-1): #-1 to prevent "out of range error"
        if (round(all_zeros[i],tol)!=round(all_zeros[i+1],tol) and round(all_zeros[i],tol) not in distinct_zeros) :
            distinct_zeros.append(all_zeros[i])
    distinct_zeros.append(all_zeros[len(all_zeros)-1])
        
    
    umin=distinct_zeros[0]
    
    if len(distinct_zeros)<2:
        umax=distinct_zeros[0]
    else:
        umax=distinct_zeros[2]
    
    return umin, umax