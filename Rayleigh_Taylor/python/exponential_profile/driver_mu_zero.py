#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Mon Jun 15 11:26:31 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
from scipy.integrate import solve_ivp

# Parameters
p = {}
p["rho_m"] = 0.1
p["rho_p"] = 1.0  # bigger than rho_m
p["h"] = 1.5
p["L"] = 1.0

# Don't change
p["g"] = 9.8
p["a"] = np.log(p["rho_p"] / p["rho_m"]) / p["h"]

# ODE system
def ode(x, y):
    A = np.array([
        [0, 1],
        [-(p["a"]**2) / 4, -p["a"]]
    ])
    return A @ y

# Solve from x = -h to x = 0
sol = solve_ivp(
    ode,
    (-p["h"], 0.0),
    [0.0, 1.0],
    method="BDF",      # similar to MATLAB ode15s
    rtol=1e-10,
    atol=1e-10,
    dense_output=True
)

# Evaluate solution at x = 0
temp = sol.sol(0.0)

constant = -p["g"] * temp[1] / temp[0]
print("constant =", constant)

k2 = -p["a"] * temp[1] / temp[0] - p["a"]**2 / 4
print("k2 =", k2)

term = (2 / p["h"]) * (2 * p["h"] - p["a"]**2) - 1 / (4 * p["h"]**2)
print("term =", term)