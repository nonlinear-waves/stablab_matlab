#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:12:37 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
import matplotlib.pyplot as plt
from get_eigenvalues import get_eigenvalues as get_eigenvalues


print("Plot of the largest 5 eigenvalues for k = 680")

# Parameters
p = {}
p["rho_m"] = 0.1
p["rho_p"] = 0.2  # bigger than rho_m
p["h"] = 2.0
p["L"] = 1.0

# Fixed
p["g"] = 9.8
p["rho_0_der"] = (p["rho_p"] - p["rho_m"]) / p["h"]

x = np.linspace(-p["h"], 0.0, 1000)
rho_0 = p["rho_m"] + p["rho_0_der"] * (x + p["h"])

p["L0"] = 1.0 / np.max(p["rho_0_der"] / rho_0)

# Compute first 5 eigenvalues
choices = np.arange(1, 6)  # MATLAB: 1:1:5

k_max = 680

left = 1.5
right = 2.193

p["k"] = k_max / p["L"]

rts = get_eigenvalues(
    p,
    choices,
    left,
    right,
    num_nodes=1000
)

# Plot
plt.figure()

bound = np.sqrt(p["g"] / p["L0"])

plt.plot(
    [1, 5],
    [bound, bound],
    '--b',
    linewidth=2
)

plt.plot(
    np.arange(1, len(rts) + 1),
    rts,
    '.-k',
    markersize=18,
    linewidth=2
)

plt.axis([
    1,
    5,
    0,
    1.1 * np.max(rts)
])

plt.xlabel('Mode index j', fontsize=18)
plt.ylabel(r'$\lambda$', fontsize=18)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()