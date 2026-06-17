#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:25:12 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
import matplotlib.pyplot as plt
from get_eigenvalues import get_eigenvalues


print("Plot of the second largest eigenvalue as a function of k")

# Parameters
p = {}
p["rho_m"] = 0.1
p["rho_p"] = 1.0  # bigger than rho_m
p["h"] = 1.5
p["L"] = 1.0

# Don't change
p["g"] = 9.8
p["a"] = np.log(p["rho_p"] / p["rho_m"]) / p["h"]

profile_ratio = lambda x: p["a"]

# Match MATLAB exactly
Atwood_number = p["rho_p"] - p["rho_m"] / (p["rho_p"] + p["rho_m"])
print(f"\nAtwood number = {Atwood_number:.4g}")

bound = np.sqrt(p["a"] * p["g"])

# Compute eigenvalues
k = 200

rts = np.zeros(5)
indices = np.arange(1, 6)  # MATLAB: 1:1:5

for j in range(len(rts)):
    rts[j] = get_eigenvalues(indices[j], k, p["a"], p["h"])

# Plot
plt.figure()

plt.plot(
    [1, 5],
    [bound, bound],
    '--b',
    linewidth=2,
    label='Bound'
)

plt.plot(
    indices,
    rts,
    '.-k',
    markersize=18,
    linewidth=2
)

plt.axis([
    1,
    5,
    3.86,
    1.001 * np.max(rts)
])

plt.xlabel('Mode index j', fontsize=18)
plt.ylabel(r'$\lambda$', fontsize=18)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()