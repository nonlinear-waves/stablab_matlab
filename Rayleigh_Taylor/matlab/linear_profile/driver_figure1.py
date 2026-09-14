#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:07:36 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
import matplotlib.pyplot as plt
from get_eigenvalues import get_eigenvalues as get_eigenvalues

print("Plot of the largest eigenvalue as a function of k")

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

choices = 1  # largest eigenvalue
k_max = 686

left = 0.01
right = 2.0

k_vals = np.arange(1, k_max + 1) / p["L"]
rts1 = np.zeros_like(k_vals, dtype=float)

for j in range(len(k_vals)):
    print(j + 1)

    p["k"] = k_vals[j]

    rts = get_eigenvalues(p, choices, left, right)

    # get_eigenvalues returns an array
    rts1[j] = np.atleast_1d(rts)[0]

    left = rts1[j]

    if j > 0:
        right = 2 * rts1[j] - rts1[j - 1]
    else:
        right = 2.0

# Save results (MATLAB save('rts1','rts1'))
np.save("rts1.npy", rts1)

# Plot
plt.figure()

plt.plot(
    k_vals,
    rts1,
    '.k',
    markersize=8
)

plt.xlabel('k', fontsize=18)
plt.ylabel(r'$\lambda_1(k)$', fontsize=18)

bound = np.sqrt(p["g"] / p["L0"])

plt.plot(
    [k_vals[0], k_vals[-1]],
    [bound, bound],
    '--b',
    linewidth=2
)

plt.axis([
    -2,
    k_vals[-1],
    0,
    1.1 * np.max(rts1)
])

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()