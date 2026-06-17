#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:11:33 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
import matplotlib.pyplot as plt
from get_eigenvalues import get_eigenvalues as get_eigenvalues


print("Plot of the second largest eigenvalue as a function of k")

# Load previously computed largest eigenvalue branch
rts1 = np.load("rts1.npy")

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

choices = 2  # second largest eigenvalue
k_max = 685

left = 0.001
right = 2.0

k_vals = np.arange(1, k_max + 1) / p["L"]
rts2 = np.zeros_like(k_vals, dtype=float)

# Real-time plot (equivalent to MATLAB's plotting inside the loop)
plt.figure()

plt.plot(
    k_vals,
    rts1[:-1],
    '.r',
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

for j in range(len(k_vals) - 1):

    print(j + 1)

    p["k"] = k_vals[j]

    rt = get_eigenvalues(
        p,
        choices,
        left,
        right
    )

    rts2[j] = np.atleast_1d(rt)[0]

    left = rts2[j]

    if j > 1:
        right = min(
            2 * rts2[j] - rts2[j - 1],
            rts1[j + 1] - 1e-6
        )
    else:
        choices = 1
        right = rts1[1] - 1e-6

    plt.plot(
        p["k"],
        rts2[j],
        '.k',
        markersize=18
    )

plt.show()

# Save results
np.save("rts2.npy", rts2)

# Final figure
plt.figure()

plt.plot(
    k_vals,
    rts2,
    '.k',
    markersize=8
)

plt.xlabel('k', fontsize=18)
plt.ylabel(r'$\lambda_2(k)$', fontsize=18)

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
    1.1 * np.max(rts2)
])

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()