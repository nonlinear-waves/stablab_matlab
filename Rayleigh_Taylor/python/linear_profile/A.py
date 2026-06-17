#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:01:14 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np

def A(x, lam, p):
    """
    Translation of MATLAB function A(x, lambda, ~, p)

    Parameters
    ----------
    x : float
    lam : float
        Eigenvalue parameter (lambda is a reserved keyword in Python).
    p : dict or object
        Must contain:
        rho_m, rho_0_der, h, k, g

    Returns
    -------
    numpy.ndarray
        2x2 matrix
    """
    rho_0 = p["rho_m"] + p["rho_0_der"] * (x + p["h"])

    fun = p["rho_0_der"] / rho_0

    out = np.array([
        [0.0, 1.0],
        [p["k"]**2 * (1.0 - p["g"] * fun / lam**2), -fun]
    ])

    return out