#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:06:35 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
from local_evans import local_evans as local_evans

def get_eigenvalues(p, choices, left, right, num_nodes=10):
    """
    Find eigenvalues by locating sign changes of the Evans function
    and refining roots with bisection.
    """

    s = {
        "rtol": 1e-8,
        "atol": 1e-8
    }

    lambda_vals = np.linspace(
        left,
        min(right, np.sqrt(p["g"] / p["L0"])),
        num_nodes
    )

    D = np.zeros_like(lambda_vals)

    for j, lam in enumerate(lambda_vals):
        D[j] = local_evans(lam, s, p)

    # Detect sign changes
    index = np.where(D[1:] * D[:-1] < 0)[0]

    choices = np.atleast_1d(choices)

    if len(index) < np.max(choices):
        raise ValueError(
            f"Only {len(index)} roots detected, "
            f"but requested choice {np.max(choices)}."
        )

    out = np.zeros(len(choices))

    for cnt, choice in enumerate(choices):

        idx = index[-choice]

        a = lambda_vals[idx]
        b = lambda_vals[idx + 1]

        fa = D[idx]

        c = 0.5 * (a + b)

        while (b - a) > 1e-8:

            c = 0.5 * (a + b)
            fc = local_evans(c, s, p)

            if np.sign(fc) == np.sign(fa):
                a = c
                fa = fc
            else:
                b = c

        out[cnt] = c

    return out