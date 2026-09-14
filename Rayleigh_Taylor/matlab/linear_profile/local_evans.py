#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 09:04:04 2026

Translated from Matlab to Python with ChatGPT

@author: blakebarker
"""

import numpy as np
from scipy.integrate import solve_ivp
from A import A as A

def local_evans(lam, s, p):
    """
    Translation of MATLAB local_evans(lambda,s,p)
    """

    def fun(x, y):
        return A(x, lam, p) @ y

    sol = solve_ivp(
        fun,
        (-p["h"], 0.0),
        [0.0, 1.0],
        method="BDF",  # stiff solver, similar role to ode23s/ode15s
        rtol=s["rtol"],
        atol=s["atol"]
    )

    temp = sol.y[:, -1]

    out = p["rho_p"] * (
        lam**2 * temp[1]
        + p["g"] * p["k"]**2 * temp[0]
    )

    return out