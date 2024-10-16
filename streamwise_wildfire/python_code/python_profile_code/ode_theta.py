#define the ODE for the eigenvalue problem (autonomous)

import numpy as np 
#import os
#dir=os.getcwd()

from python_profile_code.react import react_prime
from python_profile_code.wind import wind

from python_profile_code.parameters import alpha

def ODE_theta(y,x,c,omega, profile_poly):
    theta=y 
    
    u=profile_poly(x)
    
    return (-react_prime(u)+omega)*(np.cos(theta))**2-(c-(wind(x,alpha)))*np.cos(theta)*np.sin(theta)-(np.sin(theta))**2