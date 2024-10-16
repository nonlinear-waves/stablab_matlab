import numpy as np

from python_profile_code.parameters import alpha, vstar, l

from python_profile_code.react import react_prime

from python_profile_code.find_roots import find_roots

from scipy.integrate import odeint

from python_profile_code.ode_theta import ODE_theta

umin,umax=find_roots(vstar,l)



def get_solutions(omega,c,x_dom,profile_poly):
   
    #fixed points in angular system
    
        fp_LH_1=0.5*(-(c-alpha)+np.sqrt((c-alpha)**2+4*(-omega-react_prime( umin ))))
        fp_LH_2=0.5*(-(c-alpha)-np.sqrt((c-alpha)**2+4*(-omega-react_prime( umin ))))

        theta_LH_1=np.arctan(fp_LH_1)
        theta_LH_2=np.arctan(fp_LH_2)
    
 
    
        if np.sign(fp_LH_1)==1.0:
            theta_0= theta_LH_1
        elif np.sign(fp_LH_2)==1.0:
            theta_0= theta_LH_2 
        # else: print ("for omega= "+ str(omega)+" fp at neg infinity not a saddle!")
        
    
        fp_RH_1=0.5*(-(c+alpha)+np.sqrt((c+alpha)**2+4*(-omega-react_prime( umax ))))
        fp_RH_2=0.5*(-(c+alpha)-np.sqrt((c+alpha)**2+4*(-omega-react_prime( umax ))))

        theta_RH_1=np.arctan(fp_RH_1)
        theta_RH_2=np.arctan(fp_RH_2)
    

        if np.sign(fp_RH_1)==-1.0:
            theta_f= theta_RH_1
            theta_other=theta_RH_2
        elif np.sign(fp_RH_2)==-1.0:
            theta_f= theta_RH_2 
            theta_other=theta_RH_1
    
        sol = odeint(ODE_theta, theta_0, x_dom,args=(c,omega, profile_poly), rtol=10**(-8))
        
        
        return sol
 