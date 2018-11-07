# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

# The numpy package for numerical functions and pi
import numpy as np

def FTCS(phiOld, c, nt):
    "Linear advection of profile in phiOld using FTCS, Courant number c"
    "for nt time-steps"
    
    nx = len(phiOld)
    # new time-step array for phi
    phi = phiOld.copy()

    # FTCS for each time-step
    for it in range(nt):
        # Loop through all space using remainder after division (%)
        # to cope with periodic boundary conditions
        for j in range(nx):
            phi[j] = phiOld[j] - 0.5*c*\
                     (phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
        
        # update arrays for next time-step
        phiOld = phi.copy()
        

    return phi

#trying FTBS
    
def FTBS(phiOld, c, nt):
    nx=len(phiOld)
    
    phi = phiOld.copy()
    
    for it in range(nt):
        for j in range(nx):
            phi[j] = phiOld[j] - c*(phiOld[j]%nx - phiOld[(j-1)%nx])
            
            
        phiOld = phi.copy()
            
    return phi

    


        

def CTCS(phiOld, c, nt, nx):
    phi = np.zeros((nt,nx))  #making matrix for n time steps and j spacial steps
    
    phi[0,:] = phiOld.copy()  #initial conditions for time step t=0
    phi[1,:] = FTCS(phi[0:1], c, nt)  #initial conditions for t=1 using FTCS  

    
    for n in range(1, nt-1):
        for j in range(1, nx-1):
            phi[n+1,j] = phi[n-1,j] - c*(phi[n,j+1] - phi[n, j-1])
            
    return phi
    



def LW(phiOld, c, nt, nx):
    phi = phiOld.copy()
    phijplus = phiOld.copy()
    phijminus = phiOld.copy()
   
    for n in range(nt):
        
        
        for j in range(nx):
            phijplus[j] = 0.5*(1+c)*phiOld[j] +0.5*(1-c)*phiOld[(j+1)%nx]
            phijminus[j]= 0.5*(1+c)*phiOld[(j-1)%nx] +0.5*(1-c)*phiOld[j]
            phi[j]=phiOld[j]-c*(phijplus[j]-phijminus[j])
              
       
        phiOld = phi.copy()
    
    return phi
    
