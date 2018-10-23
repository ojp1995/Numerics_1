# Numerical schemes for simulating linear advection for outer code
# linearAdvect.py 

# If you are using Python 2.7 rather than Python 3, import various
# functions from Python 3 such as to use real number division
# rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1
#from __future__ import absolute_import, division, print_function

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
#    

#def CTCS_1(phiOld, c , nt):
#    nx=len(phiOld)
#    
#    phi = phiOld.copy()
#    
#    for it in range(nt):
#        for n in range(1,nx-1):
#            phi[n] = phiOld[n] - (c/2)*(phiOld[n+1]%nx - phiOld[n-1]%nx)
#            
#        phiOld = phi.copy()
#        
#    return phi\
        

def CTCS(phiOld, c, nt, nx):
    phi = np.zeros((nt,nx))  #making matrix for n time steps and j spacial steps
    
    phi[0,:] = phiOld.copy()  #initial conditions for time step t=0
    phi[1,:] = phiOld.copy()  #initialconditions for time step t=1 
#    
#    for j in range(1,nx-1):
#        phi[1:j] = phi[0:j] - (0.5*c)*(phi[0:j+1]%nx - phi[0:(j-1)%nx])

    
    for n in range(1, nt-1):
        for j in range(1, nx-1):
            phi[n+1,j] = phi[n-1,j] - c*(phi[n,j+1] - phi[n, j-1])
            
    return phi
    
#
#def lax_wend(phiOld, c, nt, nx):
#    
#    phi = np.zeros((nx,nt))  #making matrix for n time steps and j spacial steps
#    
#    phi[0,:] = phiOld.copy()
#    
#    phijplus = np.zeros(nx)  #creating the half time step vectors
#    phijminus = np.zeros(nx)
#    
#    for n in range (0,nx-1):
#        for j in range (0,nx-1):
#            phijplus[j] = 0.5*(1+c)*phi[n,j]%nx + 0.5*(1-c)*phi[n,j+1]%nx  #for each spacial step computing new phi j+1/2, j-1/2
#            phijminus[j] = 0.5*(1+c)*phi[n,j-1]%nx + 0.5*(1-c)*phi[n,j]%nx
#            
#            phi[n+1,j] = -c*(phijplus[j]%nx-phijminus[j]%nx)+phi[n,j]%nx  #each row on the matrix computes for each timestep
#    print(phi[:,:])        
#    return phi


#def LW(phiOld, c, nt, nx):
#    phi = phiOld.copy()
#    phijplus = phiOld.copy()
#    phijminus = phiOld.copy()
#    for n in range (nt):
#        for j in range  (nx):
#            phijplus[j] = 0.5*(1+c)*phiOld[j] + 0.5*(1-c)*phi[(j+1)%nx]
#            phijminus[j] = 0.5*(1+c)*phiOld[(j-1)%nx] + 0.5*(1-c)*phiOld[j]
#            
#            phi[j] = -c*(phijplus[j]+phijminus[j]) + phiOld[j]
#            
#        phiOld = phi.copy()
#    
#        
#     
#  
#            
#    return phi
#    
    
    

    
def LW(phiOld, c, nt, nx):
    phi = phiOld.copy()
    phijplus=phiOld.copy()
    phijminus=phiOld.copy()
   
    for n in range(nt):
        
        
        for j in range(nx):
            phijplus[j] = 0.5*(1+c)*phiOld[j] +0.5*(1-c)*phiOld[(j+1)%nx]
            phijminus[j]= 0.5*(1+c)*phiOld[(j-1)%nx] +0.5*(1-c)*phiOld[j]
            phi[j]=phiOld[j]-c*(phijplus[j]-phijminus[j])
              
       
        phiOld = phi.copy()
    
    return phi
    
    
    