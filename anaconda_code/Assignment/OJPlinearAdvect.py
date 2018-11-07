#!/usr/bin/python3

# Outer code for setting up the linear advection problem on a uniform
# grid and calling the function to perform the linear advection and plot.

### The matplotlib package contains plotting functions              ###
import matplotlib.pyplot as plt
import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application
from OJPinitialConditions import *
from OJPadvectionSchemes import *
from OJPdiagnostics import *

### The main code is inside a function to avoid global variables    ###



 
def main():
    "Advect the initial conditions using various advection schemes and"
    "compare results"
        
    # Parameters
    xmin = 0
    xmax = 1
    nx = 40
    nt = 40
    c = 0.1
        
    # Derived parameters
    dx = (xmax - xmin)/nx
    
    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)
    
    # Initial conditions
    phiOld = cosBell(x, 0.25, 0.75)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.25, 0.75)#, 0.5, 0.75)
    
    # Advect the profile using finite difference for all the time steps
    phiFTCS = FTCS(phiOld, c, nt)
    phiFTBS = FTBS(phiOld, c, nt)
    phiCTCS = CTCS(phiOld, c, nt, nx)
    phiLW = LW(phiOld, c, nt, nx)
    
    
    
    # Calculate and print out error norms
    
    
    
   
    
    
    
    ##plot for FTCS
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1,figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTCS, label='FTCS', color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])  #increased y limiy to show where LW seems to be going wrong
    plt.legend()
    plt.xlabel('$x$')
    print("FTCS l2 error norm = ", l2ErrorNorm(phiFTCS, phiAnalytic))
    print("FTCS linf error norm = ", lInfErrorNorm(phiFTCS, phiAnalytic))
    
    
    
    ##plot for FTBS
    plt.figure(2,figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTBS, label='FTBS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])  #increased y 
    plt.legend()
    plt.xlabel('$x$')
    print("FTBS l2 error norm = ", l2ErrorNorm(phiFTBS, phiAnalytic))
    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
    
    
    ##plot for CTCS
    plt.figure(3,figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, phiCTCS[nt-1,:], label='CTCS', color='green') #using second to last time step of t to plot
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])  #increased y limiy to show where LW seems to be going wrong
    plt.legend()
    plt.xlabel('$x$')
    print("CTCS l2 error norm = ", l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic))
    print("CTCS linf error norm = ", lInfErrorNorm(phiCTCS[nt-1,:], phiAnalytic))
    
    ##plot for LW
    plt.figure(4,figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, phiLW, label='Lax-Wendroff', color="orange")  #using second to last time step to plot
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.2])  #increased y limiy to show where LW seems to be going wrong
    plt.legend()
    plt.xlabel('$x$')
    
    print("Lax-Wendroff l2 error norm = ", l2ErrorNorm(phiLW, phiAnalytic))
    print("Lax-Wendroff linf error norm = ", lInfErrorNorm(phiLW, phiAnalytic))
    
#    plt.figure(2)
#    plt.plot(nt, errorFTBS)
    #input('press return to save file and continue')
    #plt.savefig('plots/mixed_different_coeff_2_initial_conditions.pdf')
            
###Run the function main defined in this file                      ###
            
main()

        
