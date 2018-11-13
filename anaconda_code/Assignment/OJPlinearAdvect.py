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
    nx = 60
    nt = 60
    c = 0.4
        
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
    
    
    l2FTCS, errorFTCS = l2ErrorNorm(phiFTCS, phiAnalytic)
    l2FTBS, errorFTBS = l2ErrorNorm(phiFTBS, phiAnalytic)
    l2CTCS, errorCTCS = l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic)
    l2LW, errorLW = l2ErrorNorm(phiLW, phiAnalytic)
    
#    
    
    
    ##plot for FTCS
    font = {'size'   : 20}
    plt.rc('font', **font)
    plt.figure(1,figsize=(10,7))
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
             linestyle='--', linewidth=2)
    plt.plot(x, phiFTBS, label='FTBS', color='red')
    plt.plot(x, phiCTCS[nt-1,:], label='CTCS', color='green') #using second to last time step of t to plot
    plt.plot(x, phiLW, label='Lax-Wendroff', color="orange")  #using second to last time step to plot
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-0.2,1.4])  #increased y limiy to show where LW seems to be going wrong
    plt.legend()
    plt.xlabel('$x$')
    
    print("FTBS l2 error norm = ", l2FTBS)
    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
    
     
    print("CTCS l2 error norm = ", l2CTCS)
    print("CSCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
    
     
    print("LW l2 error norm = ", l2LW)
    print("LW linf error norm = ", lInfErrorNorm(phiLW, phiAnalytic))
    ##error plot
    plt.figure(5)
    plt.clf()
    plt.ion()
    plt.plot(x, errorFTCS)
    plt.title('error plot of FTCS method')
    
    
    
#    ##plot for FTBS
#    plt.figure(2,figsize=(10,7))
#    plt.clf()
#    plt.ion()
#    plt.plot(x, phiOld, label='Initial', color='black')
#    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
#             linestyle='--', linewidth=2)
#    plt.plot(x, phiFTBS, label='FTBS', color='red')
#    plt.axhline(0, linestyle=':', color='black')
#    plt.ylim([-0.2,1.2])  #increased y 
#    plt.legend()
#    plt.xlabel('$x$')
#    print("FTBS l2 error norm = ", l2FTBS)
#    print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
#    
#    
#    ##plot for CTCS
#    plt.figure(3,figsize=(10,7))
#    plt.clf()
#    plt.ion()
#    plt.plot(x, phiOld, label='Initial', color='black')
#    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
#             linestyle='--', linewidth=2)
#    plt.plot(x, phiCTCS[nt-1,:], label='CTCS', color='green') #using second to last time step of t to plot
#    plt.axhline(0, linestyle=':', color='black')
#    plt.ylim([-0.2,1.2])  #increased y limiy to show where LW seems to be going wrong
#    plt.legend()
#    plt.xlabel('$x$')
#    print("CTCS l2 error norm = ", l2CTCS)
#    print("CTCS linf error norm = ", lInfErrorNorm(phiCTCS[nt-1,:], phiAnalytic))
#    
#    ##plot for LW
#    plt.figure(4,figsize=(10,7))
#    plt.clf()
#    plt.ion()
#    plt.plot(x, phiOld, label='Initial', color='black')
#    plt.plot(x, phiAnalytic, label='Analytic', color='black', 
#             linestyle='--', linewidth=2)
#    plt.plot(x, phiLW, label='Lax-Wendroff', color="orange")  #using second to last time step to plot
#    plt.axhline(0, linestyle=':', color='black')
#    plt.ylim([-0.2,1.2])  #increased y limiy to show where LW seems to be going wrong
#    plt.legend()
#    plt.xlabel('$x$')
#    
#    print("Lax-Wendroff l2 error norm = ", l2LW)
#    print("Lax-Wendroff linf error norm = ", lInfErrorNorm(phiLW, phiAnalytic))
#    
##    plt.figure(2)
##    plt.plot(nt, errorFTBS)
#    #input('press return to save file and continue')
#    #plt.savefig('plots/mixed_different_coeff_2_initial_conditions.pdf')
#            
####Run the function main defined in this file                      ###
            
main()

def convergence_exp():
    c = 0.4  #courant number constant
    xmin = 0
    xmax = 1
    l2FTBS_dx_err = np.zeros(10)
    l2CTCS_dx_err = np.zeros(10)
    l2LW_dx_err = np.zeros(10)
    errorFTBS = np.zeros(10)
    errorCTCS = np.zeros(10)
    errorLW = np.zeros(10)
    dx_it = np.zeros(10)
    
    for i in range(1,10):
        print(i)
        nx = i*10 
        dx = (xmax - xmin)/nx
        nt = nx ##keeping overall time constant
        dx_it[i] = dx
        
        # spatial points for plotting and for defining initial conditions
        x = np.arange(xmin, xmax, dx)
        
        # Initial conditions
        phiOld = cosBell(x, 0.25, 0.75)
        # Exact solution is the initial condition shifted around the domain
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.4, 0.6)#, 0.5, 0.75)
        
        # Advect the profile using finite difference for all the time steps
        phiFTBS = FTBS(phiOld, c, nt)
        phiCTCS = CTCS(phiOld, c, nt, nx)
        phiLW = LW(phiOld, c, nt, nx)
        
        
        l2FTBS_dx_err[i], errorFTBS = l2ErrorNorm(phiFTBS, phiAnalytic)
        l2CTCS_dx_err[i], errorCTCS = l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic)
        l2LW_dx_err[i], errorLW = l2ErrorNorm(phiLW, phiAnalytic)
    
    dxsquared = dx**2
    plt.figure()
    plt.loglog(dx_it, l2FTBS_dx_err, label='FTBS', color = 'red')
    plt.loglog(dx_it, l2CTCS_dx_err, label='CTCS', color = 'green')
    plt.loglog(dx_it, l2LW_dx_err, label = 'LW', color = 'orange')
#    plt.loglog(dx, dxsquared, label = 'dx^2', linestyle=':', color='black')
    plt.legend()
    
    
convergence_exp()
    
        
        
        
        
    
        
convergence_exp()
c_list = (-0.4, 0.1, 0.4, 0.9, 1, 1.1)
print(c_list)
def c_exp():
    
    # Parameters
    xmin = 0
    xmax = 1
    nx = 60
    nt = 60
        
    # Derived parameters
    dx = (xmax - xmin)/nx
    
    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)
    
    for i in range(len(c_list)):
        
        c = c_list[i]
    
        # Initial conditions
        phiOld = cosBell(x, 0.25, 0.75)
        # Exact solution is the initial condition shifted around the domain
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.25, 0.75)
        
        # Advect the profile using finite difference for all the time steps
        phiFTCS = FTCS(phiOld, c, nt)
        phiFTBS = FTBS(phiOld, c, nt)
        phiCTCS = CTCS(phiOld, c, nt, nx)
        phiLW = LW(phiOld, c, nt, nx)
        
        
        l2FTCS, errorFTCS = l2ErrorNorm(phiFTCS, phiAnalytic)
        l2FTBS, errorFTBS = l2ErrorNorm(phiFTBS, phiAnalytic)
        l2CTCS, errorCTCS = l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic)
        l2LW, errorLW = l2ErrorNorm(phiLW, phiAnalytic)
        
        font = {'size'   : 20}
        plt.rc('font', **font)
        plt.figure(figsize=(10,7))
        plt.clf()
        plt.ion()
        plt.plot(x, phiOld, label='Initial', color='black')
        plt.plot(x, phiAnalytic, label='Analytic', color='black', 
                 linestyle='--', linewidth=2)
        plt.plot(x, phiFTBS, label='FTBS', color='red')
        plt.plot(x, phiCTCS[nt-1,:], label='CTCS', color='green') #using second to last time step of t to plot
        plt.plot(x, phiLW, label='Lax-Wendroff', color="orange")  #using second to last time step to plot
        plt.axhline(0, linestyle=':', color='black')
        plt.ylim([-0.2,1.4])  #increased y limiy to show where LW seems to be going wrong
        plt.legend()
        plt.xlabel('$x$')
        plt.title('Linear Advection where c=%f'%c)
        
        print("FTBS l2 error norm = ", l2FTBS)
        print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
        
        
        print("CTCS l2 error norm = ", l2CTCS)
        print("CSCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
        
        
        print("LW l2 error norm = ", l2LW)
        print("LW linf error norm = ", lInfErrorNorm(phiLW, phiAnalytic))
       
    
        
c_exp()
        

        
        
        
