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
    u=0.4  ##wind speed, keeping constant
    c = u*(nx/nt)
        
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
    phiCTCS = CTCS(phiOld, c, nt)
    phiLW = LW(phiOld, c, nt)
    
    
    l2FTCS, errorFTCS = l2ErrorNorm(phiFTCS, phiAnalytic)
    l2FTBS, errorFTBS = l2ErrorNorm(phiFTBS, phiAnalytic)
    l2CTCS, errorCTCS = l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic)
    l2LW, errorLW = l2ErrorNorm(phiLW, phiAnalytic)
    
#    
    
    
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
    
#            
####Run the function main defined in this file                      ###
            
main()

def convergence_exp():
    n_exp_size = 30
    u = 0.4  #courant number constant
    xmin = 0
    xmax = 1
    l2FTBS_dx_err = np.zeros(n_exp_size)
    l2CTCS_dx_err = np.zeros(n_exp_size)
    l2LW_dx_err = np.zeros(n_exp_size)
    errorFTBS = np.zeros(n_exp_size)
    errorCTCS = np.zeros(n_exp_size)
    errorLW = np.zeros(n_exp_size)
    dx_it = np.zeros(n_exp_size)
    
    delta_x = np.zeros(n_exp_size)
    delta_x2 = np.zeros(n_exp_size)
    
    for i in range(0,n_exp_size):
        nx = i*10 + 10
        dx = (xmax - xmin)/nx
        nt = nx ##keeping overall time constant
        dx_it[i] = dx
        c = 0.4*(nx/nt)
     
        
        # spatial points for plotting and for defining initial conditions
        x = np.arange(xmin, xmax, dx)
        
        # Initial conditions
        phiOld = cosBell(x, 0.25, 0.75)
        # Exact solution is the initial condition shifted around the domain
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.25, 0.75)#, 0.5, 0.75)
        
        # Advect the profile using finite difference for all the time steps
        phiFTBS = FTBS(phiOld, c, nt)
        phiCTCS = CTCS(phiOld, c, nt)
        phiLW = LW(phiOld, c, nt)
        
        ##create dx anmd dx^2 lines
        delta_x[i] = dx
        delta_x2[i] = dx**2       
        
        l2FTBS_dx_err[i], errorFTBS = l2ErrorNorm(phiFTBS, phiAnalytic)
        l2CTCS_dx_err[i], errorCTCS = l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic)
        l2LW_dx_err[i], errorLW = l2ErrorNorm(phiLW, phiAnalytic)
    
    plt.figure(7, figsize=(10,7))
    plt.clf()
    plt.loglog(dx_it, l2FTBS_dx_err, label='FTBS', color = 'red')
    plt.loglog(dx_it, l2CTCS_dx_err, label='CTCS', color = 'green')
    plt.loglog(dx_it, l2LW_dx_err, label = 'LW', color = 'orange')
    plt.loglog(dx_it, delta_x, label='$\Delta x$', linestyle=':', color = 'black')
    plt.loglog(dx_it, delta_x2, label = '$\Delta x^{2}$', linestyle = ':', color = 'blue')
    plt.ylabel('$l_{2}$ error norm')
    plt.xlabel('Number of spacial steps')
    plt.legend()
    plt.title('Loglog plot of error norms as number of spacial points increase')
    
    
convergence_exp()
    
        
        
        
        
    
        

nx_list = (40, 80, 100, 120)
def c_exp():
    
    # Parameters
    xmin = 0
    xmax = 1
    
    
    for i in range(len(nx_list)):
        
        
        nx = nx_list[i]
        nt = 40
        
        u=0.4
        dx = (xmax - xmin)/nx
        c = u*(nx/nt)
        print(c)
        
        x = np.arange(xmin, xmax, dx)
    
        # Initial conditions
        phiOld = cosBell(x, 0.25, 0.75)
        # Exact solution is the initial condition shifted around the domain
        phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.25, 0.75)
        
        # Advect the profile using finite difference for all the time steps
        phiFTCS = FTCS(phiOld, c, nt)
        phiFTBS = FTBS(phiOld, c, nt)
        phiCTCS = CTCS(phiOld, c, nt)
        phiLW = LW(phiOld, c, nt)
        
        
        l2FTCS, errorFTCS = l2ErrorNorm(phiFTCS, phiAnalytic)
        l2FTBS, errorFTBS = l2ErrorNorm(phiFTBS, phiAnalytic)
        l2CTCS, errorCTCS = l2ErrorNorm(phiCTCS[nt-1,:], phiAnalytic)
        l2LW, errorLW = l2ErrorNorm(phiLW, phiAnalytic)
        
        font = {'size'   : 20}
        plt.rc('font', **font)
        plt.figure(i+2, figsize=(10,7))
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
        plt.legend(bbox_to_anchor=(0.4,1))
        plt.xlabel('$x$')
        plt.title('Linear Advection where c=%f'%c)
        
        print("FTBS l2 error norm = ", l2FTBS)
        print("FTBS linf error norm = ", lInfErrorNorm(phiFTBS, phiAnalytic))
        
        
        print("CTCS l2 error norm = ", l2CTCS)
        print("CSCS linf error norm = ", lInfErrorNorm(phiCTCS, phiAnalytic))
        
        
        print("LW l2 error norm = ", l2LW)
        print("LW linf error norm = ", lInfErrorNorm(phiLW, phiAnalytic))
       
    
        
c_exp()
        
def TV():
    
    # Parameters
    xmin = 0
    xmax = 1
    nx = 100
    nt = 100
    u=0.4  ##wind speed, keeping constant
    c = u*(nx/nt)
        
    # Derived parameters
    dx = (xmax - xmin)/nx
    
    # spatial points for plotting and for defining initial conditions
    x = np.arange(xmin, xmax, dx)
    
    # Initial conditions
    phiOld = cosBell(x, 0.25, 0.75)
    # Exact solution is the initial condition shifted around the domain
    phiAnalytic = cosBell((x - c*nt*dx)%(xmax - xmin), 0.25, 0.75)
    
    TV_FTBS = np.zeros(nx-2)
    TV_CTCS = np.zeros(nx-2)
    TV_LW = np.zeros(nx-2)
    phi_CTCS = np.zeros((nt,nx))
    phiCTCS = np.zeros((1, nx))
    for k in range(2,nt):  ##looping for each time step
        
        TVinter_FTBS = np.zeros(nt)
        TVinter_CTCS = np.zeros(nt)
        TVinter_LW = np.zeros(nt)
        
        phi_FTBS = FTBS(phiOld, c, k)
        phi_CTCS = CTCS(phiOld, c, k)
        phi_LW = LW(phiOld, c, k)
        phiCTCS = phi_CTCS[k-2,:] 
        
        for i in range(nx):  ##loop for each spacial difference
            
            
            
            TVinter_FTBS[i] = abs( phi_FTBS[(i+1)%nx] - phi_FTBS[i] )
            TVinter_CTCS[i] = abs( phiCTCS[(i+1)%nx] - phiCTCS[i] )
            TVinter_LW[i] = abs( phi_LW[(i+1)%nx] - phi_LW[i] )
            
        TV_FTBS[k-2] = sum(TVinter_FTBS)
        TV_CTCS[k-2] = sum(TVinter_CTCS)
        TV_LW[k-2] = sum(TVinter_LW)
        
    
    plt.figure()
    plt.plot(TV_FTBS, label='FTBS', color='blue')
    plt.plot(TV_CTCS, label='CTCS', color='green')
    plt.plot(TV_LW, label='LW', color='orange')
    plt.legend()
TV()
            
        
            
            
        
        
    
    
        
        
