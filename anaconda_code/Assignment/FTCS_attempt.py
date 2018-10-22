#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:22:48 2018

@author: ojp18
"""
#attempt to make it with FTCS
from numpy import zeros
from matplotlib import pyplot

u = 0.5  #constant speed
Nx = 100  #number of x points
Nt = 100  #number of time steps
dx = 0.001  #spacial step
dt = 0.001  #time step  
c = u*(dt/dx)  #courant number
phi = zeros((Nx,Nt))  #matrix of x position in row vectors with each new row being a new time step
 
phi[0,0] = 0.1  #initial conditions
phi[2,0] = 0.1
for n in range (0, Nt-1): #time loop
    for j in range (1, Nx-1): #space loop
        phi[n+1,j]=phi[n,j]-(c/2)*(phi[n,j+1]-phi[n,j-1])
        
#print (phi)
#pyplot.plot(phi[28,:])  #working plot but at specific times
for i in range(0, Nx-1):
    pyplot.figure(1)
    pyplot.plot(phi[:,i])