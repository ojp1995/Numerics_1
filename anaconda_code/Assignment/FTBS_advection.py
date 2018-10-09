#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:57:21 2018

@author: ojp18
"""

#Attempting to write script for FTBS advection scheme
import numpy as np
import matplotlib.pyplot as plt
u = 0.1  #assuming constant speed
Nx=1000


x = np.zeros(Nx)  #creating vector of x with intial value 0 for all entries 
x_next=np.zeros(Nx)
dt = 0.001  #time step
dx = 0.001  #spacial step

x[0] = 1  #intial condition


    
def compute_step (x,dt):
    x_next=np.zeros(len(x))
    x_next[0]=x[0]
    for i in range (1, len(x)):
        x_next[i] = x[i] -u*(dt/dx)*(x[i]-x[i-1])
    return x_next


for i in range(0,5000):
    x=compute_step(x,dt)
    plt.draw()
    #time.sleep(0.0005)

plt.plot(np.linspace(0,Nx*dt,Nx),x)
plt.xlabel('Position x')

