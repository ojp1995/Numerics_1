#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 15:57:21 2018

@author: ojp18
"""

#Attempting to write script for FTBS advection scheme
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
u = 0.1  #assuming constant speed
Nx=1000
dt = 0.001  #time step
dx = 1/(Nx-1)  #spacial step



def FTBS (x,dt):
    x_next=np.zeros(Nx)
    x_next[0]=x[0]
    for i in range (1, Nx):
        x_next[i] = x[i] -u*(dt/dx)*(x[i]-x[i-1])
    return x_next

#x_FTBS = FTBS(x,dt)

#fig=plt.figure()
#ax=plt.axes( xlim=(0,1), ylim=(0,1))
#line,=ax.plot([],[],lw=2)
    
x = np.zeros(Nx)  #creating vector of x with intial value 0 for all entries 
for i in range(0,Nx):
    if (i*dx>0.25 and i*dx<0.75):
        x[i] = 0.5*(1-np.cos(4*np.pi*(i*dx-0.25)))
        


for i in range(0,Nx):
    x_FTBS=FTBS(x,dt)
    plt.plot(x_FTBS)
  

#plt.plot(x_FTBS)
#plt.xlabel('Position x')
  

###animation  
#fig=plt.figure()
#ax=plt.axes( xlim=(0,1), ylim=(0,1))
#line,=ax.plot([],[],lw=2)
#
#def animate (frame):
#    global x
#    for i in range(0,500):
#        x=compute_step(x,dt)
#        line.set_data(np.linspace(0,1,Nx),x)
#        return line,
#    
#    amim = animation.FuncAnimation(fig,animate,frame=100,interval=int(dt*50000),blit=True)
#    plt.show()
