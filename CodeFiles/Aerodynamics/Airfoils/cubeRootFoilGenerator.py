# -*- coding: utf-8 -*-
"""
Generates an uncambered airfoil, with a cube root thickness distribution.  
This is used to test the structural generation algorithm.

@author: jdv5076
"""
import numpy as np

tauMax = 0.10
x_tauMax = 0.3
nAF     = 41

nS = nAF//2 + 1
xUS = np.zeros((nS,1))
xLS = np.zeros((nS,1))
dTheta = np.pi/(nS - 1)
for i in range(nS):
    xUS[i,0] = 0.5 - 0.5*np.cos(np.pi - i*dTheta)
    if xUS[i,0] <= 0.5:
        xUS[i,0] = xUS[i,0]*(x_tauMax/0.5)
    else:
        xUS[i,0] = ((1-x_tauMax)/(1-0.5))*(xUS[i,0] - 0.5) + x_tauMax
for i in range(nS):
    xLS[i,0] = 0.5 + 0.5*np.cos(np.pi - i*dTheta)
    if xLS[i,0] <= 0.5:
        xLS[i,0] = xLS[i,0]*(x_tauMax/0.5)
    else:
        xLS[i,0] = ((1-x_tauMax)/(1-0.5))*(xLS[i,0] - 0.5) + x_tauMax

x = np.concatenate([xUS,xLS[1:nS]])

y = np.zeros((nAF,1))
for i in range(nAF):
    if i<=nS-1:
        surf = 1
    else:
        surf = -1
    if x[i,0] <= x_tauMax:
        y[i,0] = 0.5*surf*tauMax*(x[i,0]/x_tauMax)**(1/3)
    else:
        y[i,0] = 0.5*surf*tauMax*((1-x[i,0])/(1-x_tauMax))**(1/3)

import matplotlib.pyplot as plt
plt.plot(x,  y,'k-')
# plt.gca().set_aspect("equal")
coords = np.concatenate([x,y],axis = 1)
print(coords)