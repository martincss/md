# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 22:37:38 2017

@author: martin
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

datos = np.loadtxt('n216_T100_t100.csv', delimiter = ',', skiprows = 1)
tiempo_tot = datos.shape[0]//3
n_part = datos.shape[1]

x_part = np.zeros((tiempo_tot, n_part))
y_part = np.zeros((tiempo_tot, n_part))
z_part = np.zeros((tiempo_tot, n_part))

for t in range(tiempo_tot):
    for i in range(n_part): 
        x_part[t,i] = datos[3*t,i]
        y_part[t,i] = datos[3*t+1,i]
        z_part[t,i] = datos[3*t+2,i]
        
#%%

fig = plt.figure()
ax3 = Axes3D(fig)
for t in range(tiempo_tot):
    ax3.set_xlim3d(0, 8.3)
    ax3.set_ylim3d(0, 8.3)
    ax3.set_zlim3d(0, 8.3)
    ax3.scatter(x_part[t,:], y_part[t,:], z_part[t,:])
    #plt.savefig('D:\\martin\\Documents\\GitHub\\md\\data\\frames\\frame_{:03d}.png'.format(t))
    plt.savefig('/home/martin/Documents/md/data/animation/frames/frame_{:03d}.png'.format(t))
    #plt.pause(0.001)
    ax3.clear()
plt.close()
    
    
    
    