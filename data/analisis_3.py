from __future__ import division
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('md_3_prueba_512_seis_rhos_mejor.csv', delimiter=',', skiprows=1)
cant_rho = 6
rho = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

r = datos[:,0]

plt.figure()
for i in range(cant_rho):
    plt.plot(r, datos[:,i+1]/2, label = 'rho = %f' % rho[i])
plt.legend()
plt.show()
