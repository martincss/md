from __future__ import division
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('md_1a_paso_1e-4.csv', delimiter=',', skiprows=1)

K = datos[:,0]
U = datos[:,1]
E = datos[:,2]

plt.figure()
plt.plot(K, 'r', U, 'b', E, 'g')
plt.show()

paso_term = 3500
K = datos[3500:,0]
U = datos[3500:,1]
E = datos[3500:,2]

U -= np.mean(U)
K -= np.mean(K)

def autocorr(x):
    result = np.correlate(x, x, mode='full')/np.sum(x**2)
    return result[result.size//2:]

rho_U = autocorr(U)
rho_K = autocorr(K)

plt.figure()
plt.plot(rho_U, 'b', rho_K, 'r')
plt.show()
