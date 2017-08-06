from __future__ import division
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('md_2_prueba.csv', delimiter=',', skiprows=1)
rho_inicial = 0.4
rho_final = 0.8
temp_inicial = 2.0
temp_final = 0.4
cant_rhos = 5
cant_temps = 10

list_rhos = np.linspace(rho_inicial, rho_final, cant_rhos)
list_temps = np.linspace(temp_inicial, temp_final, cant_temps)
 
rho = datos[:,0]
temp = datos[:,1]
E = datos[:,2]
P = datos[:,3]

# filas para rho, columnas para t
energia = np.zeros((cant_rhos, cant_temps))
presion = np.zeros((cant_rhos, cant_temps))

for i in range(cant_rhos):
    for j in range(cant_temps):
        energia[i,j] = E[i*cant_temps + j]
        presion[i,j] = P[i*cant_temps + j]
        
plt.figure()
for i in range(cant_temps):
    plt.subplot(1,2,1)
    plt.plot(list_rhos, energia[:,i])
    plt.subplot(1,2,2)
    plt.plot(list_rhos, presion[:,i])

plt.figure()
for i in range(cant_rhos):
    plt.subplot(1,2,1)
    plt.plot(list_temps, energia[i,:])
    plt.subplot(1,2,2)
    plt.plot(list_temps, presion[i,:])
                