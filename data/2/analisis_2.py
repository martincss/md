from __future__ import division
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('md_2_125_r8_t20.csv', delimiter=',', skiprows=1)
rho_inicial = 0.1
rho_final = 0.8
temp_inicial = 2.0
temp_final = 0.4
cant_rhos = 8
cant_temps = 20

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

#%%

plt.figure()
plt.suptitle('$N = 125$', fontsize = 20)
for i in [2*j for j in range(10)]: #range(cant_temps):
    plt.subplot(1,2,1)
    plt.plot(list_rhos, energia[:,i], label = 'T = %.2f' % list_temps[i])
    plt.subplot(1,2,2)
    plt.plot(list_rhos, presion[:,i], label = 'T = %.2f' % list_temps[i])
plt.subplot(1,2,1)
plt.grid(True)
plt.xlabel('$\\rho$', fontsize = 15)
plt.ylabel('Energía', fontsize = 15)
plt.title('Energía en función de la densidad', fontsize = 20)
#ax = plt.gca()
#ax.set_xticklabels(ax.get_xticks(), fontsize = 15)
#ax.set_yticklabels(ax.get_yticks(), fontsize = 15)
plt.legend()

plt.subplot(1,2,2)
plt.grid(True)
plt.xlabel('$\\rho$', fontsize = 15)
plt.ylabel('Presión', fontsize = 15)
plt.title('Presión en función de la densidad', fontsize = 20)
#ax = plt.gca()
#ax.set_xticklabels(ax.get_xticks(), fontsize = 15)
#ax.set_yticklabels(ax.get_yticks(), fontsize = 15)
plt.legend()

#%%

plt.figure()
plt.suptitle('$N = 125$', fontsize = 20)
for i in range(cant_rhos):
    plt.subplot(1,2,1)
    plt.plot(list_temps, energia[i,:], label = '$\\rho = %.2f $' % list_rhos[i])
    plt.subplot(1,2,2)
    plt.plot(list_temps, presion[i,:], label = '$\\rho = %.2f $' % list_rhos[i])

plt.subplot(1,2,1)
plt.grid(True)
plt.xlabel('Temperatura', fontsize = 15)
plt.ylabel('Energía', fontsize = 15)
plt.title('Energía en función de la temperatura', fontsize = 20)
#ax = plt.gca()
#ax.set_xticklabels(ax.get_xticks(), fontsize = 15)
#ax.set_yticklabels(ax.get_yticks(), fontsize = 15)
plt.legend()

plt.subplot(1,2,2)
plt.grid(True)
plt.xlabel('Temperatura', fontsize = 15)
plt.ylabel('Presión', fontsize = 15)
plt.title('Presión en función de la temperatura', fontsize = 20)
#ax = plt.gca()
#ax.set_xticklabels(ax.get_xticks(), fontsize = 15)
#ax.set_yticklabels(ax.get_yticks(), fontsize = 15)
plt.legend()
