from __future__ import division
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('md_3_prueba_512_seis_rhos_mejor.csv', delimiter=',', skiprows=1)
cant_rho = 6
rho = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

r = datos[:,0]
colors = {0: 'r', 2: 'b', 5: 'g'}

plt.figure()
for i in [0,2,5]: #range(cant_rho):
    plt.plot(r, datos[:,i+1]/2, label = '$ \\rho = %.1f $' % rho[i], linewidth = 3, \
             color = colors[i])
plt.grid(True)
plt.xlabel('$r$', fontsize = 15)
plt.ylabel('$g(r)$', fontsize = 15)
plt.title('Función de distribución radial', fontsize = 20)
ax = plt.gca()
ax.set_xticklabels(ax.get_xticks(), fontsize = 15)
ax.set_yticklabels(ax.get_yticks(), fontsize = 15)
plt.legend()
plt.show()
