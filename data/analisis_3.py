from __future__ import division
# coding=utf-8

import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt('md_3_prueba_posta_3.csv', delimiter=',', skiprows=1)

r = datos[:,0]
g_1 = datos[:,1]
g_2 = datos[:,2]
g_3 = datos[:,3]

plt.figure()
plt.plot(r, g_1, 'r', r, g_2, 'b', r, g_3, 'g')
plt.show()
