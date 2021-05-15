#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 16:10:51 2021

@author: kuba
"""

import math
import matplotlib.pyplot as plt
#import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import statistics as st

b = 0.
l = list(range(1,361))
k = 1/(180.*60.*60.*1000.)
r_sun = 8.5 * 3*10**16 #km
v_sun = 220. #km/s
omega = v_sun/r_sun
u_lsr = -10.
#u_lsr = 0
v_lsr = 5.2
#v_lsr = 0
w_lsr = 7.2
#w_lsr = 0
#v_sun = math.sqrt(u_lsr**2 + v_lsr**2 + w_lsr**2)

#r = random.sample(range(0., 16.), 360)


r = []
vr_list = []
vl_list = []
vb_list = []
r_gc = []
l_gc = []
b_gc = []
l_rad = []
zero = []
pm = [] #ruch własny = ruch własny w dł. gal., bo ruch własny w szerokości = 0 (b=0)
for i in range(len(l)):
    l_rad.append(math.radians(l[i]))
    
    r.append(np.random.uniform(8, 9))
    r[0] = 10. 
    
    vr_list.append(u_lsr * math.cos(b) * math.cos(l_rad[i]) - v_lsr * math.cos(b) * math.sin(l_rad[i]) - \
                   w_lsr * math.sin(b))
    
    vl_list.append(-u_lsr * math.sin(l_rad[i]) - v_lsr * math.cos(l_rad[i]) + omega*r[i] * math.cos(b))
    zero.append(0.)
    
    vb_list.append(-u_lsr * math.sin(b) * math.cos(l_rad[i]) + v_lsr * math.sin(b) * math.sin(l_rad[i]) - w_lsr * math.cos(b))

    #pm.append((vl_list[i]*60*60*24*365.25)/(180*60*60*1000))
    pm.append(((vl_list[i]*60*60*24*365.25)/(r[i]*3*10**16))*(180/math.pi)*60*60*1000)
    
    # Względem GC
    
    r_gc.append(0.)
    l_gc.append(v_sun*math.cos(b))
    b_gc.append(0.)
    

A = 0.

#B = abs(0.-st.mean(vr_list))
B = omega

print('Stałe Oorta:')
print('A:', A)
print('B:', B)

#Axes3D.scatter(np.array(vr_list), np.array(vl_list), np.array(vb_list))
#ax.plot_wireframe(np.array(vr_list), np.array(vl_list), np.array(vb_list))
    
plt.plot(l, vr_list)
plt.plot(l, vl_list)
plt.plot(l, vb_list)
plt.plot(l, zero, c='k')
plt.legend(('v_r', 'v_l', 'v_b'), loc='upper right')
plt.xlabel('l [deg]')
plt.ylabel('v [km/s]')
plt.show()

# =============================================================================
# plt.plot(l, r_gc)
# plt.plot(l, l_gc)
# plt.plot(l, b_gc)
# plt.legend(('v_r', 'v_l', 'v_b'), loc='upper right')
# plt.xlabel('l [deg]')
# plt.ylabel('v_r [km/s*kpc]')
# plt.show()
# =============================================================================

plt.scatter(l, pm, c=r)
plt.xlabel('l [deg]')
plt.ylabel('pm [marcsec/yr]')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(np.array(vr_list), np.array(vl_list), np.array(vb_list))
ax.set_xlabel('v_r [km/s]')
ax.set_ylabel('v_l [km/s]')
ax.set_zlabel('v_b [km/s]')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(np.array(r_gc), np.array(l_gc), np.array(b_gc))