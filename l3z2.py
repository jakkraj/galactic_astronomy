#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 19:11:35 2021

@author: kuba
"""

import numpy as np
import statistics as st
import math
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, CartesianRepresentation, CartesianDifferential
from astropy.coordinates import ICRS

data = np.loadtxt('l3z2data.txt', unpack=False, skiprows=0)


# 1. Środek geometryczny układu gromad oraz jego współrzędne w prostokątnym układzie galaktycznym.

ra = data[:,0] + data[:,1]/60 + data[:,2]/3600    #ra w radianach  
dec = data[:,3] + data[:,4]/60 + data[:,5]/3600   #dec w radianach
r = data[:,6]#*1000 #pc
v_r = data[:,8] #km/s


ra1 = []
dec1 = []
r1 = []
v_r1 = []

# Zamiana ra i dec na radiany
for i in range(len(data)):
    ra1.append(data[i,0] + data[i,1]/60. + data[i,2]/3600.)
    dec1.append(data[i,3] + data[i,4]/60. + data[i,5]/3600.)
    r1.append(data[i,6])
    v_r1.append(data[i,8])
    
    ra[i] = math.radians(ra[i])
    dec[i] = math.radians(dec[i])


coords = SkyCoord(ra=ra1, dec=dec1, distance=r1, radial_velocity=v_r1*u.km/u.s,\
                  unit=(u.hourangle, u.deg, u.kpc), frame='icrs')

coords_gal = coords.galactic    
    


# Współrzędne bieguna galaktyki i bieguna nieba (l_P)
ra_G = math.radians(12. + 51.4/60)
dec_G = math.radians(27. + 8./60)
l_P = math.radians(123.)

# Współrzędne centrum galaktyki

ragc = 17. + 45./60. + 40.04/3600.
decgc = -29 + (-28.1)/3600.
rgc = 8.178

coords_gc = SkyCoord(ra=ragc*u.hourangle, dec=decgc*u.deg, distance=rgc*u.kpc, frame='icrs')
xyz_gc = CartesianRepresentation(coords_gc.cartesian.xyz)
gc = coords_gc.galactic
gc_gal = CartesianRepresentation(gc.cartesian.xyz)



clusters_center = CartesianRepresentation(coords.cartesian.xyz).mean()
clusters_center_gal = CartesianRepresentation(coords_gal.cartesian.xyz)#.mean()




print('Współrzędne prostokątne środka układu gromad:', clusters_center.get_xyz())
print('===============================================')
print('Odległość Słońca od środka układu gromad:', clusters_center.norm())


# 3. Ruch Słońca ku apeksowi względem centroidu gromad opartego na prędkościach radialnych; składowe prędkości Słońca



rek = coords.ra.rad
dek = coords.dec.rad
l = coords_gal.l.rad
b = coords_gal.b.rad


A = np.array([np.cos(dek)*np.cos(rek), np.cos(dek)*np.sin(rek), np.sin(dek)]).T
B = np.array(-coords.radial_velocity).T

A1 = np.array([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)]).T
B1 = np.array(-coords.radial_velocity).T

v_sun = np.linalg.inv(A.T @ A) @ A.T @ B
v_sun_gal = np.linalg.inv(A1.T @ A1) @ A1.T @ B1




# ==================================
# Współrzędne galaktyczne



sun_relative_to_cluster_center = ICRS(x = clusters_center.x, y = clusters_center.y, z = clusters_center.z,\
                                     v_x = v_sun[0]*u.km/u.s, v_y=v_sun[1]*u.km/u.s, v_z=v_sun[2]*u.km/u.s,\
                                     representation_type=CartesianRepresentation, differential_type=CartesianDifferential)


gal = sun_relative_to_cluster_center.spherical




print('===============================================')
print('Składowe prędkości Słońca (układ równikowy): ')
print(v_sun)


print('Składowe prędkości Słońca (układ galaktyczny): ')
print(v_sun_gal)






# ================================================
# Wykresy

    
plt.plot(coords_gal.cartesian.x, coords_gal.cartesian.y, 'ro', markersize=4)
plt.plot(gc_gal.x, gc_gal.y, marker='o', markersize=5, color="blue")
plt.plot(0, 0, marker='o', markersize=5, color="yellow")
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.show()
plt.savefig('xy.png', dpi=300)

plt.plot(coords_gal.cartesian.x, coords_gal.cartesian.z, 'ro', markersize=5)
plt.plot(gc_gal.x, gc_gal.z, marker='o', markersize=5, color="blue")
plt.plot(0, 0, marker='o', markersize=5, color="yellow")
plt.xlabel('x [kpc]')
plt.ylabel('z [kpc]')
plt.show()
plt.savefig('xz.png', dpi=300)



#plt.hist(coords_gal.l, bins=50)
plt.hist(np.array(coords_gal.l), bins=50)
#plt.axis([0, 360, 0, 12])
plt.xlabel('Długość galaktyczna [stopnie]')
plt.show()
plt.savefig('hist.png', dpi=300)

