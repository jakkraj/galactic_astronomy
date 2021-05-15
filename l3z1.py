#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 00:39:04 2021

@author: kuba
"""

import numpy as np
import math
import statistics as st
from astropy import units as u
from astropy.coordinates import SkyCoord


data = np.loadtxt('l3z1data.txt', unpack=False, skiprows=0)



ra = data[:,1] + data[:,2]/60. + data[:,3]/3600.    #ra w stopniach
dec = data[:,4] + data[:,5]/60. + data[:,6]/3600.   #dec w stopniach
pm_ra = data[:,9]/1000./3600.    #ruch własny [stopnie łuku/rok]
pm_dec = data[:,11]/1000./3600.  #ruch własny [stopnie łuku/rok]
v_r = data[:,13]
par = data[:,7]/1000.

for i in range(len(ra)):
    ra[i] = math.radians(ra[i])
    dec[i] = math.radians(dec[i])
    pm_ra[i] = math.radians(pm_ra[i])
    pm_dec[i] = math.radians(pm_dec[i])



a = []
b = []
c = []

# 1. Wpółrzędne równikowe i galaktyczne punktu zbieżności.

# =============================================================================
# for x in range(len(data)):
#     a.append(pm_ra[x] * math.sin(dec[x]) * math.cos(dec[x]) * math.cos(ra[x]) - pm_dec[x] * math.sin(ra[x]))
#     b.append(pm_ra[x] * math.sin(dec[x]) * math.cos(dec[x]) * math.sin(ra[x]) - pm_dec[x] * math.cos(ra[x]))
#     c.append(pm_ra[x] * math.cos(dec[x])**2)
# 
# 
# x = []
# y = []
# xy = []
# x2 = []
# for i in range(len(a)):
#     x.append(-b[i]/a[i])
#     y.append(c[i]/a[i])
#     xy.append(x[i]*y[i])
#     x2.append(x[i]**2)
# =============================================================================
    


#Podejście 1
# =============================================================================
# H = (sum(xy)-len(a)*st.mean(x)*st.mean(y))/(sum(x2)-len(a)*st.mean(x)**2)
# G = st.mean(y) - H * st.mean(x)
# =============================================================================
    


# =============================================================================
# rek = coords.ra.rad
# dek = coords.dec.rad
# 
# A = np.array([np.cos(dek)*np.cos(rek), np.cos(dek)*np.sin(rek), np.sin(dek)]).T
# B = np.array(-coords.radial_velocity).T
# 
# v_sun = np.linalg.inv(A.T @ A) @ A.T @ B
# =============================================================================

A = np.array([pm_ra * np.sin(dec) * np.cos(dec) * np.cos(ra) - pm_dec * np.sin(ra), \
              pm_ra * np.sin(dec) * np.cos(dec) * np.sin(ra) - pm_dec * np.cos(ra)]).T

B = np.array(pm_ra * np.cos(dec)**2).T
    
result = np.linalg.inv(A.T @ A) @ A.T @ B 
    
H = result[0]
G = result[1]    
    
# =============================================================================
# #Podejście 2
# S = []
# Sx = []
# Sy = []
# Sxx = []
# Sxy = []
# 
# for i in range(len(data)):    
#     S.append(1./data[i,])
# =============================================================================
    
    
# Wpółrzędne równikowe punktu zbieżności:

ra_V = 2*math.pi + math.atan(H/G)
dec_V = math.atan(1./math.sqrt(G**2+H**2))
ra_V_deg = math.degrees(ra_V)
dec_V_deg = math.degrees(dec_V)


gal = SkyCoord(ra=ra_V_deg, dec=dec_V_deg, unit=(u.degree, u.degree))

print('Współrzędne równikowe punktu zbieżności: \n', 'RA: %.2f' % ra_V_deg, 'Dec: %.2f' % dec_V_deg)
print('Współrzędne galaktyczne punktu zbieżności: \n', gal.galactic)


# 2. Wartość prędkości przestrzennej gromady.

lam = []
s_v = []
cos_lam = []

for i in range(len(ra)):
    lam.append(math.sin(dec_V) * math.sin(dec[i]) + math.cos(dec_V) * math.cos(dec[i]) * math.cos(ra_V * ra[i]))
    cos_lam.append(math.cos(lam[i])**2)
    s_v.append(v_r[i]*math.cos(lam[i]))
    #v = v_r[i]
    v = (sum(s_v))/(sum(cos_lam))
    
    
vv = np.dot(lam,v_r) / np.dot(lam,lam)
    
    
print('Wartość prędkości przestrzennej gromady: ', v)

# 3. Odległości do poszczególnych gwiazd oraz ich średnia.

pm = []
r_st = []
r_pc = []
k = math.pi/math.radians(180.)
lam_mean = st.median(lam)

for i in range(len(pm_dec)):
    pm.append(math.sqrt((pm_dec[i]/31557600)**2 + (pm_ra[i]/31557600)**2)) #zamiana jednostek ruchów własnych z radianów na rok na radiany na sekundę
    r_st.append((v_r[i] * math.tan(lam_mean)) / (k * pm[i]))
    
for i in range(len(r_st)):
    r_pc.append(r_st[i]/(149597870*206264))

r_cl = st.mean(r_pc)
r_sdev = st.stdev(r_pc)
r_pc = [x for x in r_pc if x<r_cl + r_sdev]
r_cl = st.mean(r_pc)

# =============================================================================
# mu = np.sqrt(np.power(pm_dec/31557600,2) + np.power(pm_ra/31557600,2))
# r_i = np.divide(vv*np.sqrt(1.-np.power(lam,2)), mu*np.pi/180)
# r_i = r_i/(149597870*206264)
# r_mean = np.mean(r_i)
# r_std = np.std(r_i,ddof = 1)
# =============================================================================



#print('Odległość do poszczególnych gwiazd: ', r_st)
print('Średnia odległość: %.2f' % r_cl, 'pc.')


# 4. Średnia odległość do gromady z indywidualnych paralaks Hipparcos'a

paralax = []

for i in range(len(data)):
    paralax.append(1/par[i])

# =============================================================================
# r_par = st.median(paralax)
# r_dev = st.stdev(paralax)
# 
# paralax = [x for x in paralax if x<r_par + r_dev]
# =============================================================================
r_par = st.median(paralax)
r_dev = st.stdev(paralax)

print('Odległość do gromady (paralaksy): ', r_par)
