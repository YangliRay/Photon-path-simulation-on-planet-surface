# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 16:10:50 2022

@author: 1229290416
"""

import matplotlib.pyplot as plt
import numpy as np
import random

dl = 0.01

#radius of the earth R
R = 42

#radius of the atmosphere R_max
R_max = 50

#temperaturn as a function of height
def T(r):
    return 1

#number density as a function of height
def n(r):
    return 0.0021 * np.exp(300/r)

#collision cross section area as a function of height
def sigma(r):
    return 0.2


class photon():
    def __init__(self, x, y, z, v):
        
        #position
        self.x = x
        self.y = y
        self.z = z
        
        #velocity direction
        self.v = v
        
        #radius distance
        self.r = np.sqrt(x**2 + y**2 + z**2)
        
        #Has enter the atmosphere or not
        self.inside = False
        
    def deflect_velocity(self, theta, phi):
        
        vx = self.v[0]
        vy = self.v[1]
        vz = self.v[2]
        
        n1 = np.cos(theta)
        n2 = (np.sin(theta) * np.cos(phi))/(vz**2 + vy**2)
        n3 = (np.sin(theta) * np.sin(phi))/(vz**4 + vy**4 + 2* vz**2 * vy**2 
                                            + vx**2 * vy**2 + vz**2 * vx**2)
        
        self.v = [
            vx * n1 - (vz**2 + vy**2) * n3, 
            vy * n1 - vz * n2 + vy * vx * n3, 
            vz * n1 + vy * n2 + vz * vz * n3]
        
    def update_position(self):
        self.x = self.x + self.v[0] * dl
        self.y = self.y + self.v[1] * dl
        self.z = self.z + self.v[2] * dl
        self.r = np.sqrt(self.x**2 + self.y**2 + self.z**2)

#shot one photon from position (x, y, z) with direction v
def shot_photon(x, y, z, v):
    o = photon(x, y, z, v)
    while o.r > R:
        current_height = o.r
        o.update_position()
        
        #if inside the atmosphere
        if o.inside == True:
            test_for_collision = random.random()
            
            #photon escape the atmosphere
            if o.r > R_max:
                return [False, []]
            
            #if collision with a particle is true
            if test_for_collision <= n(current_height) * sigma(current_height) * dl:
                
                #random choose a scattering angle theta by Rayleigh phase function
                scattering_theta = random.random()
                a = 4 - 8 * scattering_theta
                b = pow((a + np.sqrt(a ** 2 + 4))/2, 1/3)
                theta = np.arccos(b - 1/b)
                
                #random choose a scattering angle phi between 0 and 2 * pi
                phi = random.random() * 2 * np.pi
                
                #update velocity
                o.deflect_velocity(theta, phi)
            
            else:
                pass
                
        elif o.r < R_max:
            o.inside = True
    
    #if reach earth's surface, return the coordinate
    if o.r <= R:
        return [True, 
                [np.arctan(o.y/o.x), np.arctan(o.z/np.sqrt(o.x**2 + o.y**2))]]


#generate photon
points = []
for i in range(0, 300000, 1):
    prob = random.random()
    position = np.sqrt(prob * (R_max**2))
    o = shot_photon(position, 0, 51, [0, 0, -1])
    
    if o[0] == True:
        
        #record the phi coordinate of the photon who reach earth surface
        points.append(o[1][1])

#plot the histogram
bin_number = 40

frequency_hist = np.histogram(points, bins = bin_number)
frequency = list(frequency_hist[0]) + [0]
phi_values = list(frequency_hist[1])

#plot frequency distribution
plt.plot(phi_values, frequency, "-", label = 'frequency', color = 'r')

for i in range(0, bin_number, 1):
    
    #convert frequncy into intensity (ratio)
    frequency[i] = frequency[i] / np.cos(phi_values[i])

#plot intensity distribution
plt.plot(phi_values, frequency, "-", label = 'Intensity', color = 'k')
plt.legend(fontsize=10)

plt.savefig("photon_intensity.jpg", dpi = 500)







