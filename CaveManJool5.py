# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 22:12:16 2019

@author: vince
"""

from pykep import planet, DEG2RAD, epoch, AU
from Math import sqrt, PI


from _Kerbol_System import Moho, Eve, Kerbin, Duna, Jool

KAU = 13599840256 #m


def plot_innerKerbol(epoch = epoch(0)):
    
    """
    Plots the Galilean Moons of Jupiter at epoch

    USAGE: plot_moons(epoch = epoch(34654, epoch.epoch_type.MJD)):
    * epoch: the epoch one wants the galilean moons to be plotted at
	
    """
    from pykep.orbit_plots import plot_planet
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    print(epoch)
    fig = plt.figure()
    ax1 = fig.gca(projection='3d')

	#plot_planet(ax,Moho,color = 'y', units = 1.0, t0 = epoch, legend=True)
	#plot_planet(ax,Eve,color = 'm', units = 1.0, t0 = epoch, legend=True)
	#lot_planet(ax=ax1,Kerbin,color = 'c', units = 1.0, t0, legend=True)
    #plot_planet(ax = ax1,Duna,color = 'r', units = 1.0, t0 = epoch, legend=True)
    plot_planet(Moho, t0=epoch, color = 'y', legend=True, units=KAU, ax=ax1)
    plot_planet(Eve, t0=epoch, color = 'm', legend=True, units=KAU, ax=ax1)
    plot_planet(Kerbin, t0=epoch, color = 'c', legend=True, units=KAU, ax=ax1)
    plot_planet(Duna, t0=epoch, color = 'r', legend=True, units=KAU, ax=ax1) 
    plot_planet(Jool, t0=epoch, color = 'g', legend=True, units=KAU, ax=ax1) 
    ax1.set_xlim3d(-3,3)
    ax1.set_ylim3d(-3,3)
    ax1.set_zlim3d(-3,3)
    plt.show()
    
    """
    Toying with basic physics and limited input data
    
    What we got :
        > Ap
        > Pe
        > SMA
        > mu_c
        > altitude |r|
        > speed |v|
        > times
        
        And that's pretty much it.
        
        e = (Rp*Vp^2)/(mu_c) - 1
        
        e = 1 - 2*Pe/(Pe+Ap)
        
        Kepler : M - Mo = n (t - to)
        
        n = sqrt(mu_c/a^3)
        
        cos E = (e+cos(nu))/(1+e*cos(nu)) nu = true anomaly, E excentric anomaly
        
        if e small, nu ~ M + 2e*sin(M) + 1.25e^2*sin(2M)
        
        We don't have M, but we have t and to and mu_c
        
    """
    

    
class CavemanNavigator:
    
    ''' 
        goal
        From caveman data predict r2,v2, at t2 from available data.
    '''
    
    def __init__(self):
        self.Pe = 10157
        self.Ap = 10496
        self.SMA = (self.Pe+self.Ap)/2
        self.e = 1 - 2*self.Pe/(self.Pe+self.Ap)
        self.t0 = 0
        self.t1 = 0
        self.r0 = 0
        self.r1 = 0
        self.v0 = 0
        self.v1 = 0
        self.N = sqrt(Kerbin.mu_c/(self.SMA^3))
        
    def setT0(self, year, day, hour, minute, second):
        h2s = 3600
        d2s = 6*h2s
        y2s = 2556.5*h2s
        self.t0 = year*y2s + day*d2s + hour*h2s + minute*60 + second 
    
    def setT1(self, year, day, hour, minute, second):
        h2s = 3600
        d2s = 6*h2s
        y2s = 2556.5*h2s
        self.t0 = year*y2s + day*d2s + hour*h2s + minute*60 + second 
        
    def EccAnom(self):
        K = PI/180
        max_i = 50
        i = 0
        delta = 0.0001
        
        
    
plot_innerKerbol(epoch(100))