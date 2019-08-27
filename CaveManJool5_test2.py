# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 22:12:16 2019

@author: vince
"""

from pykep import planet, DEG2RAD, epoch, AU
from math import sqrt, pi, acos
from astropy import units as u

from _Kerbol_System import Moho, Eve, Kerbin, Duna, Jool

KAU = 13599840256 #m



    

    
class CavemanNavigator:
    
    ''' 
        goal
        From caveman data predict r2,v2, at t2 from available data.
        we need to use the launch to get an info on arg Pe
    '''
    
    def __init__(self):
        self.mu_c = 1.1723328e9
        self.Pe = 12267590602
        self.Ap = 29130431165
        self.SMA = (self.Pe+self.Ap)/2
        self.e = 1 - 2*self.Pe/(self.Pe+self.Ap)
        self.t0 = self.setTime(1,64,4,33,58)
        self.t1 = self.setTime(1,123,4,23,30)
        self.tL = self.setTime(1,39,3,13,30)
        self.r0 = 12432026500
        self.r1 = 13484050000
        self.v0 = 11348.2
        self.v1 = 10707.2
        self.N = sqrt(self.mu_c/(self.SMA**3))

        self.rL = 13338420256
        self.TAN_L = self.calcTAN(self.rL)
        self.TAN_t0 = self.calcTAN(self.r0)
        self.TAN_t1 = self.calcTAN(self.r1)
        
    def KerbinTAN(self,time):
        h2s = 3600
        y2s = 2556.5*h2s
        KTAN = (time/y2s)*360 #will need modulo but hey.
        return KTAN
        
        # return 
    
    def setTime(self, year, day, hour, minute, second):
        h2s = 3600
        d2s = 6*h2s
        y2s = 2556.5*h2s
        time = (year-1)*y2s + (day-1)*d2s + hour*h2s + minute*60 + second 
        return time
        
    def calcTAN(self,radius):
        TrueAnomaly = acos(1/self.e*((self.SMA/radius*(1-self.e**2))-1))*180/pi
        return TrueAnomaly
       
    def print(self):
        ArgPe = self.KerbinTAN(self.tL) - (self.TAN_L)
        print("TAN at Launch    : ",self.TAN_L)
        print("KTAN at Launch   : ",self.KerbinTAN(self.tL))
        print("Arg Pe           : ",ArgPe)
        print("TAN at t0        : ",90+self.TAN_t0)
        print("TAN at t1        : ",90+self.TAN_t1)
 
TestNav = CavemanNavigator()
TestNav.print()       
#plot_innerKerbol(epoch(100))