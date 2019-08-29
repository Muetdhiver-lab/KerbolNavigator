# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 22:12:16 2019

@author: vince
"""

from pykep import planet, DEG2RAD, epoch, AU, propagate_lagrangian, par2ic, lambert_problem
from math import sqrt, pi, acos, cos, atan
from numpy.linalg import norm
from _Kerbol_System import Moho, Eve, Kerbin, Duna, Jool
from decimal import *
from matplotlib import pyplot as plt
import numpy as np
from scipy import array
from numpy import deg2rad,dot,cross,linalg,rad2deg

KAU = 13599840256 #m


    
class CavemanNavigator:
    
    ''' 
        goal
        From caveman data predict r2,v2, at t2 from available data.
        we need to use the launch to get an info on arg Pe
    '''
    
    def __init__(self):
        # Ugly, but bag it all here
        self.mu_c = 1.1723328e9*1e9
        self.r_c = 261600000
        self.Pe = 13339382571+self.r_c
        self.Ap = 71600522497+self.r_c
        self.SMA = (self.Pe+self.Ap)/2
        self.e = 1 - 2*self.Pe/(self.Pe+self.Ap)
        print(self.e)
        self.tL = self.setTime(1,229,4,12,23)
        self.t0 = self.setTime(1,231,0,34,15)
        self.t1 = self.setTime(1,323,3,27,15)
        self.t0_M = self.t0-self.tL
        self.t1_M = self.t1-self.tL
        self.rnL = 13338420256+self.r_c
        self.rn0 = 13341111200+self.r_c
        self.rn1 = 19701742109+self.r_c
        self.rL = [0,0,0]
        self.r0 = [0,0,0]
        self.r1 = [0,0,0]
        self.vn0 = 12044.7
        self.vn1 = 9487.6
        self.vL = [0,0,0]
        self.v0 = [0,0,0]
        self.v1 = [0,0,0]
        self.vL_data = [0,0,0]
        self.rL_data = [0,0,0]        
        self.N = sqrt(self.mu_c/(self.SMA**3))
        self.TAN_L = self.calcTAN(self.rnL,-1)
        self.TAN_t0 = self.calcTAN(self.rn0,1)
        self.TAN_t1 = self.calcTAN(self.rn1,1)
        self.EAN_L = self.calcEAN(self.TAN_L)
        self.EAN_t0 = self.calcEAN(self.TAN_t0)
        self.EAN_t1 = self.calcEAN(self.TAN_t1)        
        self.ArgPe = self.calcArgPe()

        
        
        
    def KerbinTAN(self,time):
        h2s = 3600
        y2s = 2556.5*h2s
        KTAN = (time/y2s)*2*pi +3.14 #will need modulo but hey.
        return KTAN
        

    
    def setTime(self, year, day, hour, minute, second):
        h2s = 3600
        d2s = 6*h2s
        y2s = 2556.5*h2s
        time = (year-1)*y2s + (day-1)*d2s + hour*h2s + minute*60 + second 
        print(time)
        return time
        
    def calcTAN(self,radius,sign):
        itRes = 1/self.e*((self.SMA/radius*(1-self.e**2))-1)
        print(itRes)
        if itRes > 1:
            itRes = 1
        TrueAnomaly = acos(itRes)
        if sign < 0: #sign 
            TrueAnomaly = -1*TrueAnomaly
        return TrueAnomaly

    def calcArgPe(self):
        #This stuff is likely wrong in some ways
        ArgPe = self.KerbinTAN(self.tL) - (self.TAN_L)
        return ArgPe
    
    def calcEAN(self,TAN):
        EAN = acos((self.e+cos(TAN))/(1+self.e*cos(TAN)))
        if TAN < 0:
            EAN =  -1*EAN
        return EAN

    
    def printVect(self,Vector,Name):
        print(Name,"         "," Norm : ", '{:12}'.format("{0:.0f}".format(norm(Vector))),"[","{0:.0f}".format(Vector[0]),",","{0:.0f}".format(Vector[1]),",","{0:.0f}".format(Vector[2]),"]")

    def deltasAtPoint(self,nR_pred,nV_pred,nR_data,nV_data,Name_r,Name_v):
        delta_nR = nR_data-nR_pred
        delta_nV = nV_data-nV_pred
        print("Delta for radius  ",Name_r,"  : ",'{:12}'.format("{0:.0f}".format(delta_nR)),"[m]  ","rel. err. = ","{0:.5f}".format(100*delta_nR/nR_data),"%")
        print("Delta for speed   ",Name_v,"  : ",'{:12}'.format("{0:.0f}".format(delta_nV)),"[m/s]","rel. err. = ","{0:.5f}".format(100*delta_nV/nV_data),"%")
        
    def calcOrbit(self):
        errEAN_L = 0  *  pi/180
        errArgPe = 0  *  pi/180  
        v0_pred = [0,0,0]
        r0_pred = [0,0,0]
        v1_pred = [0,0,0]
        r1_pred = [0,0,0]
        # we try to first convert keplerian elements to carthesian, then propagate
        print("------------------------------------------------------")
        self.rL_data, self.vL_data = par2ic([self.SMA,self.e,0,0,self.ArgPe+errArgPe,self.EAN_L+errEAN_L],self.mu_c)
        self.printVect(self.rL_data,"rL pred")
        self.printVect(self.vL_data,"vL pred")
        r0_pred, v0_pred = propagate_lagrangian(self.rL_data,self.vL_data,self.t0_M,self.mu_c)
        self.printVect(r0_pred,"r0 pred")
        self.printVect(v0_pred,"v0 pred")
        self.deltasAtPoint(norm(r0_pred),norm(v0_pred),self.rn0,self.vn0,"r0","v0")
        r1_pred, v1_pred = propagate_lagrangian(self.rL_data,self.vL_data,self.t1_M,self.mu_c)
        self.printVect(r1_pred,"r1 pred")
        self.printVect(v1_pred,"v1 pred")
        self.deltasAtPoint(norm(r1_pred),norm(v1_pred),self.rn1,self.vn1,"r0","v0")
        '''
        print("------------------------------------------------------")
        self.r0, self.v0 = par2ic([self.SMA,self.e,0,0,self.ArgPe+errArgPe,self.EAN_t0+errEAN_L],self.mu_c)
        print("r0       : ",self.r0," norm : ", norm(self.r0))
        print("v0       : ",self.v0," norm : ", norm(self.v0))      
        self.r1, self.v1 = propagate_lagrangian(self.r0,self.v0,self.t1_M-self.t0_M-600000,self.mu_c)
        print("r1       : ",self.r1," norm : ", norm(self.r1))
        print("v1       : ",self.v1," norm : ", norm(self.v1))
        '''
        
    def print(self):
        ArgPe = self.KerbinTAN(self.tL) - (self.TAN_L)
        print("Eccentricity     : ",self.e)
        print("------------------------------------------")
        print("TAN at Launch    : ",self.TAN_L*180/pi)
        print("EAN at Launch    : ",self.EAN_L*180/pi)        
        print("KTAN at Launch   : ",self.KerbinTAN(self.tL)*180/pi)
        print("Arg Pe           : ",ArgPe*180/pi)
        print("TAN at t0        : ",self.TAN_t0*180/pi)
        print("EAN at t0        : ",self.EAN_t0*180/pi)  
        print("TAN at t1        : ",self.TAN_t1*180/pi)
        print("EAN at t1        : ",self.EAN_t1*180/pi)  



    def correctionBurn(self):
        plt.rcParams['savefig.dpi']=100
        ManT = 323       # manoevre time
        ManT_W = 5     # manoevre window
        dy2s = 6*3600
        start_epochs = np.arange(ManT*0.25,(ManT+ManT_W)*0.25,0.25)
        ETA = 1000
        ETA_W = 2000
        duration = np.arange(ETA*0.25,(ETA+ETA_W)*0.25,0.25)


        
        #these are Kdays, to *0.25 to Edays (for eph function).
        #Kerbin = planet.keplerian(epoch(0), (13599840256 ,0,0,0,0,3.14), 1.1723328e18, 3.5316000e12,600000, 670000 , 'Kerbin')
        #Duna   = planet.keplerian(epoch(0), (20726155264 ,0.051,deg2rad(0.06) ,0,deg2rad(135.5),3.14), 1.1723328e18, 3.0136321e11,320000, 370000 , 'Duna')
        r2,v2 = Jool.eph(epoch(322.5*0.25)) #check jool position
        print(norm(r2))
        print(norm(v2))
        r1,v1 = propagate_lagrangian(self.rL_data,self.vL_data,322.5*dy2s-self.tL,self.mu_c)        
        print(norm(r1))
        print(norm(v1))       
        
        data = list()
        v_data = list()
        for start in start_epochs:
            row = list()
            v_row = list()
            for T in duration:
                #Need to replace the kerbin start point by the ship at time t using
                r1,v1 = propagate_lagrangian(self.rL_data,self.vL_data,(start-1)*dy2s,self.mu_c)

                #r1,v1 = Kerbin.eph(epoch(start))
                r2,v2 = Jool.eph(epoch(start+T))
                l = lambert_problem(r1,r2,T*60*60*24,Kerbin.mu_central_body) #K day = 6h
                DV1 = np.linalg.norm( array(v1)-array(l.get_v1()[0]) )
                v_DV1 = array(v1)-array(l.get_v1()[0]) 
                #DV2 = np.linalg.norm( array(v2)-array(l.get_v2()[0]) )
                #DV1 = max([0,DV1-4000])
                #DV = DV1+DV2
                DV = DV1
                #DV = sqrt(dot(DV1, DV1) + 2 * Kerbin.mu_self / Kerbin.safe_radius) - sqrt(Kerbin.mu_self / Kerbin.safe_radius )
                v_row.append(v_DV1)
                row.append(DV)
            data.append(row)
            v_data.append(v_row)
    

    
        minrows = [min(l) for l in data]
        i_idx = np.argmin(minrows)
        j_idx = np.argmin(data[i_idx])
        best = data[i_idx][j_idx]
        v_best = v_data[i_idx][j_idx]
        
        progrd_uv   = array(v1) / linalg.norm(v1)

        plane_uv   = cross(v1, r1)
        plane_uv   = plane_uv / linalg.norm(plane_uv)
        radial_uv   = cross(plane_uv, progrd_uv)
        EJBK      = sqrt(dot(v_best, v_best) + 2 * Kerbin.mu_central_body / norm(r1)) - sqrt(Kerbin.mu_central_body / norm(r1) )

        progrd_v = dot(progrd_uv, v_best)
        radial_v = dot(radial_uv, v_best)
        
        #print(rad2deg(atan(radial_v/progrd_v)))


        print("TransX escape plan - Kerbin escape")
        print("--------------------------------------")
        print("Best DV: " + str(best))
        print("Launch epoch (K-days): " +  str(start_epochs[i_idx]*4))
        print("Duration (K-days): " +  str(duration[j_idx]*4))
        print("Prograde:            %10.3f m/s" % np.round(dot(progrd_uv, v_best), 3))
        print("Radial:              %10.3f m/s" % np.round(dot(radial_uv, v_best), 3))
        print("Planar:              %10.3f m/s" % np.round(dot(plane_uv, v_best), 3))
        print("Hyp. excess velocity:%10.3f m/s" % np.round(sqrt(dot(v_best, v_best)), 3))
        #print("Earth escape burn:   %10.3f m/s" % np.round(EJBK, 3))


        duration_pl, start_epochs_pl = np.meshgrid(duration, start_epochs)
        plt.contour(start_epochs_pl*4,duration_pl*4,array(data),levels = list(np.linspace(best,3000,12)))
        #plt.imshow(array(data).T, cmap=plt.cm.rainbow, origin = "lower", vmin = best, vmax = 5000, extent=[0.0, 850, 10, 470.0], interpolation='bilinear')

        #plt.colorbar(im);
        plt.colorbar()
        plt.show()

TestNav = CavemanNavigator()
TestNav.print()       
TestNav.calcOrbit()
TestNav.correctionBurn()

#plot_innerKerbol(epoch(100))