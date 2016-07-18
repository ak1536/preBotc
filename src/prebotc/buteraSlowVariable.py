 
'''
Created on Jul 11, 2016

@author: andrewkennedy
'''

import numpy as np
from scipy.integrate import odeint 
import matplotlib.pyplot as plt
from scipy.special import expit
from integration import integrate
from utils import findClosest
from buteraSingleCell import Rhs as ButeraRHS


class Rhs(ButeraRHS):
    def __init__(self, N):
        super(Rhs, self).__init__(N)
        s = self
        s.numVarsPerCell = 4
        s.alpha = 3.0       # 1/s
        s.beta = 0.1        # 1/mV
        s.thetaApce = -30.0   # mV
        
   
    def rhs(self, V, h, n, Apce):
        dVdt, dhdt, dndt = super(Rhs, self).rhs(V, h, n)
        s = self
        gV = expit(s.beta *(V - s.thetaApce))
        dApcedt = s.alpha * (gV - Apce)
        return dVdt, dhdt, dndt, dApcedt
    
    def integrate(self, initial, tmin, tmax, **kwargs):
        return integrate(initial, self, tmin, tmax, giveTime=True, **kwargs)
