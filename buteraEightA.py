'''
Created on Jul 11, 2016

@author: andrewkennedy
'''

import numpy as np
from scipy.integrate import odeint 
import scipy.special as spec
import matplotlib.pyplot as plt
from scipy.special import expit


class Rhs(object):
    def __init__(self, N):
        s = self
        s.N = N    # Note that all units are in terms of milli (10**-3)
        s.A = np.zeros((N , N))
        s.gNa = 2.8    # nS
        s.VNa = 50    # 50mV
        s.gk = 11.2    # 11.2nS
        s.Vk = -85.0   # -85mV
        s.gNap = 2.8    # 2.8nS
        s.gl = 2.8    # 2.8nS
        s.Vl = -65.0   # -65mV 
        s.Iapp = 0.0    # Zero unless otherwise noted
        s.C = 21.0   # 21pF
   
    def rhs(self, V , h, n):
        s = self
        tauhV = 0.010 / (np.cosh((V + 48.0) / (12.0)))     # Conversion from ms to s
        taunV = 10/ (np.cosh((V + 29.0) / (-8.0))) # Conversion from ms to s
        hinfV = expit((V + 48.0)/ (6.0))
        minfNaPV = spec.expit((V + 40.0)/ (-6.0))
        minfNaV = spec.expit((V + 34.0)/ (-5.0))
        ninfNaV = spec.expit((V + 29.0)/ (-4.0)) 
        
        # m n h need expit functions
        # V h n are the function parameters
        # dVdt = [- Inaph - INa - IK - IL - Itonice==0 + Iapp==0 ] / C
        Iapp = 182.000000
        dVdt = (
                (-s.gNap * minfNaPV * h * (V - s.VNa))
                + ((-s.gNa * (minfNaV**3) * (1 - n)) * (V - s.VNa))
                + ((-s.gk * (ninfNaV**4)) * (V - s.Vk))
                + (-s.gl * (V - s.Vl))
                + Iapp) / s.C *1000 #Conversion to mV because C is defined in terms of Farad Volts
        dhdt = (hinfV - h) / tauhV 
        dndt = (ninfNaV - n) / taunV
        return dVdt, dhdt, dndt
    
    def __call__(self, X, t):
        N = self.N
        V = X[(0*N):(1*N)]
        h = X[(1*N):(2*N)]
        n = X[(2*N):(3*N)]
        dVdt, dhdt, dndt = self.rhs(V, h, n)
        return np.hstack((dVdt.reshape(N, ), dhdt.reshape(N, ), dndt.reshape(N, ))) 
    def integrate(self, tmin, tmax, Nvals=1000):
        N = self.N
        Tvals = np.linspace(tmin, tmax, Nvals)
        Xvals = odeint(self,np.zeros((3*N, )),Tvals) 
        return Xvals, Tvals
    
        
        

N = 64
V = np.zeros((N , 1))
h = np.zeros((N , 1))
n = np.zeros((N , 1))
r = Rhs(N)
A = np.random.uniform(low = 0, high = 1, size = (N, N))
for i in range(N):
    for j in range(i, N):
        if i == j :
            A[i, j] = 0
        else:
            A[i, j] = A[j, i]
r.A = A
nominal = r.gNa
r.gNa = np.random.uniform(low = nominal*0.9, high = nominal*1.1, size = (N,))
        
Xvals, Tvals = r.integrate(0, 1)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = "3d")
X = Xvals.reshape(Tvals.size, 3, N)
Vi = X[500, 0, :]
# Vi = Xvals[500, :N]
ax.scatter(r.gNa, r.A.sum(1), Vi)
ax.set_xlabel('gNa [nS]')
ax.set_ylabel('degree')
ax.set_zlabel('V [mV]')

fig, ax = plt.subplots()
ax.plot(Tvals, X[: , 0, :])
ax.set_xlabel('Time')
ax.set_ylabel('MilliVolts')
plt.show()







