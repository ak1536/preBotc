'''
Created on Jul 11, 2016

@author: andrewkennedy
'''
import numpy as np
from scipy.integrate import odeint 

class Rhs(object):
    def __init__(self, N):
        s = self
        s.N = N    # Note that all units are in terms of milli (10**-3)
        s.A = np.zeros((N , N))
        s.gNa = 0.000028    # 28nS
        s.VNa = 50.0    # 50mV
        s.gk = 0.0000112    # 11.2nS
        s.Vk = -0.000085    # -85mV
        s.gNap = 0.0000028    # 2.8nS
        s.gNaph = 0.0000028    # 2.8nS  Check this value
        s.gl = 0.0000028    # 2.8nS
        s.Vl = -60.0    # -65mV 
        s.Iapp = 0.0    # Zero unless otherwise noted
        s.C = 0.000000021    # 21pF
        s.eps = 0.1    # The maximal tau value is located at V == thetax ==thetah == thetam????
            
    def rhs(self, V , h, n, s):
        s = self
        sV = 1 / (1 + np.exp( - (V+0.034) / (-0.005)))    # May not be accurate
        tauV = 1 / (s.eps * np.cosh((V + 0.034) / (-0.010)))    # Thetam and Sigmam are used for thetah and sigmah (where h is set equal to x)
        hinfV = 1 / (1 + np.exp((V + 0.034) / (-0.005)))    # Thetam and Sigmam are used for thetah and sigmah (where h is set equal to x)
        mV = 1 / (1 + np.exp( - (V + 0.037) / -0.006))    # May not be accurate
        

        dVdt = ((-s.gNaph * mV * h * (V - s.VNa)) + (-s.gNa * mV**3 * (1 - s.N)) * (V - s.VNa) + (-s.gk * s.N**4) * (V - s.Vk) + (-s.gl * (V - s.Vl))) / s.C   
        dhdt = (hinfV - h) / tV
        dsdt = ((1 - s) * sV - k * s)/ tau    #ds equation (3)
        
        
        
        return dVdt, dhdt, dsdt, dtdt
    def __call__(self, X, t):
        N = self.N
        V = X[(0*N):(1*N)]
        h = X[(1*N):(2*N)]
        s = X[(2*N):(3*N)]
        t = X[(3*N):(4*N)]
        dVdt, dhdt = self.rhs(V, h)
        return np.hstack((dVdt.reshape(N, ), dhdt.reshape(N, ))) 
    def integrate(self, tmin, tmax, Nvals=1000):
        N = self.N
        Tvals = np.linspace(tmin, tmax, Nvals)
        Xvals = odeint(self,np.zeros((4*N, )),Tvals) 
        return Xvals, Tvals
    
        
        

N = 64
V = np.zeros((N , 1))
h = np.zeros((N , 1))
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
        
Xvals, Tvals = r.integrate(0, 10)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = "3d")
Vi = Xvals[500, :N]
ax.scatter(r.gNa, r.A.sum(1), Vi)

fig, ax = plt.subplots()
ax.plot(Tvals, Xvals)
plt.show()

# Read about artificial neural network and specifically as they apply to control
# 

  # Check to see if that is the correct "N" for line 44
        # dVdt = [- Inaph - INa - IK - IL - Itonice==0 + Iapp==0 ] / C





