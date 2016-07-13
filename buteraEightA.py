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
        s.N = N
        s.A = np.zeros((N, N))
        s.ENa = 50          # mV
        s.gkBar = 11.2      # nS
        s.Ek = -85.0        # mV
        s.gNaPBar = 2.8     # nS
        s.gNaBar = 28.0     # nS
        s.gLBar = 2.8       # nS
        s.EL = -65.0        # mV 
        s.Iapp = 10.0        # pA; Zero unless otherwise noted
        s.C = 21.0          # pF
        s.tauBarh = 10.     # s
        s.tauBarn = 0.010   # s
        s.thetah = -48.0    # mV
        s.thetan = -29.0    # mV
        s.thetamNaP = -40.0 # mV
        s.thetamNa = -34    # mV
        s.sigmah = 6.0      # mV
        s.sigman = -4       # mV
        s.sigmamNaP = -6.0  # mV
        s.sigmamNa = -5.0   # mV
        
   
    def rhs(self, V, h, n):
        s = self
        tauhV = s.tauBarh / (np.cosh((V - s.thetah) / 2 / s.sigmah))     
        taunV = s.tauBarn / (np.cosh((V - s.thetan) / 2 / s.sigman)) 
        hinfV = expit((s.thetah - V) / s.sigmah)
        minfNaPV = expit((s.thetamNaP - V) / s.sigmamNaP)
        minfNaV = expit((s.thetamNa - V) / s.sigmamNa)
        ninfV = expit((s.thetan - V) / s.sigman)
        
        # m n h need expit functions
        # V h n are the function parameters
        # dVdt = [- Inaph - INa - IK - IL - Itonice==0 + Iapp==0 ] / C
        dVdt = ( 
                - s.gNaPBar * minfNaPV * h * (V - s.ENa)         # -I_{NaP} = I_{NaP-h} for model 1. (\neq I_{Na})
                - s.gNaBar * minfNaV**3 * (1 - n) * (V - s.ENa)  # -I_{Na}
                - s.gkBar * n**4 * (V - s.Ek)                    # -I_K
                - s.gLBar * (V - s.EL)                           # -I_L
                + s.Iapp
                ) / s.C  * 1000 # Conversion to mV. Naturally in terms of V/s without this.
        dhdt = (hinfV - h) / tauhV 
        dndt = (ninfV - n) / taunV
        return dVdt, dhdt, dndt
    
    def __call__(self, X, t):
        N = self.N
        V = X[(0*N):(1*N)]
        h = X[(1*N):(2*N)]
        n = X[(2*N):(3*N)]
        dVdt, dhdt, dndt = self.rhs(V, h, n)
        out = np.hstack((dVdt.reshape(N, ), dhdt.reshape(N, ), dndt.reshape(N, )))
        return out 
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
nominal = r.gNaBar
r.gNaBar = np.random.uniform(low = nominal*0.9, high = nominal*1.1, size = (N,))
        
Xvals, Tvals = r.integrate(0, 20)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = "3d")
X = Xvals.reshape(Tvals.size, 3, N)
Vi = X[500, 0, :]
# Vi = Xvals[500, :N]
ax.scatter(r.gNaBar, r.A.sum(1), Vi)
ax.set_xlabel('gNaBar [nS]')
ax.set_ylabel('degree')
ax.set_zlabel('V [mV]')

fig, ax = plt.subplots()
ax.plot(Tvals, X[:, 0, :])
ax.set_xlabel('Time')
ax.set_ylabel('MilliVolts')
plt.show()







