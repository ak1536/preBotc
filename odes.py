import numpy as np
from scipy.integrate import odeint

class Rhs(object):
    def __init__(self, N):
        s = self
        s.N = N
        s.gsyn = 0.3
        s.Vsyn = 0
        s.A = np.zeros((N , N))
        s.gNa = 3.0
        s.VNa = 50.0 
        s.gl = 2.4
        s.Vl = -65.0
        s.Iapp = 24.0
        s.C = 0.21
        s.eps = 0.1
            
    def rhs(self, V , h):
        s = self
        sV = 1 / (1 + np.exp( - (V+40) / 5))
        tV = 1 / (s.eps * np.cosh((V + 44) / 12))
        hinfV = 1 / (1 + np.exp((V + 44) / 6))
        mV = 1 / (1 + np.exp( - (V + 37) / 6))
        
        Isyn = ((s.gsyn * (s.Vsyn - V)) / s.N) * np.dot(s.A , sV)
        dVdt = (-s.gNa * mV * h * (V - s.VNa) - s.gl * (V - s.Vl) + Isyn + s.Iapp) / s.C
        dhdt = (hinfV - h) / tV
        return dVdt, dhdt
    def __call__(self, X, t):
        N = self.N
        V = X[:N]
        h = X[N:]
        dVdt, dhdt = self.rhs(V, h)
        return np.hstack((dVdt.reshape(N, ), dhdt.reshape(N, ))) 
    def integrate(self, tmin, tmax, Nvals=1000):
        N = self.N
        Tvals = np.linspace(tmin, tmax, Nvals)
        Xvals = odeint(self,np.zeros((2*N, )),Tvals) 
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





