 
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
from buteraEightA import Rhs as ButeraRHS


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

  
if __name__ == '__main__':
            
    N = 1
    V = np.zeros((N , 1))
    h = np.zeros((N , 1))
    n = np.zeros((N , 1))
    Apce = np.zeros((N , 1))
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
    EL = -65.0
    r.EL = EL
    X0 = np.hstack([
                    np.zeros(N,),
                    np.zeros(N,),
                    np.zeros(N,),
                    np.ones(N,) * 0.15,
                    ])
    states, times = r.integrate(X0, 0, 32, progressBar='IVP ')
    i = findClosest(times, 20)
    X = states.reshape(times.size, 4, N)
        
        
    # # Plot PCE-fittable surface.
    # from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1, projection = "3d")
    # Vi = X[500, 0, :]
    # # Vi = states[500, :N]
    # ax.scatter(r.gNaBar, r.A.sum(1), Vi)
    # ax.set_xlabel('gNaBar [nS]')
    # ax.set_ylabel('degree')
    # ax.set_zlabel('V [mV]')
    
    
    # Plot voltage trajectories.
    fig, Ax = plt.subplots(nrows=4)
    fig.suptitle('$E_L = %s [mV]$' % EL)
    for k in range(4):
        ax = Ax[k]
        ax.plot(times[i:], X[i:, k, :])
        ax.set_ylabel(['V [mV]', 'h', 'n', 'A'][k])
        for a in Ax:
            a.set_xlim((min(times[i:]), max(times[i:])))
        if k != 3:
            ax.set_xticks([])
    
    fig.subplots_adjust(hspace=0)
    
    fig, ax = plt.subplots()
    ax.plot(X[i:, 1, :], X[i:, 3, :])
    ax.set_ylabel('A')
    ax.set_xlabel('h')
    
    fig.savefig('buteraA_bursts-EL%s.png' % EL)
    plt.show()
        
        
        
        
        
        

        
        
        
        
