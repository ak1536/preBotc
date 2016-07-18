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


class Rhs(object):
    def __init__(self, N):
        s = self
        s.N = N
        s.A = np.zeros((N, N))
        s.ENa = 50.0          # mV
        s.gkBar = 11.2      # nS
        s.Ek = -85.0        # mV
        s.gNaPBar = 2.8     # nS  # TODO: Is this or gNaBar a typo in the paper?
        s.gNaBar = 28.0     # nS
        s.gNaP = 2.8        # nS
        s.gNa = 28.0        # nS
        s.gK = 11.2         # nS
        s.gL = 2.8          # nS
        s.gt = 0.35         # nS
        s.tauh = 10.0       # s
        s.taun = 0.01       # s
        s.taus = 0.005      # s 
        s.thetamp = -40.0   # mV
        s.thetam = -34.0    # mV
        s.thetan = -29.0    # mV
        s.thetah = -48.0    # mV
        s.thetas = -10.0    # mV
        s.sigmamp = -6.0    # mV
        s.sigmam = -5.0     # mV
        s.sigman = -4.0     # mV
        s.sigmah = 6.0      # mV
        s.sigmas = -5.0     # mV
        s.gLBar = 2.8       # nS
        s.EL = -65.0        # mV 
        s.Iapp = 0.0        # pA; Zero unless otherwise noted
        s.C = 21.0          # pF
        s.tauBarh = 10.0     # s
        s.tauBarn = 0.010   # s
        s.thetah = -48.0    # mV
        s.thetan = -29.0    # mV
        s.thetamNaP = -40.0 # mV
        s.thetamNa = -34.0    # mV
        s.sigmah = 6.0      # mV
        s.sigman = -4.0       # mV
        s.sigmamNaP = -6.0  # mV
        s.sigmamNa = -5.0   # mV
        s.alpha = 3.0       # 1/s
        s.beta = 0.1        # 1/mV
        s.thetaApce = -30.0   # mV
        
        
        
    def rhs(self, V, h, n, Apce):
        s = self
#         tauhV = s.tauBarh / (np.cosh((V - s.thetah) / 2 / s.sigmah))     
#         taunV = s.tauBarn / (np.cosh((V - s.thetan) / 2 / s.sigman)) 
#         hinfV = expit((s.thetah - V) / s.sigmah)
#         ninfV = expit((s.thetan - V) / s.sigman)
#         minfNaPV = expit((s.thetamNaP - V) / s.sigmamNaP)
#         minfNaV = expit((s.thetamNa - V) / s.sigmamNa)
        gV = expit(s.beta *(V - s.thetaApce))
        betah = (1 /(2 * s.tauh)) *  np.exp((V - s.thetah)/(2 * s.sigmah))
        alphah = (1 /(2 * s.tauh)) *  np.exp(-1 *(V - s.thetah)/(2 * s.sigmah))
        betan = (1 /(2 * s.taun)) *  np.exp((V - s.thetan)/(2 * s.sigman))
        alphan = (1 /(2 * s.taun)) *  np.exp(-1 *(V - s.thetan)/(2 * s.sigman))
        mbarp = expit((s.thetamp - V)/ s.sigmamp)
        mbar = expit((s.thetam - V)/ s.sigmam)
#         sbar = expit((s.thetas - V)/ s.sigmas)
   
        # m n h need expit functions
        # V h n are the function parameters
        # dVdt = [- Inap - INa - IK - IL - gtV] / C
        dVdt = ( 
                - s.gNaP * mbarp * h * (V - s.ENa)         # -I_{NaP} = I_{NaP-h} for model 1. (\neq I_{Na})
                - s.gNa * mbar**3 * (1 - n) * (V - s.ENa)  # -I_{Na}
                - s.gK * n**4 * (V - s.Ek)                    # -I_K
                - s.gL * (V - s.EL)
                - s.gt * V                           # -I_L
                ) / s.C  * 1000 # Conversion to mV. Naturally in terms of V/s without this.
        # dhdt = alphah*(1 - h) - betah*h
        dhdt = alphah * (1 - h) - betah * h 
        dndt = alphan * (1 - n) - betan * n
        dApcedt = s.alpha * (gV - Apce)
        return dVdt, dhdt, dndt, dApcedt
    
    def __call__(self, X, t):
        N = self.N
        V = X[(0*N):(1*N)]
        h = X[(1*N):(2*N)]
        n = X[(2*N):(3*N)]
        Apce = X[(3*N):(4*N)]
        dVdt, dhdt, dndt, dApcedt = self.rhs(V, h, n, Apce)
        out = np.hstack((dVdt.reshape(N, ), dhdt.reshape(N, ), dndt.reshape(N, ), dApcedt.reshape(N, )))
        return out
    
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
    ax = Ax[0]
    bx = Ax[1]
    cx = Ax[2]
    dx = Ax[3]
    ax.plot(times[i:], X[i:, 0, :])
    ax.set_ylabel('MilliVolts')
    ax.set_title('$E_L = %s [mV]$' % EL)
    for a in Ax:
        a.set_xlim((min(times[i:]), max(times[i:])))
    ax.set_xticks([])
    ax.set_ylim(-65, 10)

#     f, (ax, bx, cx, dx) = plt.subplots(4, sharex=True, sharey=False)
#     ax.plot(times[i:], X[i:, 0, :])
#     ax.set_ylim(-65, 10)
#     bx.plot(times[i:], X[i:, 2, :])
#     bx.set_ylim(0, 1)
#     cx.plot(times[i:], X[i:, 1, :])
#     cx.set_ylim(0.525, 0.55)
#     dx.plot(times[i:], X[i:, 3, :])
#     dx.set_ylim(0.10, 0.135)
#     # Fine-tune figure; make subplots close to each other and hide x ticks for
#     # all but bottom plot.
#     f.subplots_adjust(hspace=0)
#     plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
# 
#     
#     bx.plot(times[i:], X[i:, -1, :])
#     bx.set_ylabel('A')
#     bx.set_xlabel('Time (s)')
#     
    fig.subplots_adjust(hspace=0)
    fig, ax = plt.subplots()
    ax.plot(X[i:, 1, :], X[i:, 3, :])
    ax.set_ylabel('A')
    ax.set_xlabel('h')
#     
#     fig, ax = plt.subplots()
#     ax.plot(times[i:], X[i:, 0, :])
#     ax.set_ylabel('V')
#     ax.set_xlabel('t')
#     
#     fig, ax = plt.subplots()
#     ax.plot(times[i:], X[i:, 1, :])
#     ax.set_ylabel('h')
#     ax.set_xlabel('t')
#     
#     fig, ax = plt.subplots()
#     ax.plot(times[i:], X[i:, 2, :])
#     ax.set_ylabel('n')
#     ax.set_xlabel('t')
#     
#     fig, ax = plt.subplots()
#     ax.plot(times[i:], X[i:, 3, :])
#     ax.set_ylabel('A')
#     ax.set_xlabel('t')
#     
    
#     fig.savefig('buteraA_bursts-EL%s.png' % EL)
    plt.show()
        
        
        
        
        
        
