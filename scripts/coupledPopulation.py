'''
Created on Jul 19, 2016

@author: andrewkennedy
'''
import numpy as np
import matplotlib.pyplot as plt
from prebotc.utils import findClosest
from prebotc.buteraMultipleCells import Rhs
from progressbar import ProgressBar, Bar, ETA

N = 12
ntonic = 4
assert ntonic < N
V = np.zeros((N, 1))
h = np.zeros((N, 1))
n = np.zeros((N, 1))
S = np.zeros((N, 1))
# 
wastedTime = 30
tmax = wastedTime + 20

gtonice = [0] * (N - ntonic)
gtonice.extend([2.0]*ntonic)  # Just ntonic cells get tonic input.
gsynebar = 1  # Strength of connections in the network.
     
r = Rhs(N)
r.gsyn = gsynebar * (np.ones((N, N)) - np.eye(N))  # All-all connection without self-loops.
r.gt = np.array(gtonice)
r.EL = -60.0  # Cells are more excitable when this gets closer to 0 (less negative).

X0 = np.hstack([
                np.zeros(N,),
                np.zeros(N,),
                np.zeros(N,),
                np.zeros(N,),
                ])

states, times = r.integrate(X0, 0, tmax, progressBar=True)
i = findClosest(times, wastedTime)
X = states.reshape(times.size, 4, N) # variables are V, h, n, s
    
    
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
fig, axes = plt.subplots(nrows=2)
fig.suptitle('$E_L = %s [mV]$' % r.EL)
fig.suptitle(r'$\bar{g}_{syn,e}=%f$, $g_{tonic,e}=%f$ for last %d of %d cells' % (gsynebar, gtonice[-1], ntonic , N))

ax = axes[0]
ax.plot(times[i:], X[i:, 0, :], color='black')  # Plot all the voltages.
ax.set_ylabel('V [mV]')
ax.set_xlim((min(times[i:]), max(times[i:])))

ax = axes[1]
ax.set_title(r'$\bar{g}_{syn,e}=%f$, $g_{tonic,e}=%f$ for last %d of %d cells' % (gsynebar, gtonice[-1], ntonic , N))
ax.plot(times[i:], X[i:, 1, :], color='red')
ax.set_ylabel('$h$')
ax.set_xlim((min(times[i:]), max(times[i:])))
    

fig.savefig('../doc/buteraB_population-EL%s-N%d.png' % (r.EL, N))
plt.show()
