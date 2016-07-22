'''
Created on Jul 19, 2016

@author: andrewkennedy
'''
import numpy as np
import matplotlib.pyplot as plt
from prebotc.utils import findClosest
from prebotc.buteraMultipleCells import Rhs
from progressbar import ProgressBar, Bar, ETA

N = 2
V = np.zeros((N , 1))
h = np.zeros((N , 1))
n = np.zeros((N , 1))
S = np.zeros((N , 1))
# 
# gsynebarVals = [0, 1, 2, 10, 12, 15]
# gsynebarVals = [15, 12, 10, 2, 1, 0]
gsynebarVals = [12, 10]
# gtoniceVals = [0.2, 0.3, 0.4, 0.7, 1.]
gtoniceVals = [0.4, 0.7]
# gsynebarVals = [0, 12]
# gtoniceVals = [0.2, .3]
wastedTime = 30
tmax = wastedTime + 20
fig, Ax = plt.subplots(nrows=len(gsynebarVals), ncols=len(gtoniceVals))

pbar = ProgressBar(widgets=['grid plot ', ETA(), ' ', Bar()],
                   maxval=len(gsynebarVals)*len(gtoniceVals))
pbar.start()
subplotNumber = -1
for k, gsynebar in enumerate(gsynebarVals):
    for j, gtonice in enumerate(gtoniceVals):
        subplotNumber += 1
        pbar.update(subplotNumber)
        r = Rhs(N)
        r.gsyn = gsynebar * (np.ones((N, N)) - np.eye(N))
        r.gt = gtonice
         
        # nominal = r.gNaBar
        # r.gNaBar = np.random.uniform(low = nominal*0.9, high = nominal*1.1, size = (N,))
        EL = -60.0
        r.EL = EL
        X0 = np.hstack([
                        np.zeros(N,),
                        np.zeros(N,),
                        np.zeros(N,),
                        np.zeros(N,),
                        ])
        states, times = r.integrate(X0, 0, tmax, progressBar=False)
        i = findClosest(times, wastedTime)
#         i = 0
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
        fig.suptitle('$E_L = %s [mV]$' % EL)
        ax = Ax[k][j]
        ax.set_title(r'$\bar{g}_{syn,e}=%f$, $g_{tonic,e}=%f$' % (gsynebar, gtonice))
        ax.plot(times[i:], X[i:, 0, :])
        ax.set_ylabel('V [mV]')
        ax.set_xlim((min(times[i:]), max(times[i:])))
        ax.set_xticks([])
        ax.set_yticks([])
    
#         fig.subplots_adjust(hspace=0)
pbar.finish()

fig.savefig('../doc/buteraB_grid-EL%s.png' % EL)
plt.show()