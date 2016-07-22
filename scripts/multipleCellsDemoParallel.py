'''
Created on Jul 19, 2016

@author: andrewkennedy
'''
import numpy as np
import matplotlib.pyplot as plt
from prebotc.utils import findClosest
from prebotc.buteraMultipleCells import Rhs
from progressbar import ProgressBar, Bar, ETA
from joblib import Parallel, delayed

EL = -60.0

def job(datum):
    gsynebar, gtonice, k, j = datum
    print datum
    
    N = 2
    wastedTime = 60
    tmax = wastedTime + 15
    r = Rhs(N)
    r.gsyn = gsynebar * (np.ones((N, N)) - np.eye(N))
    r.gt = gtonice
     
    # nominal = r.gNaBar
    # r.gNaBar = np.random.uniform(low = nominal*0.9, high = nominal*1.1, size = (N,))
    r.EL = EL
    X0 = np.hstack([
                    np.zeros(N,),
                    np.zeros(N,),
                    np.zeros(N,),
                    np.zeros(N,),
                    ])
    states, times = r.integrate(X0, 0, tmax, progressBar=True)
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
    
    return times, X, datum, i

    
def main():
    # gsynebarVals = [0, 1, 2, 10, 12, 15]
    # gsynebarVals = [15, 12, 10, 2, 1, 0]
    gsynebarVals = [12, 10]
    # gtoniceVals = [0.2, 0.3, 0.4, 0.7, 1.]
    gtoniceVals = [0.4, 0.7]
    # gsynebarVals = [0, 12]
    # gtoniceVals = [0.2, .3]
    
    fig, Ax = plt.subplots(nrows=len(gsynebarVals), ncols=len(gtoniceVals))
    
    data = []
    for k, x in enumerate(gsynebarVals):
        for j, y in enumerate(gtoniceVals):
            data.append((x, y, k, j))
    
    result = Parallel(n_jobs=-1)(delayed(job)(datum) for datum in data)
    
    def plot(times, X, datum, i):
        # Plot voltage trajectories.
        gsynebar, gtonice, k, j = datum
        fig.suptitle('$E_L = %s [mV]$' % EL)
        ax = Ax[k][j]
        ax.set_title(r'$\bar{g}_{syn,e}=%f$, $g_{tonic,e}=%f$' % (gsynebar, gtonice))
        ax.plot(times[i:], X[i:, 0, :])
        ax.set_ylabel('V [mV]')
        ax.set_xlim((min(times[i:]), max(times[i:])))
        ax.set_xticks([])
        ax.set_yticks([])
        
    for res in result:
        plot(*res)
    
    fig.savefig('../doc/buteraB_grid_Parallel-EL%s.png' % EL)
    plt.show()
    
if __name__ == '__main__':
    main()
    