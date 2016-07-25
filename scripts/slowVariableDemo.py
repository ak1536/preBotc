import numpy as np
import matplotlib.pyplot as plt
from prebotc.utils import findClosest
from prebotc.buteraSlowVariable import Rhs

N = 2
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
EL = -60.0
r.EL = EL
X0 = np.hstack([
                np.zeros(N,),
                np.zeros(N,),
                np.zeros(N,),
                np.zeros(N,),
                np.ones(N,) * 0.15,
                ])
states, times = r.integrate(X0, 0, 32, progressBar='IVP ')
i = findClosest(times, 20)
X = states.reshape(times.size, 5, N)
    
    
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
fig1, Ax = plt.subplots(nrows=5)
fig1.suptitle('$E_L = %s [mV]$' % EL)
for k in range(5):
    ax = Ax[k]
    ax.plot(times[i:], X[i:, k, :])
    ax.set_ylabel(['V [mV]', 'h', 'n', 'S', 'A'][k])
    for a in Ax:
        a.set_xlim((min(times[i:]), max(times[i:])))
    if k != 3:
        ax.set_xticks([])
fig1.subplots_adjust(hspace=0)

fig2, ax = plt.subplots()
ax.plot(X[i:, 1, :], X[i:, 4, :])
ax.set_ylabel('A')
ax.set_xlabel('h')

fig1.savefig('../doc/buteraB-network_VhnSA-EL%s.png' % EL)
fig2.savefig('../doc/buteraB-network_h_A_ring-EL%s.png' % EL)
plt.show()