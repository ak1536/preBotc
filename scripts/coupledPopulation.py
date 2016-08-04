'''
Created on Jul 19, 2016

@author: andrewkennedy
'''
import numpy as np
import matplotlib.pyplot as plt
from prebotc.utils import findClosest
from prebotc.buteraMultipleCells import Rhs
import matplotlib.animation as animation
from progressbar import ProgressBar, ETA, Bar

N = 12
ntonic = 1
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
# Cells are more excitable when this gets closer to 0 (less negative).
r.EL = np.random.uniform(low=-61.5, high=-61, size=(N,))

def PbarETA(label, maxval=1):
    return ProgressBar(maxval=maxval, widgets=[label, ' | ', ETA(), Bar()]).start()
    


X0 = np.hstack([
                np.zeros(N,),
                np.zeros(N,),
                np.zeros(N,),
                np.zeros(N,),
                ])

states, times = r.integrate(X0, 0, tmax, progressBar=True, rtol=1e-3, atol=1e-3)
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

fig.savefig('../doc/buteraB_population-N%d.png' % (N,))


## Make a movie.
fig, ax = plt.subplots()

fps = 15
framenumbers = range(int((max(times) - min(times))*fps))

# Find the indices into the trajectory corresponding to 
frameTimeIndices = []
lasti = 0
for t in PbarETA('Finding frame times.')(np.linspace(min(times), max(times), len(framenumbers))):
    nexti = lasti + findClosest(times[lasti:], t)
    frameTimeIndices.append(nexti)
    lasti = frameTimeIndices[-1]
              

pbar = PbarETA('Animating.', maxval=max(framenumbers))
def updateFrame(frame):
    pbar.update(frame)
    frame = frameTimeIndices[frame]
    ax.cla()
    h = X[:, 1, :]
    ax.scatter(r.EL, h[frame])
    ax.set_ylim(h.min(), h.max())
    ax.set_xlabel('$E_L$ [mV]')
    ax.set_ylabel('$h$')
    ax.set_title('$t=%.1f$ [s]' % times[frame])
anim = animation.FuncAnimation(fig, updateFrame, framenumbers)
Writer = animation.writers['ffmpeg']
writer = Writer(fps=fps, metadata=dict(artist='Tom Bertalan'), bitrate=1800)
anim.save('ELhet-N%d-%dtonic.mp4' % (N, ntonic), writer=writer)
pbar.finish()


## Show plots and animation.
#plt.show()
