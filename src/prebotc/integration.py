'''A wrapper function for integrating systems of ODEs.

Created on Jun 22, 2014

@author: bertalan@princeton.edu
'''
import logging
import numpy as np
from scipy.integrate import ode


def integrate(initial, dxdt, tmin=0, tmax=800, giveTime=False, minSteps=1000,
              progressBar=None, backward=False, backend='vode',  # dop853, dopri5, vode
              **kwargs
              ):
    """Solve an initial value problem.
 
    Based on http://stackoverflow.com/a/12935589/1224886
    
    Parameters
    ==========
    initial : array_like
        An array of size N indicating the inital values.
 
    dxdt : function
        Parameters
        ----------
        X : array_like
            state, array of size N
        t : float
            time; required only if giveTime==True
            
        Returns
        -------
        dXdt : array_like
            ODEs evaluated at (X, t).
 
    tmin : float, optional
        starting time
 
    tmax : float, optional
        ending time
 
    giveTime : bool, optional
        Whether the RHS function dxdt should take a second time parameter.
        For nonautonomous systems.
 
    minSteps : int, optional
        Minimum number of time steps to take. More intermediate steps might be taken,
        but not returned.
 
    progressBar : bool, optional
        Whether to display a progress bar (True or str), or not (False),
        or let this depend on whether the loglevel is below INFO (None, default).
        If a string is provided, it will be included in the bar.
 
    backward : bool, optional
        Whether to integrate backward, from tmax to tmin.
        nonautonomous systems should still be treated correctly, being given
        the correct time.
        Results will be in order of decreasing t.
        
    backend : str, optional.
        See documentation for scipy.integrate.ode.
 
    **kwargs
        Other keyword arguments are passed on to scipy,integrate.ode.set_integrator;
        See the documentation for scipy.integrate.ode.
 
    Returns
    =======
    X : ndarray
        the trajectory, shape (nstepSaved, N)
    T : ndarray, shape (nstepSaved,)
        the times for each trajectory point
    nstepSaved will probably be minSteps,
    but there might be a few extra rows to guarantee that the last time
    in the T array is >= tmax.
 
 
    >>> def diff(X):
    ...     x = X[0]
    ...     y = X[1]
    ...     return np.array([x ** 2 / y ** 2, x ** 1 / y ** 3])
    ...
    >>> X, T = integrate([0, 1], diff, minSteps=100000)
    >>> X.shape
    (100000, 2)
    >>> T.shape
    (100000,)
    >>> print X[-1,:]
    [ 0.  1.]
    >>> assert abs(800.0 - T[-1]) < 0.1
    """
    
    
    ## SET UP THE ANNOTATED RHS FUNCTION.
    # Handle autonomous and nonautonomous differnetly, for convenience in the former case.
    if giveTime:
        def dxdtTimed(t, y):
            return dxdt(y, t)
    else:
        def dxdtTimed(t, y):
            return dxdt(y)
    # If backwards integration is called for, multiple RHS output by -1.
    if backward:
        # TODO: Unittest backwards integration.
        def dxdtSigned(t, y):
            return -dxdtTimed(t, y)
    else:
        def dxdtSigned(t, y):
            return dxdtTimed(t, y)
    
    
    ## SET UP THE SOLVER OBJECT
    # The solver should take at least minSteps steps.
    maximumDt = float(tmax - tmin) / minSteps
    solver = ode(dxdtSigned).set_integrator(backend, **kwargs)
    
    
    ## SET UP PROGRESSBAR.
    # If the loglevel wants at least as much output as INFO, we'll add a progress bar.
    logger = logging.getLogger(__name__)
    if logger.getEffectiveLevel() <= logging.INFO:
        if progressBar is None:
            progressBar = True
    else:
        if progressBar is None:
            progressBar = False
    if progressBar:
        from progressbar import ProgressBar, Bar, ETA
        if isinstance(progressBar, str):
            barLabel = progressBar
        else:
            barLabel = 'IVP '
        pbar = ProgressBar(maxval=(tmax-tmin),
                           widgets=[barLabel, Bar(), ETA()])
        pbar.start()
        def updatePbar(t):
            if t <= tmax:
                pbar.update(t - tmin)
            pbar.widgets[0] = '%s (t=%f) ' % (barLabel.strip(), t)
        finishPbar = lambda : pbar.finish()
    else:
        updatePbar = lambda t : None
        finishPbar = lambda : None
        
    
    ## DO THE INTEGRATION.
    solver.set_initial_value(initial, tmin)
    # Unlike scipy.odeint, the ode solvers do not return trajectories,
    # but instead return a final point.
    # Solvers like dopri5 r dop853 will accept a solout callback function
    # which can be used to collect all the intermediate steps taken between
    # calls to .integrate. But they don't appear to handle stiff problems well,
    # and claim that our problems are stiff.
    # So, we store history in lists, to be compacted to arrays upon return.
    T = []
    X = []
    while solver.successful() and solver.t < tmax:
        solver.integrate(solver.t + maximumDt, step=True)
        t = solver.t
        updatePbar(t)
        T.append(t)
        X.append(solver.y)
    if solver.t >= tmax:
        finishPbar()
    
    return np.array(X), np.array(T)
