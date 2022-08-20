
import numpy as np
from numba import jit
from tqdm import tqdm
from IPython.display import clear_output

# Compile with jit for faster matmul. WARNING: Jit is extremely temperamental
@jit(nopython=True)
def matrix_deriv(u_arr, dx, D_const, N):

    derivmat = np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            if i!=j:
                derivmat[i,j] = -(2*(np.pi**2)/(N**2 * dx**2))*(-1)**(i-j)/np.sin(np.pi*(i-j)/N)**2
            else:
                derivmat[i,j] = -((np.pi**2)/(3 * dx**2))*(1+2/N**2)
            
    return D_const*derivmat@u_arr

def forward_euler(x_range, f_IC, dt_FE, T_FE):
    
    """
    Performs forward Euler numerical integration for the 1D diffusion IC-BVP problem.

    Parameters
    ----------
    x_rng : ndarray
            N linearly spaced values on the interval [0,L] with shape (N,).
    f_IC  : ndarray
            Function defining initial condition, with shape (N,).
    dt_FE : float
            Time step increment to be used for integration.
    T_FE  : float
            Max time over which to integrate.
    Returns
    -------
    u_FE  : ndarray
            Numerical solution, with shape (N,M), where M = int(np.ceil(T_FE/dt_FE)).
    """
    