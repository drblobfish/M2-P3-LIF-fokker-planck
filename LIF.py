import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

def integral(p,h):
    return np.sum(p)*h

def gaussian_initial_cond(v0,sigma0,VF,Vmax,h):
    N = int((VF-Vmin)/h)
    v = np.linspace(Vmin,VF,N+1)
    p = np.exp(-(v[1:N]-v0)**2/sigma0)
    p /= integral(p,h)
    return p

def fokker_plank_solve(p0,VF,Vmin,VR,a0,a1,b,h,tau,T):
    N = int((VF-Vmin)/h)
    nb_iter = int(T/tau)
    VR_index = int(np.round((VR-Vmin)*N/(VF-Vmin)))
    v = np.linspace(Vmin,VF,N+1)
    p = p0.copy()

    Nhl = np.empty(nb_iter)

    # at each step
    for i in range(nb_iter):
        Nh = a0 * p[N-2] / (h - a1 * p[N-2])
        Nhl[i] = Nh
        aNh = a0 + a1 * Nh
        M = np.exp(- (v - b * Nh)**2 / (2*(a0 + a1 * Nh)))
        Mmidpoints = 1/(0.5 * (1/M[:-1] + 1/M[1:]))

        diagonal = np.zeros(N-1)
        diagonal[1:N-2] = 1 + (tau/h**2) * aNh * (Mmidpoints[1:N-2] + Mmidpoints[2:N-1])/M[2:N-1]
        diagonal[0] = 1 + (tau/h**2) * aNh * Mmidpoints[1]/M[1]
        diagonal[N-2] = 1 + (tau/h**2) * aNh * Mmidpoints[N-2]/M[N-1]

        subdiag = - (tau/h**2) * aNh * Mmidpoints[1:N-1]/M[1:N-1]
        surdiag = - (tau/h**2) * aNh * Mmidpoints[1:N-1]/M[2:N]

        p[N-2] -= (tau/h) * Nh
        p[VR_index] += (tau/h) * Nh

        sparray = sp.diags_array([subdiag,diagonal,surdiag],offsets=[-1,0,1],format="csr")
        p = spla.spsolve(sparray,p)

    return p,Nhl

if __name__=="__main__":
    VF = 2
    Vmin = -2
    VR = 1
    a0 = 1
    a1 = 0
    b = 1.5

    h = 6/384
    tau = 1/20_000
    T = 0.0405

    v0 = 1.5
    sigma0 = 0.005
    N = int((VF-Vmin)/h)
    p0 = gaussian_initial_cond(v0,sigma0,VF,Vmin,h)
    p = p0.copy()

    plt.plot(np.linspace(Vmin,VF,N+1)[1:N],p)
    plt.show()

    p,Nhl = fokker_plank_solve(p,VF,Vmin,VR,a0,a1,b,h,tau,T)
    print(integral(p,h))

    plt.plot(np.linspace(Vmin,VF,N+1)[1:N],p)
    plt.show()
    plt.plot(Nhl)
    plt.show()
