import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def integral(p,h):
    return np.sum(p)*h

def gaussian_initial_cond(v0,sigma0,VF,Vmin,h,N=None):
    if N == None:
        N = int(np.round((VF-Vmin)/h))
    else :
        h = (VF-Vmin)/N
    v = np.linspace(Vmin,VF,N+1)
    p = np.exp(-(v[1:N]-v0)**2/sigma0)
    p /= integral(p,h)
    return p

def stationnary_initial_cond(N_inf,a0,b,VF,Vmin,VR,h,N=None):
    if N == None:
        N = int(np.round((VF-Vmin)/h))
    else :
        h = (VF-Vmin)/N
    v = np.linspace(Vmin,VF,N+1)
    VR_index = int(np.round((VR-Vmin)*N/(VF-Vmin)))
    max_v_VR = v.copy()
    max_v_VR[:VR_index+1] = VR
    return 2 * N_inf * np.exp(0.5 * (b * N_inf - (v - b * N_inf)**2)/a0) * (np.exp(-0.5 * max_v_VR/a0) - np.exp(-0.5 * VF/a0)),v

p,v = stationnary_initial_cond(0.1924,1,1.5,2,-4,1,0.02)
p.shape
v.shape
plt.plot(v,p)
plt.show()

def fokker_plank_solve(p0,VF,Vmin,VR,a0,a1,b,h,tau,T,N=None,nb_iter=None,harmonic_mean = False):
    if N == None:
        N = int(np.round((VF-Vmin)/h))
    else :
        h = (VF-Vmin)/N
    if nb_iter == None:
        nb_iter = int(T/tau)
    else :
        tau = T/nb_iter
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
        if harmonic_mean :
            Mmidpoints = 1/(0.5 * (1/M[:-1] + 1/M[1:]))
        else :
            v_midpoint = 0.5 * (v[:-1] + v[1:])
            Mmidpoints  = np.exp(- (v_midpoint - b * Nh)**2 / (2*(a0 + a1 * Nh)))

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

    return p,Nhl,v

class UpdateDist:
    def __init__(self,fig,ax,p,VF,Vmin,VR,a0,a1,b,h,tau,T,nb_step):
        self.fig = fig
        self.ax = ax
        self.p = np.zeros(N+1)
        self.p[1:N] = p
        self.VF = VF
        self.Vmin = Vmin
        self.VR = VR
        self.a0 = a0
        self.a1 = a1
        self.b = b
        self.h = h
        self.tau = tau
        self.T = T
        self.nb_step = nb_step
        self.N = int((VF-Vmin)/h)
        self.v = np.linspace(Vmin,VF,N+1)
        self.ln, = ax.plot(self.v, self.p)
        self.ax.set_ylim(0,2)
        self.ax.set(title = f"a = {a0}; b = {b}; VF = {VF}; VR = {VR}")
        self.time_text = ax.text(0.05, 0.9, 't = 0', transform=ax.transAxes)
    def start(self):
        return self.ln,self.time_text
    def __call__(self,frame):
        self.p[1:N],_,_ = fokker_plank_solve(self.p[1:N],self.VF,self.Vmin,self.VR,self.a0,self.a1,self.b,self.h,self.tau,self.T/self.nb_step)
        self.ln.set_data(self.v, self.p)
        self.time_text.set_text(f"t = {frame:.3f}")
        return self.ln,self.time_text

if __name__=="__main__":
    VF = 2
    Vmin = 0
    VR = 1
    a0 = 1
    a1 = 0
    b = 1.5

    h = 0.02
    tau = 1e-4
    T = 0.0405

    v0 = 1.5
    sigma0 = 0.005
    N = int((VF-Vmin)/h)
    p0 = gaussian_initial_cond(v0,sigma0,VF,Vmin,h)
    p = p0.copy()

    # plt.plot(np.linspace(Vmin,VF,N+1)[1:N],p)
    # plt.show()

    # p,Nhl,v = fokker_plank_solve(p,VF,Vmin,VR,a0,a1,b,h,tau,T)
    # p_full = np.empty(N+1)
    # p_full[1:N] = p
    # plt.plot(v,p_full)
    # plt.show()
    # plt.plot(Nhl)
    # plt.show()

    fig, ax = plt.subplots()
    nb_step = 50
    ud = UpdateDist(fig,ax,p,VF,Vmin,VR,a0,a1,b,h,tau,T,nb_step)
    ani = FuncAnimation(fig, ud,init_func=ud.start, frames=np.linspace(0, T, nb_step), blit=True,repeat= False)
    plt.show()
    # ani.save("movie.mp4")
