from LIF import *

h = 0.02
tau = 1e-3

# fig 2 top
VF = 2
Vmin = -4
VR = 1
a0 = 1
a1 = 0
b = 1.5

v0 = 0
sigma0 = 0.25

N_stat1 = 0.1924
N_stat2 = 2.319

p_init = gaussian_initial_cond(v0,sigma0,VF,Vmin,h)
p_stat1,_ = stationnary_initial_cond(N_stat1,a0,b,VF,Vmin,VR,h)
p_stat2,_ = stationnary_initial_cond(N_stat2,a0,b,VF,Vmin,VR,h)

p2,Nl,_,traj = fokker_plank_solve(p_init,VF,Vmin,VR,a0,a1,b,h,tau,5,return_traj = True)
time = np.linspace(0,5,Nl.shape[0])

entro1 = np.empty(traj.shape[0])
for i in range(traj.shape[0]):
    entro1[i] = entropy(traj[i,:],p_stat1,h)

entro2 = np.empty(traj.shape[0])
for i in range(traj.shape[0]):
    entro2[i] = entropy(traj[i,:],p_stat2,h)

if True :
    fig,axs = plt.subplots(1,3,layout="constrained",figsize=(15,5))
    axs[0].plot(v[1:-1],p_init,label="t=0")
    axs[0].plot(v[1:-1],traj[100,:],label="t=1")
    axs[0].plot(v[1:-1],traj[400,:],label="t=4")
    axs[0].plot(v[1:-1],p_stat1,"--k",label="p_stationary stable")
    axs[0].plot(v[1:-1],p_stat2,":k",label="p_stationary unstable")
    axs[0].legend()
    axs[0].set_xlabel("v")
    axs[0].set_ylabel("p(v,t)")
    axs[1].plot(time,entro1)
    axs[1].set_xlabel("t")
    axs[1].set_ylabel("S(t)")
    axs[1].set_title("Entropy w.r.t stable distribution")
    axs[2].plot(time,entro2)
    axs[2].set_xlabel("t")
    axs[2].set_ylabel("S(t)")
    axs[2].set_title("Entropy w.r.t unstable distribution")
    fig.savefig("figs/fig5.pdf")
    fig.savefig("figs/fig5.png")

plt.show()
