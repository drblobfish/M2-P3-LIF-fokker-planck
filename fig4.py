from LIF import *

h = 0.02
tau = 1e-3

# fig 2 top
VF = 2
Vmin = -4
VR = 1
a0 = 1
a1 = 0
b = 0

v0 = 0
sigma0 = 0.25

N_stat = 0.1377

p_init = gaussian_initial_cond(v0,sigma0,VF,Vmin,h)
p_stat,_ = stationnary_initial_cond(N_stat,a0,b,VF,Vmin,VR,h)

p2,Nl,_,traj = fokker_plank_solve(p_init,VF,Vmin,VR,a0,a1,b,h,tau,5,return_traj = True)
time = np.linspace(0,5,Nl.shape[0])

entro = np.empty(traj.shape[0])
for i in range(traj.shape[0]):
    entro[i] = entropy(traj[i,:],p_stat,h)

fig,axs = plt.subplots(1,3)

axs[0].plot(v[1:-1],p_init,label="t=0")
axs[0].plot(v[1:-1],traj[100,:],label="t=1")
axs[0].plot(v[1:-1],traj[400,:],label="t=4")
axs[0].plot(v[1:-1],p_stat,"--k",label="p_stationary")
axs[0].legend()
axs[0].set_xlabel("v")
axs[0].set_ylabel("p(v,t)")

axs[1].plot(time,Nl)
axs[1].set_xlabel("t")
axs[1].set_ylabel("N(t)")

axs[2].plot(time,entro)
axs[2].set_xlabel("t")
axs[2].set_ylabel("S(t)")

plt.show()
