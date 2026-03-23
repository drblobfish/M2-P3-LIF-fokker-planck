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

N_stable = 0.1924
N_instable = 2.319

p_stable,v = stationnary_initial_cond(N_stable,a0,b,VF,Vmin,VR,h)
p_instable,_ = stationnary_initial_cond(N_instable,a0,b,VF,Vmin,VR,h)

p_instable1,_,_ = fokker_plank_solve(p_instable,VF,Vmin,VR,a0,a1,b,h,tau,0.2)
plt.plot(v[1:-1],p_instable)
plt.show()

p_instable1,_,_ = fokker_plank_solve(p_instable,VF,Vmin,VR,a0,a1,b,h,tau,0.5)
p_instable2,_,_ = fokker_plank_solve(p_instable,VF,Vmin,VR,a0,a1,b,h,tau,1)
p_instable3,_,_ = fokker_plank_solve(p_instable,VF,Vmin,VR,a0,a1,b,h,tau,4)
p_instable4,_,_ = fokker_plank_solve(p_instable,VF,Vmin,VR,a0,a1,b,h,tau,4.5)
_,Nl_instable,_ = fokker_plank_solve(p_instable,VF,Vmin,VR,a0,a1,b,h,tau,5)
time = np.linspace(0,5,Nl_instable.shape[0])
_,Nl_stable,_ = fokker_plank_solve(p_stable,VF,Vmin,VR,a0,a1,b,h,tau,5)
p_stable1,_,_ = fokker_plank_solve(p_stable,VF,Vmin,VR,a0,a1,b,h,tau,0.005)
p_stable2,_,_ = fokker_plank_solve(p_stable,VF,Vmin,VR,a0,a1,b,h,tau,0.05)
p_stable3,_,_ = fokker_plank_solve(p_stable,VF,Vmin,VR,a0,a1,b,h,tau,0.5)

p_instable1_full = np.zeros_like(v)
p_instable2_full = np.zeros_like(v)
p_instable3_full = np.zeros_like(v)
p_instable4_full = np.zeros_like(v)
p_instable1_full[1:-1] = p_instable1
p_instable2_full[1:-1] = p_instable2
p_instable3_full[1:-1] = p_instable3
p_instable4_full[1:-1] = p_instable4

p_stable1_full = np.zeros_like(v)
p_stable2_full = np.zeros_like(v)
p_stable3_full = np.zeros_like(v)
p_stable1_full[1:-1] = p_stable1
p_stable2_full[1:-1] = p_stable2
p_stable3_full[1:-1] = p_stable3

fig,axs = plt.subplots(2,2)
axs[0,0].plot(time,Nl_instable)
axs[0,0].set_xlabel("t")
axs[0,0].set_ylabel("N(t)")
axs[0,1].plot(time,Nl_stable)
axs[0,1].set_xlabel("t")
axs[0,1].set_ylabel("N(t)")
axs[1,0].plot(v,p_instable1_full,label = "t=0.5")
axs[1,0].plot(v,p_instable2_full,label = "t=1")
axs[1,0].plot(v,p_instable3_full,label = "t=4")
axs[1,0].plot(v,p_instable4_full,label = "t=4.5")
axs[1,0].legend()
axs[1,0].set_xlabel("v")
axs[1,0].set_ylabel("p(v,t)")
axs[1,1].plot(v,p_stable1_full,label = "t=0.005")
axs[1,1].plot(v,p_stable2_full,label = "t=0.05")
axs[1,1].plot(v,p_stable3_full,label = "t=0.5")
axs[1,1].legend()
axs[1,1].set_xlabel("v")
axs[1,1].set_ylabel("p(v,t)")

plt.show()

# fig.savefig("figs/fig3.pdf")
# fig.savefig("figs/fig3.png")

