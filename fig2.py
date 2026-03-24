from LIF import *

h = 0.02
tau = 1e-3

# fig 2 top
VF = 2
Vmin = -4
VR = 1
a0 = 1
a1 = 0
b = 3


v0 = -1
sigma0 = 0.5

p1,_,v = fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,2.95)

p2,_,_= fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,3.15)

p3,_,_= fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,3.35)

p4,Nl,_= fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,3.38)
t = np.linspace(0,3.38,Nl.shape[0])

p_full1 = np.zeros_like(v)
p_full1[1:-1] = p1
p_full2 = np.zeros_like(v)
p_full2[1:-1] = p2
p_full3 = np.zeros_like(v)
p_full3[1:-1] = p3
fig,axs = plt.subplots(1,2,layout="constrained",figsize=(10,5))
axs[1].plot(v,p_full1,label = "t=2.95")
axs[1].plot(v,p_full2,label = "t=3.15")
axs[1].plot(v,p_full3,label = "t=3.35")
axs[1].legend()
axs[1].set_xlabel("v")
axs[1].set_ylabel("p(v,t)")
axs[0].plot(t,Nl)
axs[0].set_xlabel("t")
axs[0].set_ylabel("N(t)")
fig.savefig("figs/fig1_top.pdf")
fig.savefig("figs/fig1_top.png")

# fig 2 bottom
VF = 2
Vmin = 0
VR = 1
a0 = 1
a1 = 0
b = 1.5

v0 = 1.5
sigma0 = 0.005

p1,_,v = fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,0.0325)

p2,_,_= fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,0.0365)

p3,_,_= fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,0.0405)

p4,Nl,_= fokker_plank_solve(
    gaussian_initial_cond(v0,sigma0,VF,Vmin,h),
    VF,Vmin,VR,a0,a1,b,h,tau,0.043)
t = np.linspace(0,0.043,Nl.shape[0])

p_full1 = np.zeros_like(v)
p_full1[1:-1] = p1
p_full2 = np.zeros_like(v)
p_full2[1:-1] = p2
p_full3 = np.zeros_like(v)
p_full3[1:-1] = p3
fig,axs = plt.subplots(1,2,layout="constrained",figsize=(10,5))
axs[1].plot(v,p_full1,label = "t=0.0325")
axs[1].plot(v,p_full2,label = "t=0.0365")
axs[1].plot(v,p_full3,label = "t=0.0405")
axs[1].legend()
axs[1].set_xlabel("v")
axs[1].set_ylabel("p(v,t)")
axs[0].plot(t,Nl)
axs[0].set_xlabel("t")
axs[0].set_ylabel("N(t)")
fig.savefig("figs/fig1_bottom.pdf")
fig.savefig("figs/fig1_bottom.png")
