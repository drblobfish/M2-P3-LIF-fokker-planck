from LIF import *
import pandas as pd

VF = 2
Vmin = -4
VR = 1
a0 = 1
a1 = 0
b = 0.5

h = 6/384
tau = 1/20_000
T = 0.5

v0 = 0
sigma0 = 0.25

"""
# fixed tau
tau = 0.5/10_000
N_space = np.geomspace(24,768,6).astype(int)
h_space = (VF-Vmin)/N_space
norm1l = []
norminfl = []
for N in N_space :
    h = (VF-Vmin)/N
    print(f"N = {N}, h = {h}, tau = {tau}")
    print(" - h")
    ph,Nhl,vh = fokker_plank_solve(
            gaussian_initial_cond(v0,sigma0,VF,Vmin,0,N=N),
            VF,Vmin,VR,a0,a1,b,0,tau,T,N=N)
    print(" - h/2")
    ph2,Nhl,vh2 = fokker_plank_solve(
            gaussian_initial_cond(v0,sigma0,VF,Vmin,0,N=2*N),
            VF,Vmin,VR,a0,a1,b,0,tau,T,N=2*N)
    norm1 = h*np.linalg.norm(ph-ph2[1::2],ord=1)
    norm1l.append(norm1)
    norminf = np.linalg.norm(ph-ph2[1::2],ord=np.inf)
    norminfl.append(norminf)
    plt.plot(vh[1:N],ph)
    plt.plot(vh2[1:2*N],ph2)
    plt.show()

norm1l = np.array(norm1l)
norminfl = np.array(norminfl)
order1 = np.repeat(np.nan,6)
order1[:-1] = np.log2(norm1l[:-1]/norm1l[1:])
orderinf = np.repeat(np.nan,6)
orderinf[:-1] = np.log2(norminfl[:-1]/norminfl[1:])
df1 = pd.DataFrame({"h":h_space,"norm1":norm1l,"order1":order1,"norminf":norminfl,"orderinf":orderinf})

"""
print("---------------")
print("Fixed h")
# fixed h
nb_iter_space = np.geomspace(1000,32_000,6).astype(int)
tau_space = 0.5/nb_iter_space
N = 384
h = 6/N
norm1l = []
norminfl = []
order1l = []
orderinfl = []
for nb_iter in nb_iter_space :
    tau = T/nb_iter
    print(f"N = {N}, h = {h}, tau = {tau}")
    print(" - h")
    ph,Nhl,vh = fokker_plank_solve(
            gaussian_initial_cond(v0,sigma0,VF,Vmin,0,N=N),
            VF,Vmin,VR,a0,a1,b,0,tau,T,N=N)
    print(" - h/2")
    ph2,Nhl,vh2 = fokker_plank_solve(
            gaussian_initial_cond(v0,sigma0,VF,Vmin,0,N=2*N),
            VF,Vmin,VR,a0,a1,b,0,tau,T,N=2*N)
    print(" - h/4")
    ph4,Nhl,vh4 = fokker_plank_solve(
            gaussian_initial_cond(v0,sigma0,VF,Vmin,0,N=4*N),
            VF,Vmin,VR,a0,a1,b,0,tau,T,N=4*N)
    norm1 = h*np.linalg.norm(ph-ph2[1::2],ord=1)
    norm1l.append(norm1)
    norm1denom = 0.5*h*np.linalg.norm(ph2-ph4[1::2],ord=1)
    order1l.append(np.log2(norm1/norm1denom))

    norminf = np.linalg.norm(ph-ph2[1::2],ord=np.inf)
    norminfl.append(norminf)
    norminfdenom = np.linalg.norm(ph2-ph4[1::2],ord=np.inf)
    orderinfl.append(np.log2(norminf/norminfdenom))
    plt.plot(vh[1:N],ph)
    plt.plot(vh2[1:2*N],ph2)
    plt.plot(vh4[1:4*N],ph4)
    plt.show()

df2 = pd.DataFrame({"tau":h_space,"norm1":norm1l,"order1":order1l,"norminf":norminfl,"orderinf":orderinfl})
