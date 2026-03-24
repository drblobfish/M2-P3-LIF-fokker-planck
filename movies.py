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
T = 3.38
v0 = -1
sigma0 = 0.5

p = gaussian_initial_cond(v0,sigma0,VF,Vmin,h)

fig, ax = plt.subplots()
nb_step = 100
ud = UpdateDist(fig,ax,p,VF,Vmin,VR,a0,a1,b,h,tau,T,nb_step)
ani = FuncAnimation(fig, ud,init_func=ud.start, frames=np.linspace(0, T, nb_step), blit=True,repeat= False)

ani.save("movies/movie_blowup_top.mp4")


tau = 1e-4
# fig 2 bottom
VF = 2
Vmin = 0
VR = 1
a0 = 1
a1 = 0
b = 1.5
T = 0.043
v0 = 1.5
sigma0 = 0.005
p = gaussian_initial_cond(v0,sigma0,VF,Vmin,h)

fig, ax = plt.subplots()
nb_step = 100
ud = UpdateDist(fig,ax,p,VF,Vmin,VR,a0,a1,b,h,tau,T,nb_step)
ani = FuncAnimation(fig, ud,init_func=ud.start, frames=np.linspace(0, T, nb_step), blit=True,repeat= False)

ani.save("movies/movie_blowup_bottom.mp4")
