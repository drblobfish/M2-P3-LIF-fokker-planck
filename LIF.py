import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

VF = 2
Vmin = 0
VR = 1
a0 = 1
a1 = 0
b = 1.5

N = 100
tau = 0.001


h = (VF-Vmin)/N
VR_index = int(np.round((VR-Vmin)*N/(VF-Vmin)))
v = np.linspace(Vmin,VF,N+1)

# p = np.zeros(N-1)
# p[N//2] = 1 # intial distribution

v0 = 1.5
sigma0 = 0.005
p = np.exp(-(v[1:N]-v0)**2/sigma0)
p /= np.sum(p)

# at each step
plt.plot(v[1:N],p)
plt.show()
for i in range(100):
    Nh = a0 * p[N-2] / (h - a1 * p[N-2])
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
    print(np.sum(p))
    if i%10 == 0 :
        plt.plot(v[1:N],p)
        plt.show()


