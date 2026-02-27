import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

VF = 2
Vmin = -4
VR = 1
a0 = 1
a1 = 0
b = 3

N = 10
tau = 0.01
h = (VF-Vmin)/N

v = np.linspace(Vmin,VF,N+1)
p = np.zeros(N+1)
p[N//2] = 1

Nh = a0 * p[N-1] / (h - a1 * p[N-1])
aNh = a0 + a1 * Nh
M = np.exp(- (v - b * Nh)**2 / (2*(a0 + a1 * Nh)))
Mmidpoints = 1/(0.5 * (1/M[:-1] + 1/M[1:]))

diagonal = np.zeros(N+1)
diagonal[2:N-1] = 1 + (tau/h**2) * aNh * (Mmidpoints[1:N-2] + Mmidpoints[2:N-1])/M[2:N-1]
diagonal[1]

sparray = sp.diags_array([np.ones(5),np.ones(4)],offsets=[0,1])
