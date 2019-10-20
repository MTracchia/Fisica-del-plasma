
import numpy as np
import matplotlib.pylab as plt
import sys

def RK4(y, f, t, dt):
    k1 = f(t, y)
    k2 = f(t+0.5*dt, y+0.5*k1*dt)
    k3 = f(t+0.5*dt, y+0.5*k2*dt)
    k4 = f(t+1.0*dt, y+1.0*k3*dt)
    return y + (k1 + 2*k2 + 2*k3 +k4)*dt/6

F = lambda a: ( lambda x, phi: np.sqrt(2*(np.exp(phi)-1 + a*(np.sqrt(1-2*phi/a)-1))) )

a_s = np.array(sys.argv[1:], dtype=np.float64)
print(a_s)
Vo_s = [0.1, 1, 10, 50]
N = 1001
X = np.linspace(0, 10, N)
dx = X[1] - X[0]

for a in a_s:
    plt.figure()
    for Vo in Vo_s:
        PHI = np.zeros_like(X)
        PHI[0] = -Vo
        for i in range(N-1):
            PHI[i+1] = RK4(PHI[i], F(a), X[i], dx)
        plt.plot(X, PHI/Vo, "-", lw=3)
    k = np.sqrt(1.0-1.0/a)
    plt.plot(X, -np.exp(-k*X), "--")
    plt.xlabel("$x/L$")
    plt.ylabel(r"$\phi/V_o$")
    plt.legend([r"$eV_o=%.1fT_e$" %(Vo) for Vo in Vo_s], loc = "lower right")
plt.show()
