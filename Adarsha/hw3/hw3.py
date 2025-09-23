import numpy as np
import matplotlib.pyplot as plt

# given values
k_s = 2.35
EA =50000000
z=25
l=2500
del_W =-7
del_U=70
limit_W = -259
tolerance = 1e-4


#initial values
w = [0]
u=[0]
U =[0]
W = [0]
N = [0]
k=[]
i=0
k2=0
while W[i] > limit_W:
    w_use = w[i]
    u_use = u[i]
    n_use = N[i]
    del_w_use=del_W
    del_u_use=del_U
    force_check = np.array([[U[i]+del_u_use],
                  [W[i]+del_w_use]])
    k2=0
    while True:
        beta=(z+w_use)/l
        k_t = (EA / l) * np.array([
        [1, -beta],
        [-beta, beta**2 + (k_s*l/EA) +(n_use/EA)]])
        del_disp=np.linalg.inv(k_t) @ np.array([[del_u_use],
                  [del_w_use]])
        disp_updated = del_disp + np.array([[u_use],
                                    [w_use]])
        w_use = disp_updated[1,0]
        u_use = disp_updated[0,0]
        N_cache = EA*(-u_use/l + (z*w_use)/l**2 + 0.5*(w_use/l)**2)
        n_use = N_cache
        cal_force = np.array([
        [-n_use],
        [n_use*(z+w_use)/l+k_s*w_use]])
        g = cal_force - force_check
        del_u_use = -g[0,0]
        del_w_use = -g[1,0]
        k2=k2+1
        if ((abs(g[0,0])< tolerance and abs(g[1,0])< tolerance) ):
           break
    w.append(w_use)
    u.append(u_use)
    N.append(n_use)
    W.append(W[i] + del_W)
    U.append(U[i] + del_U)
    i += 1
    


plt.figure(figsize=(10, 6))
plt.plot([-x for x in w], [-y for y in W], 
         label="Solution -W vs -w", 
         color='red', linestyle='--', marker='o')
plt.xlabel("-w (mm)")
plt.ylabel("-W (N)")
plt.title("Problem 1: -W vs -w")
plt.legend()
plt.grid(True)
plt.tight_layout()


plt.figure(figsize=(10, 6))
plt.plot([x for x in u], [y for y in U], 
         label="Solution U vs u", 
         color='blue', linestyle='--', marker='s')
plt.xlabel("u (mm)")
plt.ylabel("U (N)")
plt.title("Problem 1: U vs u")
plt.legend()
plt.grid(True)
plt.tight_layout()


plt.show()

