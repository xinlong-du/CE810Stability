import numpy as np
import matplotlib.pyplot as plt

# given values
k_s = 1.35
EA =50000000
z=25
l=2500
del_W =-7
limit_W = -91
tolerance = 1e-4


#initial values
w = [0]
W = [0]
N = [0]
k=[]
i=0
k2=0
while W[i] > limit_W:
    w_use = w[i]
    n_use = N[i]
    del_w_use=del_W
    while True:
        k_cache = (EA/l**3)*(z+w_use)**2 + n_use/l + k_s
        w_cache = w_use + del_w_use/ k_cache
        N_cache = EA* ((z * w_cache) / l**2 + 0.5 * (w_cache / l)**2)
        W_cache = N_cache *(z+w_cache)/l + k_s*w_cache
        g = W_cache -(W[i]+del_W)
        w_use = w_cache
        n_use = N_cache
        del_w_use = -g
        
        k2=k2+1
        if (abs(g)< tolerance):
            break
    k_value = k_cache
    k.append(k_value)
    w.append(w_cache)
    N.append(N_cache)
    W.append(W[i] + del_W)
    i += 1

# Exact solution
w_act = np.linspace(0, -58, 500)
W_act = (EA/(l/z)**3) * ((w_act/z) + 1.5*(w_act/z)**2 + 0.5*(w_act/z)**3) + k_s*w_act


plt.figure(figsize=(10, 6))
plt.plot([-x for x in w_act], [-y for y in W_act], label="Exact Solution", color='blue', linewidth=2)
plt.plot([-x for x in w], [-y for y in W], label="Incremental Solution", color='red', linestyle='--', marker='o')

plt.xlabel("-w ")
plt.ylabel("-W ")
plt.title("Problem 4: Comparison of Exact and Iterative-Incremental Solutions ")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

