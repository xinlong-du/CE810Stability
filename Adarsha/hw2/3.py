import numpy as np
import matplotlib.pyplot as plt

# given values
k_s = 1.35
EA =50000000
z=25
l=2500
del_W =-7
limit_W = -91



#initial values
w = [0]
W = [0]
N = [0]
k=[]
i=0

while W[i] > limit_W:
    k_value = (EA/l**3)*(z+w[i])**2 + N[i]/l + k_s
    k.append(k_value)
    w.append(w[i] + del_W / k[i])
    N.append(EA * ((z * w[i+1]) / l**2 + 0.5 * (w[i+1] / l)**2))
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
plt.title("Problem 3: Comparison of Exact and Incremental Solutions ")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()