#Problem 4
#incremental iterative solution
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Exact solution
# -------------------------
EA = 5e7
z = 25
L = 2500
Ks = 1.35
w_exact = np.linspace(0, 60, 300)
W_exact = (EA/L**3)*(-z**2*w_exact + 1.5*z*w_exact**2 - 0.5*w_exact**3) - Ks*w_exact

# -------------------------
# Incremental-iterative solution
# -------------------------
delW = -7
W_final = -91
tol = 1e-4

def axial_force(w):
    return EA * ((z/L) * (w/L) + 0.5*(w/L)**2)

def tangent_stiffness(w, N):
    return (EA/L) * ((z+w)/L)**2 + (N/L) + Ks

def equilibrium_force(w, N):
    return N * ((z+w)/L) + Ks*w

w = 0.0
N = 0.0
W_applied = 0.0

W_list = [W_applied]
w_list = [w]

while W_applied > W_final:
    W_applied += delW
    # Predictor
    kt = tangent_stiffness(w, N)
    dw = delW / kt
    w += dw
    N = axial_force(w)
    
    # Corrector iterations
    W_eq = equilibrium_force(w, N)
    g = W_eq - W_applied
    while abs(g) > tol:
        kt = tangent_stiffness(w, N)
        dw = -g / kt
        w += dw
        N = axial_force(w)
        W_eq = equilibrium_force(w, N)
        g = W_eq - W_applied

    W_list.append(W_applied)
    w_list.append(w)

w_arr = np.array(w_list)
W_arr = np.array(W_list)

# -------------------------
# Plot comparison
# -------------------------
plt.figure(figsize=(8,5))

# Exact solution line
plt.plot(w_exact, -W_exact, 'g-', linewidth=3, label='Exact Solution')

# Incremental-iterative solution
plt.plot(-w_arr, -W_arr, 'r.-', markersize=10, linewidth=3, label='Incremental-Iterative Solution')

plt.xlabel('Displacement -w (mm)')
plt.ylabel('Applied Load -W (N)')
plt.title('Exact vs Incremental-Iterative Load–Displacement Curve')
plt.grid(True)
plt.legend()
plt.show()
