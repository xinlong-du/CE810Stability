#Problem 3

# In[21]:


#Problem 3
# Exact solution plot

import numpy as np
import matplotlib.pyplot as plt

# Parameters
EA = 5 * 10**7   # N
z = 25           # mm
L = 2500         # mm
Ks = 1.35        # N/mm  
w = np.linspace(0, 60, 300)  # deflection range (mm)

# Load expression
W = (EA/L**3)*(-z**2*w+1.5*z*w**2-0.5*w**3) - Ks*w

# Plot
plt.figure(figsize=(7,5))
plt.plot(w, -W, 'g-', linewidth=2)

plt.xlabel(r'Deflection $-w$ (mm)')
plt.ylabel(r'Load $-W$ (N)')
plt.title('Exact Solution: Load vs Deflection')
plt.grid(True)
plt.show()


# In[25]:


#Problem 3
#incremental solution
import numpy as np

# Parameters
EA = 5e7       # N
z = 25         # mm
L = 2500       # mm
Ks = 1.35      # N/mm
delW = -7      # N (incremental load)

# Initialize lists to store results
w_values = [0]    # displacement
N_values = [0]    # axial force
W_values = [0]    # external load
k_values = []     # tangent stiffness

# Loop until W reaches -91 N
step = 0
while W_values[-1] > -91:
    w_prev = w_values[-1]
    N_prev = N_values[-1]
    
    # Tangent stiffness
    kt = (EA / L) * ((z + w_prev) / L)**2 + (N_prev / L) + Ks
    k_values.append(kt)
    
    # Incremental displacement
    delta_w = delW / kt
    w_new = w_prev + delta_w
    
    # Updated axial force
    N_new = EA * ((z / L) * (w_new / L) + 0.5 * (w_new / L)**2)
    
    # External load
    W_new = W_values[-1]+delW
    
    # Store results
    w_values.append(w_new)
    N_values.append(N_new)
    W_values.append(W_new)
    
    # Print results
    step += 1
    print(f"Step {step}: w = {w_new:.3f} mm, N = {N_new:.2f} N, kt = {kt:.3f} N/mm, W = {W_new:.2f} N")
    
# In[25]:
# ==========================
# Exact vs Incremental Solution: Load vs Deflection
# ==========================
EA = 5 * 10**7   # N
z = 25           # mm
L = 2500         # mm
Ks = 1.35        # N/mm  
w_exact = np.linspace(0, 58, 300)  # deflection range (mm)

# Exact load expression
W_exact = (EA/L**3)*(-z**2*w_exact+1.5*z*w_exact**2-0.5*w_exact**3) - Ks*w_exact

# ==========================
# Incremental solution ( from above)
# ==========================
w_inc = np.array([0, -2.090, -4.529, -7.493, -11.337, -16.955, -27.550, -45.913,
                  -48.771, -51.057, -52.996, -54.698, -56.226, -57.617])
W_inc = np.array([0,-7, -14, -21, -28, -35, -42, -49,
                  -56, -63, -70, -77, -84, -91])

# Reverse signs to match exact solution
w_inc_pos = -w_inc
W_inc_pos = -W_inc

# ==========================
# Plot
# ==========================
plt.figure(figsize=(8,5))

# Exact solution line
plt.plot(w_exact, -W_exact, 'g-', linewidth=3, label='Exact Solution')

# Incremental solution points and connecting line
plt.plot(w_inc_pos, W_inc_pos, 'r.-', markersize=10, linewidth=3, label='Incremental Solution')

plt.xlabel(r'Deflection $-w$ (mm)')
plt.ylabel(r'Load $-W$ (N)')
plt.title('Exact vs Incremental Solution: Load vs Deflection')
plt.grid(True)
plt.legend()
plt.show()
