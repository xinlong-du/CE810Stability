import numpy as np
import matplotlib.pyplot as plt

# Define Constants 
EA = 5e7 #N
z = 25 #mm
l = 2500 #mm
Ks = 2.35 #N/mm
TOL = 1e-4 #N 

# Initialize force, displaement, and normal force vectors
F = np.array([0, 0])
del_F = np.array([70, -7])
D = np.array([0, 0])
N = 0
w_points = [0]
u_points = [0]
W_points = [0]
U_points = [0] 

# Begin incremental solution 
for i in range(37):
    # Calculate stiffness matrix
    beta = (z + D[1])/l
    K_t = EA / l * np.array([[1, -beta],[-beta, beta**2+Ks*l/EA]]) + np.array([[0,0],[0,N/l]])
    
    # Update displacements
    del_D = np.linalg.solve(K_t, del_F)
    D = D + del_D

    # Update internal axial force
    N = EA*(-(D[0]/l) + (z/l)*(D[1]/l) + 0.5 * (D[1]/l)**2)
    
    # Increment external force vector 
    F += del_F  

    # Calculate internal force vector and check equilibrium 
    F_int = N/l * np.array([-l, z+D[1]]) + np.array([0, Ks*D[1]])
    F_unbalanced = F - F_int

    # Begin newton iteration
    while np.linalg.norm(F_unbalanced) > TOL:
        # Calculate stiffness matrix
        beta = (z + D[1])/l
        K_t = EA / l * np.array([[1, -beta],[-beta, beta**2+Ks*l/EA]]) + np.array([[0,0],[0,N/l]])
        
        # Update displacements with unbalanced force vector
        del_D = np.linalg.solve(K_t, F_unbalanced)
        D = D + del_D
 
        # Update internal axial force
        N = EA*(-(D[0]/l) + (z/l)*(D[1]/l) + 1/2 * (D[1]/l)**2)

        # Calculate internal force vector and check equilibrium 
        F_int = N/l * np.array([-l, z+D[1]]) + np.array([0, Ks*D[1]])
        F_unbalanced = F - F_int
        
    # Store current force and displacement to plotting arrays
    w_points.append(-D[1])
    u_points.append(D[0])
    W_points.append(-F[1])
    U_points.append(F[0])

# Plot the solutions 
plt.figure(figsize = (10,5))
plt.subplot(1,2,1)
plt.plot(w_points, W_points, "--",marker = "x", color = "blue")
plt.xlim(0, 200)   
plt.ylim(0, 265) 
plt.ylabel("-W (N)")
plt.xlabel("-w (mm)")
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(u_points, U_points, "--",marker = "x", color = "blue")
plt.xlim(-0.25, 5)   
plt.ylim(0, 2650) 
plt.ylabel("U (N)")
plt.xlabel("u (mm)")
plt.grid(True)
plt.show()