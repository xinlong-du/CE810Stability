import numpy as np
import matplotlib.pyplot as plt


# PROBLEM 1
# constants (given)
EA = 5e7          # N
z = 25            # mm
l = 2500          # mm
k_s = 2.35        # N/mm
deltaW = -7       # N (increment)
deltaU = 70       # N (increment)

# assume the initial displacements:
w = 0.0           # mm
u = 0.0           # mm
N = 0.0           # N
W_e = 0.0    # N
U_e = 0.0    # N

# tolerance for convergence test of the unbalanced force
tolerance = 1e-4

# for plotting, keep a list of values
u_list, U_list = [], []
w_list, W_list = [], []

# to count the number of load steps:
step = 0

# use a while loop to perform the incremental-iterative solution until W_target = -259 N
while W_e > -259:
    step += 1
    U_e += deltaU
    W_e += deltaW

    # predictor: use last converged displacement
    w_iter = w
    u_iter = u

    # to ensure the loop begins, assume a large enough residual:
    residual_norm = 1e6
    iteration = 0

    # Newton-Raphson iterations
    while abs(residual_norm) > tolerance:
        iteration += 1

        # changes in geometry:
        a = l - u_iter
        b = z - w_iter
        L = np.sqrt(a**2 + b**2)

        # axial force in bar:
        N = EA * ((-u_iter / L) + (z / L) * (w_iter / L) + 0.5 * (w_iter / L) ** 2)

        # q_i = internal force vector = [qi1, qi2]
        qi1 = -N
        qi2 = N * (z + w_iter) / L + k_s * w_iter
        q_i = np.array([qi1, qi2])

        # q_e = external force vector
        q_e = np.array([U_e, W_e])

        # unbalanced force vector = residual = g = q_i - q_e
        residual = q_i - q_e
        residual_norm = np.linalg.norm(residual)

        # tangent stiffness (∂q_i / ∂p)
        # derivatives of L with respect to u and w:
        dL_du = -a / L
        dL_dw = -b / L
        # derivatives of N with respect to u and w:
        dN_du = (EA / l) * dL_du
        dN_dw = (EA / l) * dL_dw
        # qi1 = -N
        dqi1_du = -dN_du
        dqi1_dw = -dN_dw
        # qi2 = N*(z+w)/L +k_s*w
        dqi2_du = dN_du * (z + w_iter) / L + N * (0 * L - (z + w_iter) * dL_du) / (L**2)
        dqi2_dw = dN_dw * (z + w_iter) / L + N * (1 * L - (z + w_iter) * dL_dw) / (L**2) + k_s

        Kt = np.array([[dqi1_du, dqi1_dw],
                       [dqi2_du, dqi2_dw]])

        # Newton-Raphson update for displacement
        delta_p = np.linalg.solve(Kt, -residual)
        u_iter += delta_p[0]
        w_iter += delta_p[1]

        print(
            f"Step {step}, Iteration {iteration}: "
            f"u={u_iter:.6f}, w={w_iter:.6f}, "
            f"N={N:.2f}, residual={residual_norm:.4e}"
        )

        # update converged value
        u, w = u_iter, w_iter
        u_list.append(u)
        w_list.append(w)
        U_list.append(U_e)
        W_list.append(W_e)

# plot
# Ue vs. u
plt.figure(figsize=(8, 6))
plt.plot(u_list, U_list, 'o-', color='red')
plt.title('Ue vs u')
plt.xlabel('u (mm)')
plt.ylabel('Ue (N)')
plt.grid(True)

# -We vs. -w
plt.figure(figsize=(8, 6))
plt.plot(w_list, W_list, 'o-', color='darkblue')
plt.title('-We vs -w')
plt.xlabel('-w (mm)')
plt.ylabel('-We (N)')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.grid(True)

plt.show()
