import numpy as np
import matplotlib.pyplot as plt


# PROBLEM 1.
# range of x
x1 = np.linspace(0, 2.5, 400)

# function: where y1 = -(W/EA)*(l/z^3) and x1 = -w/z
y1 = x1 - 1.5 * x1**2 + 0.5 * x1**3

# plot
plt.figure(figsize=(8, 6))
plt.plot(x1, y1, color='teal')
plt.axhline(0, color='gray', lw=0.5)
plt.axvline(0, color='gray', lw=0.5)
plt.xlabel(r"$-\frac{w}{z}$")
plt.ylabel(r"$-\frac{W}{EA} \cdot \frac{l}{z^3}$")
plt.ylim(-0.3,0.5)
plt.xticks(np.linspace(0,2.0,3))
plt.yticks(np.linspace(-0.2,0.2,3))
plt.title("Problem 1. Load-Deflection Relationship of Figure (a)")
plt.grid(True)
plt.tight_layout()
plt.show()




# PROBLEM 2.
# range of x
x2 = np.linspace(0, 2.5, 400)

# function: same as Problem 1, with k_sw added. k_s = (EAz^2)/(2l^3)
# after simplification:
y2 = 1.5 * x2 - 1.5 * x2**2 + 0.5 * x2**3

# plot: the local max of the force-displacement curve in Problem 1 is avoided
plt.figure(figsize=(8, 6))
plt.plot(x2, y2, color='darkblue')
plt.axhline(0, color='gray', lw=0.5)
plt.axvline(0, color='gray', lw=0.5)
plt.xlabel(r"$-\frac{w}{z}$")
plt.ylabel(r"$-\frac{W}{EA} \cdot \frac{l}{z^3}$")
plt.ylim(0,2.2)
plt.xticks(np.linspace(0,2.0,3))
plt.yticks(np.linspace(0,2,3))
plt.title("Problem 2. Load-Deflection Relationship of Figure (b)")
plt.grid(True)
plt.tight_layout()
plt.show()




# PROBLEM 3.
# constants (given):
EA = 5e7          # N
z = 25            # mm
l = 2500          # mm
k_s = 1.35        # N/mm
deltaW = -7       # (increment)

# assume the initial conditions:
w = 0.0           # mm
N = 0.0           # N
W_actual = 0.0    # N
W_target = 0.0    # N

# for plotting, keep a list of values:
w_list = []
W_list = []

# for output: use a table format to keep track of the individual values at each increment
print(f"{'Step':<5} {'W_target (N)':<15} {'W_actual (N)':<15} {'w (mm)':<10} {'N (N)':<12} {'k_t (N/mm)':<15}")

# to count the number of increments:
step = 0

# use a while loop to perform the incremental solution until W_target = -91 N
while W_target > -91 :
    step += 1

    # apply the compressive load
    W_target += deltaW

    # calculate the tangent stiffness, k_t, using Eq. (1.5.7) from class notes
    k_t = ((EA / l) * ((z + w) / l) ** 2) + (N / l) + k_s

    # calculate the displacement at this point on the curve
    delta_w = deltaW / k_t

    # use the previously calculated value in the following calculations for N and the next k_t.
    w += delta_w

    # calculate the internal force, N, using Eq. (1.5.3) from class notes
    N = EA * (((z / l) * (w / l)) + 0.5 * (w / l) ** 2)

    # using Eq. (1.5.6) from class notes, calculate the actual load at this increment
    W_actual = ((N * (z + w)) / l) + (k_s * w)

    # for plotting, add these results to the list of values
    w_list.append(w)
    W_list.append(W_target)

    # print results in table
    print(f"{step:<5} {W_target:<15.2f} {W_actual:<15.3f} {w:<10.5f} {N:<12.2f} {k_t:<15.2f}")

# for plotting the exact solution on the same figure, use the following (W_exact equation is from class notes):
w_exact = np.linspace(min(w_list), 0, 300)
W_exact = ((EA / l ** 3) * ((z ** 2) * w_exact + 1.5 * z * w_exact ** 2 + 0.5 * w_exact ** 3)) + k_s * w_exact

# plot
plt.figure(figsize=(8, 6))
plt.plot(w_list, W_list, marker='o', linestyle='-', color='red')
plt.title('Problem 3. Incremental vs. Exact Solution of Figure (b)')
plt.xlabel('-w (mm)')
plt.ylabel('-W (N)')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xticks(np.linspace(-60, 0, 7))
plt.yticks(np.linspace(-90, 0, 10))
plt.plot(w_exact, W_exact, color='darkblue')
plt.legend(['Incremental', 'Exact'])
plt.grid(True)
plt.tight_layout()
plt.show()




# PROBLEM 4
# assume the initial conditions:
w = 0.0           # mm
N = 0.0           # N
W_target = 0.0    # N

# tolerance for convergence test of the unbalanced force
tolerance = 1e-4

# for plotting, keep a list of values:
w_list = []
W_list = []

# to count the number of load steps:
step = 0

# use a while loop to perform the incremental-iterative solution until W_target = -91 N
while W_target > -91:
    step += 1
    W_target += deltaW

    # use the incremental solution as the predictor
    w_iter = w

    # to ensure the loop begins, assume a large enough residual:
    residual = 1e6
    iteration = 0
    first_residual = None

    # Newton-Raphson iterations
    while abs(residual) > tolerance:
        iteration += 1

        # calculate the internal force N from Eq. (1.5.3)
        N = EA * ((z / l) * (w_iter / l) + 0.5 * (w_iter / l) ** 2)

        # calculate the load W_actual from Eq. (1.5.6)
        W_actual = (N * (z + w_iter)) / l + k_s * w_iter

        # check equilibrium (calculate unbalanced force)
        residual = W_actual - W_target

        # calculate tangent stiffness k_t using Eq. (1.5.7)
        k_t = (EA / l) * ((z + w_iter) / l) ** 2 + N / l + k_s

        # Newton-Raphson update for displacement
        delta_w = -residual / k_t
        w_iter += delta_w

        # Update converged value
        w = w_iter

        print(
            f"Iter {iteration}: w_iter={w_iter:.6f}, N={N:.2f}, W_actual={W_actual:.4f}, residual={residual:.4f}, k_t={k_t:.4f}, delta_w={delta_w:.6f}")

    # save values for plotting
    w_list.append(w)
    W_list.append(W_actual)

# plot
plt.figure(figsize=(8, 6))
plt.plot(w_list, W_list, 'o-', label='Newton-Raphson', color='darkorange')

# plot exact solution from class for comparison
w_exact = np.linspace(min(w_list), 0, 300)
W_exact = ((EA / l ** 3) * ((z ** 2) * w_exact + 1.5 * z * w_exact ** 2 + 0.5 * w_exact ** 3)) + k_s * w_exact
plt.plot(w_exact, W_exact, color='darkblue')

plt.title('Problem 4. Incremental-Iterative Solution vs. Exact Solution of Figure (b)')
plt.xlabel('-w (mm)')
plt.ylabel('-W (N)')
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.xticks(np.linspace(-60, 0, 7))
plt.yticks(np.linspace(-90, 0, 10))
plt.grid(True)
plt.legend(['Incremental-Iterative', 'Exact'])
plt.tight_layout()
plt.show()