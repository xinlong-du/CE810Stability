# CE 810 Problem 1_ HW3
# Pragun_Shrestha

import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Problem / Control Parameters
# -----------------------------
AxialStiff = 50_000_000.0   # [N]
OffsetZ    = 25.0           # [mm]
LengthL    = 2500.0         # [mm]
SpringK    = 2.35           # [N/mm]

Tolerance  = 1e-4           # convergence tolerance
MaxNewton  = 100            # max Newton iterations per load step

LoadStepU  = 70.0           # incremental axial load [N]
LoadStepW  = -7.0           # incremental transverse load [N]
TargetWe   = -259.0         # transverse load limit to stop incrementing

# -----------------------------
# Helper Functions
# -----------------------------
def beta_val(w_def):
    return (OffsetZ + w_def) / LengthL

def compute_strain(u_def, w_def):
    return -u_def / LengthL + (OffsetZ / LengthL) * (w_def / LengthL) + 0.5 * (w_def / LengthL)**2

def axial_force(u_def, w_def):
    return AxialStiff * compute_strain(u_def, w_def)

def tangent_matrix(w_def, N_val):
    b = beta_val(w_def)
    k_linear = (AxialStiff / LengthL) * np.array([
        [1.0,        -b],
        [-b,   b**2 + (SpringK * LengthL / AxialStiff)]
    ])
    k_geo = np.array([[0.0, 0.0],
                      [0.0, N_val / LengthL]])
    return k_linear + k_geo

def residual_vec(u_def, w_def, U_load, W_load):
    N_val = axial_force(u_def, w_def)
    r1 = -U_load - N_val
    r2 = -W_load + N_val * ((OffsetZ + w_def)/LengthL) + SpringK * w_def
    return np.array([r1, r2], dtype=float)

# -----------------------------
# Incremental-Iterative Solver
# -----------------------------
def run_solver(u_start=0.0, w_start=0.0, U_start=0.0, W_start=0.0, verbose=True):
    u_curr, w_curr = u_start, w_start
    U_curr, W_curr = U_start, W_start

    step_history = {"step": [], "u": [], "w": [], "N": [], "U": [], "W": [], "iters": []}
    step_count = 0

    while W_curr > TargetWe:
        step_count += 1
        # increment loads
        U_curr += LoadStepU
        W_curr += LoadStepW
        load_increment = np.array([LoadStepU, LoadStepW], dtype=float)

        # predictor (first guess at zero deflection & zero axial force)
        kt_pred = tangent_matrix(w_def=0.0, N_val=0.0)
        dq_pred = np.linalg.solve(kt_pred, load_increment)
        u_iter = u_curr + dq_pred[0]
        w_iter = w_curr + dq_pred[1]

        # Newton-Raphson corrector
        for n_iter in range(MaxNewton):
            res = residual_vec(u_iter, w_iter, U_curr, W_curr)
            if np.linalg.norm(res) <= Tolerance:
                break
            N_iter = axial_force(u_iter, w_iter)
            kt = tangent_matrix(w_iter, N_iter)
            delta_q = -np.linalg.solve(kt, res)
            u_iter += delta_q[0]
            w_iter += delta_q[1]

        # update state
        u_curr, w_curr = u_iter, w_iter
        N_val = axial_force(u_curr, w_curr)

        # store history
        step_history["step"].append(step_count)
        step_history["u"].append(u_curr)
        step_history["w"].append(w_curr)
        step_history["N"].append(N_val)
        step_history["U"].append(U_curr)
        step_history["W"].append(W_curr)
        step_history["iters"].append(n_iter + 1)

        if verbose:
            print(f"Step {step_count:3d}: U={U_curr:.3f}, W={W_curr:.3f}, u={u_curr:.6f}, w={w_curr:.6f}, N={N_val:.6f}")

    return step_history

# -----------------------------
# Plotting
# -----------------------------
def plot_history(hist_data):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), dpi=130)
    # -W vs -w in red
    ax[0].plot(-np.array(hist_data["w"]), -np.array(hist_data["W"]), marker='o', color='red')
    ax[0].set_xlabel("-w [mm]"); ax[0].set_ylabel("-W [N]"); ax[0].grid(ls=':')
    # U vs u in red
    ax[1].plot(hist_data["u"], hist_data["U"], marker='s', color='red')
    ax[1].set_xlabel("u [mm]"); ax[1].set_ylabel("U [N]"); ax[1].grid(ls=':')
    fig.tight_layout()
    plt.show()

# -----------------------------
# Run example
# -----------------------------
if __name__ == "__main__":
    hist = run_solver(verbose=True)
    plot_history(hist)
