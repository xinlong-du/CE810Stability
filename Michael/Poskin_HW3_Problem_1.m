% Poskin CE 810 HW3 Problem 1

clear; clc; close all

%------------INITIAILIZE VARIABLES----------------
W_i = 0; W_f = -259; delta_W = -7; % Force W range and step size
W = W_i:delta_W:W_f; n = size(W,2);

U_i = 0; delta_U = 70; U_f = U_i + delta_U*(n-1); % Force U range and step size
U = U_i:delta_U:U_f;

EA = 5*10^7; z = 25; L = 2500; k_s = 2.35; % Axial stiffness, initial height, length and string stiffness

w = zeros(1,n); u = zeros(1,n); N = 0; % Displacement and constitutive load are initially 0

%--------INCREMENTAL ITERATIVE SOLUTION-----------
% For each load step calculate tangential stiffness (k_t), used to 
% approximate displacement. Iteration with Newton-Rapson method is used to remove residual 
% force (g) from equilibrium. 
for i = 2:n
    g = [delta_U; delta_W];
    u(i) = u(i-1);
    w(i) = w(i-1);
    while norm(g) > 10^-4
        beta = (z+w(i))/L;
        k_t = EA/L*[1, -beta; -beta, beta^2 + k_s*L/EA] + [0, 0; 0, N/L];

        delta_p = k_t^-1*g;
        u(i) = u(i) + delta_p(1);
        w(i) = w(i) + delta_p(2);

        N = EA*(-u(i)/L + (z*w(i) + 0.5*w(i)^2)/L^2);
        g = [U(i); W(i)] - N/L*[-L; z+w(i)] - [0; k_s*w(i)];
    end
end

%-------------------PLOTTING----------------------
% Plot -w v. -W, numerical solution
figure(1)
plot(-w,-W)
title('Poskin HW 3 Problem 1: Nonlinear Buckling 2-DOF, -w v. -W')
xlabel('-w [mm]')
ylabel('-W [N]')

% Plot u v. U, numerical solution
figure(2)
plot(u,U)
title('Poskin HW 3 Problem 1: Nonlinear Buckling 2-DOF, u v. U')
xlabel('u [mm]')
ylabel('U [N]')
