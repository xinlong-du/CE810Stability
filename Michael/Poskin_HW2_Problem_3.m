% Poskin CE 810 HW2 Problem 3

clear; clc; close all

%------------INITIAILIZE VARIABLES----------------
W_i = 0; W_f = -91; delta_W = -7; % Force range and step size
W = W_i:delta_W:W_f; n = size(W,2);

EA = 5*10^7; z = 25; L = 2500; k_s = 1.35; % Axial stiffness, initial height, and string stiffness

w = zeros(1,n); N = 0; % Displacement and constitutive load are initially 0

%--------INCREMENTAL ITERATIVE SOLUTION-----------
% For each load step calculate tangential stiffness (k_t), used to 
% approximate displacement.
for i = 2:n
    w(i) = w(i-1);
    k_t = EA/L*((z + w(i))/L)^2 + N/L + k_s;
    delta_w = k_t^-1*delta_W;
    w(i) = w(i) + delta_w;
    N = EA*((z*w(i) + 0.5*w(i)^2)/L^2);
end

%-------------------PLOTTING----------------------
% Plots -w v. -W, numerical solution
figure(1)
plot(-w,-W)
hold on

% Plots -w v. -W, analytical solution
w_anl = linspace(0, w(n),1000);
W_anl = EA/L^3*(z^2*w_anl + 3/2*z*w_anl.^2 + 1/2*w_anl.^3) + k_s*w_anl; % Analytical solution
plot(-w_anl,-W_anl)

title('Poskin HW 2 Problem 3: Nonlinear Buckling')
xlabel('-w [mm]')
ylabel('-W [N]')
legend('incremental solution', 'analytical solution')