% Isaiah Lee/
% CE 810
% HW 3
% Due Date 9/24/2025
% Problem 1
% 1. Incremental-iterative solution 2-DOF system
%           Updates/Changes
% Date                  What was Changed
% 3/9/23                Code Created and Completed
% 9/22/25               Code copied from CE810HW2
% 9/24/25               Code completed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Problem 1
EA = 5*10^7;                   %Set value for EA Newtons
z = 25;                        %z set equal to 25 mm
L = 2500;                      %l set equal to 2500 mm
Ks = 2.35;                     %spring stiffness set equal to 2.35 N/mm
Delta_w = -7;                  %Sets step value for change in W
Delta_u = 70;                  %Sets step value for change in U
Delta = [Delta_u; Delta_w];

%First step define variable using arrays as necessary
j = 1;                           %Sets location in matrix
W_knot = 0;                      %Sets start point of 0
zeros_2_1 = zeros(2,1);          %Empty 2 by 1 matrix of zeros
zeros_2_2 = zeros(2,2);          %Empty 2 by 2 matrix of zeros
k = cell(1,38);                  %Create an array with 38 elements for storing k matrices

for x = 1:38                     %Each element of k should be an empty 2 by 2 matrix
    k{1,x} = zeros_2_2;
end
delta = cell(1,38);              %Create an array with 38 elements for storing displacement matrices
for x = 1:38                     %Each element of displacement should be an empty 2 by 1 matrix
    delta{1,x} = zeros_2_1;
end
w = zeros(1, 38);             %Create an empty 1 x 38 matrix of zeroes for w
W = zeros(1, 38);             %Create an empty 1 x 38 matrix of zeroes for W
u = zeros(1, 38);             %Create an empty 1 x 38 matrix of zeroes for w
U = zeros(1, 38);             %Create an empty 1 x 38 matrix of zeroes for U
N = zeros(1, 38);             %Create an empty 1 x 38 matrix of zeroes for N
g_v = zeros(1, 38);           %Create an empty 1 x 38 matrix of zeroes for g vertical
g_h = zeros(1, 38);           %Create an empty 1 x 38 matrix of zeroes for g horizontal


for j = 1:37                         %Do 37 repetitions (259/7)

   W(j+1) = -7*j;
   U(j+1) = 70*j;
   k{1,j}= EA/L*[1 -(z+w(j))/L; -(z+w(j))/L ((z+w(j))/L)^2 + Ks*L/EA] + [0 0; 0 N(j)/L];   %Calculate step value of tangent stiffness
   delta{1,j+1} = k{1,j}^-1*Delta;                                  %Solve for changes to u, w
   u(j+1) = u(j) + delta{1,1+j}(1);                                 %Calculate step value of u
   w(j+1) = w(j) + delta{1,1+j}(2);                                 %Calculate step value of w
   N(j+1) = EA*(-u(j+1)/L+(z/L)*w(j+1)/L+1/2*(w(j+1)/L)^2);         %Calculate step value of N
   g_v(j) = N(j+1)*((z+w(j+1))/L)+Ks*w(j+1)-W(j+1);                 %Calculate step value of g_v
   g_h(j) = -N(j+1)-U(j+1);                                         %Calculate step value of g_h

   while sqrt((g_v(j))^2+(g_h(j))^2) >= 1*10^-4          %Begin iteration loops until error below desired value
       k{1,j}= EA/L*[1 -(z+w(j+1))/L; -(z+w(j+1))/L ((z+w(j+1))/L)^2 + Ks*L/EA] + [0 0; 0 N(j+1)/L];   %Calculate iteration value of tangent stiffness
       du_dw = k{1,j}^-1*[-g_h(j); -g_v(j)];                            %Calculate du and dw matrix
       du = du_dw(1);                                                   %Calculate horizontal component
       dw = du_dw(2);                                                   %Calculate vertical component
       u(j+1) = u(j+1) + du;                                            %Calculate iteration value of u
       w(j+1) = w(j+1) + dw;                                            %Calculate iteration value of w
       N(j+1) = EA*(-u(j+1)/L+(z/L)*w(j+1)/L+1/2*(w(j+1)/L)^2);         %Calculate iteration value of N
       g_v(j) = N(j+1)*((z+w(j+1))/L)+Ks*w(j+1)-W(j+1);                 %Calculate iteration value of g_v
       g_h(j) = -N(j+1)-U(j+1);                                         %Calculate iteration value of g_h
       
   end

end

W_abs = abs(W);                             %Take absolute value to plot on negative scale
w_abs = abs(w);                             %Take absolute value to plot on negative scale


figure(1)
clf
plot(w_abs,W_abs,'k')
grid on
ylim([0 300])
xlim([0 200])
xlabel('-w')
ylabel('-W_e')

figure(2)
clf
plot(u,U,'k')
grid on
ylim([0 3000])
xlim([-0.5 5.5])
xlabel('u')
ylabel('U_e')
