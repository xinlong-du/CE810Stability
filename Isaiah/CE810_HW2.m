% Isaiah Lee/
% CE 810
% HW 2
% Due Date 9/10/2025
% Problem 1-4
% 1. Plot Normalized Load-Deflection for problem a
% 2. Plot Normalized Load-Deflection for problem b
% 3. Incremental solution for problem b
% 4. Incremental-iterative solution for problem b
%           Updates/Changes
% Date                  What was Changed
% 3/9/23                Code Created and Completed
% 9/10/25               Code copied from CE810HW1, Code changed
%                       and completed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem 1
j = 1;                          %Sets location in matrix
w = 0;                          %Initial normalized value for -w/z

while w >= -3                                %Repeat until w/z = 3

   W(j) = -w - 3/2*w^2-1/2*w^3;             %Calculate normalized W,-W/EA*(l/z)^3 
   
   j = j+1;                                 %Increment w and repeat
   w = w-0.001;
   
end

x = linspace(0,3,3001);

figure(1)
clf
plot(x,W,'k')
grid on
ylim([-0.5 0.5])
xlim([0 3])
xlabel('-w/z')
ylabel('-W/EA*(l/z)^3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Problem 2
j = 1;                          %Sets location in matrix
w = 0;                          %Initial normalized value for -w/z

while w >= -3                                %Repeat until w/z = 3

   W(j) = -3/2*w - 3/2*w^2-1/2*w^3;             %Calculate normalized W,-W/EA*(l/z)^3 
   
   j = j+1;                                 %Increment w and repeat
   w = w-0.001;
   
end

x = linspace(0,3,3001);

figure(2)
clf
plot(x,W,'k')
grid on
ylim([0 2])
xlim([0 3])
xlabel('-w/z')
ylabel('-W/EA*(l/z)^3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Problem 3
EA = 5*10^7;                   %Set value for EA Newtons
z = 25;                        %z set equal to 25 mm
L = 2500;                      %l set equal to 2500 mm
Ks = 1.35;                     %spring stiffness set equal to 1.35 N/mm
delta_w = -7;                  %Sets step value for change in W

%First step define variables
j = 1;                           %Sets location in matrix
k = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for k
w = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for w
W = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for W
N = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for N

%Incremental solution
for j = 1:13                         %Do 13 repetitions (91/7)

   k(j) = EA/L*(z/L)^2+EA/L^3*(2*z*w(j)+(w(j))^2)+N(j)/L+Ks;        %Calculate step value of tangent stiffness
   w(j+1) = w(j) + delta_w/k(j);                                    %Calculate step value of w
   N(j+1) = EA*((z/L)*w(j+1)/L+1/2*(w(j+1)/L)^2);                   %Calculate step value of N
   W(j+1) = -7*j;
   
   j = j+1;                                 %Increment w and repeat
   
end

W_abs = abs(W);                             %Take absolute value to plot on negative scale
w_abs = abs(w);                             %Take absolute value to plot on negative scale

%Exact solution
w_ex = 0;                           %Starting w value of exact solution
j=1;

while w_ex >= -57.5                         %Repeat until w reaches -57.5

    W_ex(j) = EA/L^3*(z^2*w_ex+3/2*z*w_ex^2+1/2*w_ex^3)+Ks*w_ex;     %Solve for W

    j = j+1;                                 %Increment w and repeat
    w_ex = w_ex-0.1;

end

x = linspace(0,57.5,575);
W_ex_abs = abs(W_ex);                         %Take absolute value to plot on negative scale

figure(3)
clf
plot(w_abs,W_abs,'k')
grid on
ylim([0 120])
xlim([0 60])
xlabel('-w')
ylabel('-W')

hold on
plot(x,W_ex_abs,'b--')
grid on
ylim([0 120])
xlim([0 60])
xlabel('-w')
ylabel('-W')
legend ('Incremental', 'Exact')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Problem 4
EA = 5*10^7;                   %Set value for EA Newtons
z = 25;                        %z set equal to 25 mm
L = 2500;                      %l set equal to 2500 mm
Ks = 1.35;                     %spring stiffness set equal to 1.35 N/mm
delta_w = -7;                  %Sets step value for change in W

%First step define variables
j = 1;                           %Sets location in matrix
W_knot = 0;                      %Sets start point of 0
k = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for k
w = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for w
W = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for W
N = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for N
g = zeros(1, 14);                %Create an empty 1 x 14 matrix of zeroes for g


for j = 1:13                         %Do 13 repetitions (91/7)

   W(j+1) = -7*j;
   k(j) = EA/L*(z/L)^2+EA/L^3*(2*z*w(j)+(w(j))^2)+N(j)/L+Ks;        %Calculate step value of tangent stiffness
   w(j+1) = w(j) + delta_w/k(j);                                    %Calculate step value of w
   N(j+1) = EA*((z/L)*w(j+1)/L+1/2*(w(j+1)/L)^2);                   %Calculate step value of N
   g(j) = N(j+1)*((z+w(j+1))/L)+Ks*w(j+1)-W(j+1);                   %Calculate step value of g

   while abs(g(j)) >= 1*10^-4                                       %Begin iteration loops until error below desired value
       k(j) = EA/L*(z/L)^2+EA/L^3*(2*z*w(j+1)+(w(j+1))^2)+N(j+1)/L+Ks;  %Calculate iteration value of tangent stiffness
       dw = -g(j)/k(j);                                                 %Determine displacement induced by unbalance force
       w(j+1) = w(j+1) + dw;                                            %Calculate iteration value of w
       N(j+1) = EA*((z/L)*w(j+1)/L+1/2*(w(j+1)/L)^2);                   %Calculate iteration value of N
       g(j) = N(j+1)*((z+w(j+1))/L)+Ks*w(j+1)-W(j+1);                   %Calculate iteration value of g
       
   end

   j = j+1;                                 %Increment j and repeat

end

W_abs = abs(W);                             %Take absolute value to plot on negative scale
w_abs = abs(w);                             %Take absolute value to plot on negative scale

figure(4)
clf
plot(w_abs,W_abs,'k')
grid on
ylim([0 120])
xlim([0 60])
xlabel('-w')
ylabel('-W')

hold on
plot(x,W_ex_abs,'b--')
grid on
ylim([0 120])
xlim([0 60])
xlabel('-w')
ylabel('-W')
legend ('Newton-Raphson', 'Exact')