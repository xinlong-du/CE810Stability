% Poskin CE 810 HW2 Problem 2

clear; clc; close all

w_norm = linspace(0,-2.3,1000); % Normalized displacement (w/z)

% Normalized load equation from notes (W/EA*(L/z)^3), k_s = EAz^3/2l^3
W_norm = 3/2*w_norm + 3/2*w_norm.^2 + 1/2*w_norm.^3;

plot(-w_norm,-W_norm)
title('Poskin HW2 Problem 2: Normalized Load-Deflection Curve')
xlabel('-^{w}/_{z}')
ylabel('-(^{W}/_{EA})(^{L}/_{z})^3')