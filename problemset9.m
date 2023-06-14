%%% ECE211 Problem Set 9 Seyeon Park

%% Q1.
%% 1. Euler Method:
% Substituting f(x(t), u(t), t) and g(x(t), u(t), t) with their LTI system 
% representations:
% x(t + dt) = x(t) + dt · (Ax(t) + Bu(t))
% y(t) = Cx(t) + Du(t)

% Representing the discrete-time variables with the notation used in the 
% Euler method:
% x(t + dt) = xde[n + 1]
% x(t) = xde[n]
% u(t) = ude[n]
% y(t) = yde[n]

% Substituting these variables into the Euler method equations:
% xde[n + 1] = xde[n] + dt · (Axde[n] + Bude[n])
% yde[n] = Cxde[n] + Dude[n]

% Thus
% Ade = eye(size(A)) + dt * A
% Bde = dt * B
% Cde = C
% Dde = D

%% 2. Midpoint Method
% Substitute f(x, u, t) into the equation:
% xdm[n + 1] 
% = xdm[n] + dt * (A * (xdm[n] + (dt/2) * (A * xdm[n] + B * [Im Om; Om Om]udm[n]) + B * [Om Om; Om Im]udm[n], ndt + dt/2))
% = xdm[n] + dt * (A * xdm[n] + A * (dt/2) * (A * xdm[n] + B * [Im Om; Om Om]udm[n]) + B * [Om Om; Om Im]udm[n])
% = (I + dt * A + ((dt)^2)/2 * A^2) * xdm[n] + ((dt)^2/2) * A * B * [Im Om; Om Om]udm[n] + dt * B * [Om Om; Om Im]udm[n]

% Comparing with the standard discrete-time LTI model equation:
% xdm[n + 1] = Adm * xdm[n] + Bdm * udm[n]
% ydm[n] = Cdm * xdm[n] + Ddm * udm[n]

% Adm = I + dt * A + ((dt)^2/2) * A^2
% Bdm = ((dt)^2/2) * A * B * [Im Om; Om Om] + dt * B * [Om Om; Om Im]
% Cdm = C
% Ddm = D


%% Q2.
% 1. Euler Method:
% In the Euler method, we have Ade = eye(size(A)) + dt * A. 
% If we consider the Taylor series expansion of e^(At) and truncate it 
% after the linear term, we get e^(At) ≈ I + At. 
% Comparing this approximation to the matrix Ade, we can see that the 
% linear term in the Taylor series approximation corresponds to dt * A. 
% Therefore, we can view Ade as an approximation of e^(Adt) for small dt.

% 2. Midpoint Method:
% If we consider the Taylor series expansion of e^(At) and 
% truncate it after the quadratic term, we get 
% e^(At) ≈ I + At + (1/2) * (At)^2. 
% Comparing this approximation to the matrix Adm, we can see that the 
% quadratic term in the Taylor series approximation corresponds to Adm. 
% Therefore, we can view Adm as an approximation of e^(Adt) 
% for small dt.

clc; clear; close all;

%% Q3.
% (a) 
% Define the given matrices
E = [3 1; 2 1];
A0 = [-0.3 0.4; -0.4 -0.3];

eigenvalues_A0 = eig(A0);

A = E * A0 * inv(E);
eigenvalues_A = eig(A);

disp("Eigenvalues of A0:");
disp(eigenvalues_A0.');

disp("Eigenvalues of A:");
disp(eigenvalues_A.');


% (b)
% Define symbolic variable s
syms s t

% Compute Laplace transform of e^(At)
L = laplace(expm(A*t), t, s);

% Compute inverse Laplace transform to get e^(At) as a symbolic expression
E = ilaplace(L, s, t);

% Convert symbolic expression to a function
E_func = matlabFunction(E);

% (c)
% Set initial condition xd[0]
xd_0 = [2; 1];

% Set sampling rate fs
fs = 10;

% Compute time increment dt
dt = 1/fs;

% Compute the number of time steps
num_steps = 101; % 0 to 100 inclusive

% Initialize matrix to store xd[n]
xd_matrix = zeros(2, num_steps);

% Compute exact trajectory using the matrix exponential
xd_matrix(:, 1) = xd_0; % Set xd[0]

% Compute xd[n] for 1 <= n <= 100
for n = 1:num_steps-1
    t = n * dt; % Compute time at step n
    xd_matrix(:, n+1) = expm(A * t) * xd_0; % Compute xd[n] using the matrix exponential
end

% Display the matrix with columns xd[n]
disp("Matrix xd[n]:");
disp(xd_matrix);

% (d)
% Initialize matrices to store xde[n] and xdm[n]
xde_matrix = zeros(2, num_steps);
xdm_matrix = zeros(2, num_steps);

% Set initial conditions
xde_matrix(:, 1) = xd_0; % xde[0] = xd[0]
xdm_matrix(:, 1) = xd_0; % xdm[0] = xd[0]

Ade = eye(size(A)) + dt * A;
Adm = eye(size(A)) + dt * A + ((dt)^2/2) * A^2;
% Compute xde[n] and xdm[n] for 1 <= n <= 100
for n = 1:num_steps-1
    % Compute xde[n+1] and xdm[n+1]
    xde_matrix(:, n+1) = Ade * xde_matrix(:, n); % xde[n+1] = Adexde[n]
    xdm_matrix(:, n+1) = Adm * xdm_matrix(:, n); % xdm[n+1] = Admxdm[n]
end

% Display the matrices with columns xde[n] and xdm[n]
disp("Matrix xde[n]:");
disp(xde_matrix);

disp("Matrix xdm[n]:");
disp(xdm_matrix);


% Display the matrices with columns xde[n] and xdm[n]
disp("Matrix xde[n]:");
disp(xde_matrix);

disp("Matrix xdm[n]:");
disp(xdm_matrix);

% (e)
% Plot the trajectories in state space
figure;
hold on;

% Plot xd[n]
plot(xd_matrix(1, :), xd_matrix(2, :), 'b-', 'LineWidth', 1.5);

% Plot xde[n]
plot(xde_matrix(1, :), xde_matrix(2, :), 'r--', 'LineWidth', 1.5);

% Plot xdm[n]
plot(xdm_matrix(1, :), xdm_matrix(2, :), 'g:', 'LineWidth', 1.5);

% Set plot labels and legend
xlabel('x_1');
ylabel('x_2');
legend('xd[n]', 'xde[n]', 'xdm[n]');
title('Trajectories in State Space');

hold off;

% (f)
% Compute the maximum absolute error between xd and xde matrices
error_xde = max(abs(xd_matrix - xde_matrix), [], 'all');

% Compute the maximum absolute error between xd and xdm matrices
error_xdm = max(abs(xd_matrix - xdm_matrix), [], 'all');

% Display the maximum absolute errors
disp("Maximum absolute error (xd vs xde):");
disp(error_xde);

disp("Maximum absolute error (xd vs xdm):");
disp(error_xdm);

% (g)
% As xde has bigger maximum error, xde (Midpoint) is better, but it's not a
% meaningful difference.
% xde makes more error in smaller x_1. The linear approximation made by the 
% exact discretization method (xde) may lead to more errors. 
