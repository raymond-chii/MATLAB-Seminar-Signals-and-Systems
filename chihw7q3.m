clc
clear
close all
%%%% Lei(Raymond) Chi HW7q3
syms x y z t a; 
F = [(3 + a^2)*x*z, z*exp(y), exp(y) - (x^2)*exp(pi*a)];
r = [(1 - 2*cos(t))*cos(3*t), (1 - 2*cos(t))*sin(3*t), sin(t)]; 
F_r = subs(F, [x, y, z], r); 
r_pr = diff(r, t); 
int_result = int(F_r * r_pr.', t, -pi, pi); 
int_function = matlabFunction(int_result, 'Vars', {a});

a_vals = linspace(-3, 1, 100);

int_results = arrayfun(int_function, a_vals);

figure;
fplot(int_function, [-3, 1])
xlabel('a')
ylabel('Integral of F(r(t)) · r''(t)')
title('Integral of F(r(t)) · r''(t) function of a')



% Bonus
curl_F = curl(F, [x, y, z]);
is_conservative = isequal(curl_F, [0, 0, 0]); 
a_for_conservative = solve(is_conservative, a);

if isempty(a_for_conservative)
    disp('no value of a that F is conservative.')
else
    disp(['F is conservative for a = ', char(a_for_conservative)])
end