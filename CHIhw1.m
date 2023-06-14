clc
clear
close all
diary chihw1.m
% Lei(Raymond) Chi HW1


%scalars
a = abs(sin(pi/3) + 1i/sec(-5*pi/3));
l = nthroot(8, 3);
u = sqrt(2*(sum(1:80)/factorial(6)));
m = (imag(floor(log(sqrt(66).^7j)))).^2;

%matrix
A = [a;l;u;m]
F =[a l;u m]
T = transpose(F)
B = inv(T*F);
C = [T F;F T]

%cruelty
mean_B = mean(mean(B))
means_C = mean(C, 2)

%odd types
T+F %makes sense since they have same dimensions 
T+1 %makes sense all entry plus one 
C+A %makes sense every column in C is add with A

%Not what it seems 
k=3;
vec = linspace(0, 10, k)
means_squares = sumsqr(vec) / k

%different Ks in a vector 
k=[5, 10, 300, 1000000]
for i=1:4
means_squares(i) = sumsqr(linspace(0, 10, k(i))) / k(i);
end
means_squares
% 0-10 were broken down more and more as k increases
% and it will eventually converages to 33.333
diary off
