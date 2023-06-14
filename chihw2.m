clc
clear
close all
%%%% Lei(Raymond) Chi HW2

%% space
Z = exp(0:63);
A = reshape(Z, [8,8])';
B = A(3:6,3:6);
m = nthroot(prod(B(:)),8);

newB = B;
newB(end, end) = 0;
newB = circshift(newB, [2, 2]);
A(3:6, 3:6) = A(3:6, 3:6) - newB;
A(:, 1) = flipud(A(:, 1));
A(4:5, :) = [];

%% speed

X = linspace(0, 1000*pi, 10000);

% 1. for loops and no preallocation.
tic
for a = 1:10000
    for b = 1:10000
        C(a, b) = sin(abs(X(a)+1i*X(b)))/abs(X(a)+1i*X(b));
    end
end
toc
% 477.642150 seconds.

surf(X(1:50), X(1:50), C(1:50, 1:50))
title('pain')
xlabel('pain horizontally')
ylabel('pain vertically')
zlabel('pain up and down')

%2.for loops, but preallocating C with zeroes.
C = zeros(10000, 10000);
tic
for a = 1:10000
    for b = 1:10000
C(a, b) = sin(abs(X(a)+1i*X(b)))/abs(X(a)+1i*X(b));
    end
end
toc 
% 16.599701 seconds.

surf(X(1:50), X(1:50), C(1:50, 1:50))
title('pain x 2')
xlabel('pain x 2 horizontally')
ylabel('pain x 2 vertically')
zlabel('pain x 2 up and down')

%3.meshgrid.
[I, J] = meshgrid(X);
tic
C = sin(abs(I+1i.*J))./abs(I+1i.*J);
toc 
% 2.572164 seconds.

surf(X(1:50), X(1:50), C(1:50, 1:50))
title('pain x 3')
xlabel('pain x 3 horizontally')
ylabel('pain x 3 vertically')
zlabel('pain x 3 up and down')

%4. broadcasting.
tic
C = sin(abs(X+1i.*(X.')))./abs(X+1i.*(X.'));
toc 
% 3.001460 seconds.

surf(X(1:50), X(1:50), C(1:50, 1:50))
title('pain x 4')
xlabel('pain x 4 horizontally')
ylabel('pain x 4 vertically')
zlabel('pain x 4 up and down')

% meshgrid is the most op, and for looping without preallocation is just
% painful just to wait

%% the long sssss
n = 100;
t1 = linspace(0, 6.66, n);
v1 = exp(-t1.^2);

n = 10000;
t2 = linspace(0, 6.66, n);
v2 = exp(-t2.^2);

%approx the vecs
dv1 = diff(v1);
dv2 = diff(v2);
error_v1 = (dv1 + 2 * t1(1:99) .* v1(1:99)).^2;
MSE_v1 = double(mean(error_v1));
error_v2 = (dv2 + 2 * t2(1:9999) .* v2(1:9999)).^2;
MSE_v2 = double(mean(error_v2));

% integrals with cumsum
v1_sum = cumsum(v1)*(6.66/100)*2/sqrt(pi);
v2_sum = cumsum(v2)*(6.66/10000)*2/sqrt(pi);

% integrals with cumtrapz
v1_trapz = cumtrapz(v1)*(6.66/100)*2/sqrt(pi);
v2_trapz = cumtrapz(v2)*(6.66/10000)*2/sqrt(pi);

% mean squared error with cumsum
error_cumsum1 = (v1_sum - erf(t1)).^2;
MSE_sum_v1 = mean(error_cumsum1)
error_cumsum2 = (v2_sum - erf(t2)).^2;
MSE_sum_v2 = mean(error_cumsum2)

% mean squared error wit cumtrapz
error_cumtrapz1 = (v1_trapz - erf(t1)).^2;
MSE_trapz_v1 = mean(error_cumtrapz1)
error_cumtrapz2 = (v2_trapz - erf(t2)).^2;
MSE_trapz_v2 = mean(error_cumtrapz2)

plot(t2, v2_trapz)
tittle('FInally finished')
