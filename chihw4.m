clc
clear
close all
%%%% Lei(Raymond) Chi HW4

n = 5; % number of vectors
m = 7; % length of each vector
A = rand(m, n) + 1i * rand(m, n);

G = gram_schmidt(A)
O = orthonormal(G)

g = rand(m, 1) + 1i * rand(m, 1);

proj = ortho_proj(g , G)

x = linspace(0, 2*pi, 1000);
y = sin(x); 

u = [0, pi/2, pi, 3*pi/2, 2*pi];

[x, u] = ndgrid(x, u);
for j = 1:5
    i = 1/sqrt(2*pi)*exp(-(x-u(:,j)).^2/2);
    plot(x(:,j), i)
    hold on
end

xlabel('x')
ylabel('i(x,u)')
legend('u=0', 'u=pi/2', 'u=pi', 'u=3pi/2', 'u=2pi')




U = gram_schmidt(i);
c = ortho_proj(y', U);
y_est = (U' * c)';

figure;

% Upper plot: Original and estimated sinusoid
subplot(2, 1, 1);
plot(x, y, 'b-', 'LineWidth', 1, 'DisplayName', 'Original');
hold on;
plot(x, y_est, 'r--', 'LineWidth', 1, 'DisplayName', 'Estimated');
xlabel('x');
ylabel('y');
legend('show');
xlim([0, pi/2]);  % Adjust the plot range
title('Original and Estimated Sinusoid');

% Lower plot: Orthonormal basis functions
subplot(2, 1, 2);
plot(x, U, 'LineWidth', 2);
xlabel('x');
ylabel('y');
title('Orthonormal Basis Functions');
legend('Vector 1', 'Vector 2', 'Vector 3', 'Vector 4', 'Vector 5');

function G = gram_schmidt(A)
[m, n] = size(A);
G = zeros(m, n);
    for i = 1:n
        v = A(:,i);
        for j = 1:i-1
            q = G(:,j);
            v = v - dot(q,v)*q;
        end
        G(:,i) = v / norm(v);
    end
end

function O = orthonormal(A)
[m, n] = size(A);
O = true(m,n); 
for j = 1:n
    if abs(norm(A(:,j)) - 1) > eps
        O(:,j) = false;
        return
    end
end
for j = 1:n
    for k = 1:j-1
        if abs(dot(A(:,k), A(:,j))) > eps
            O(:,k) = false;
            O(:,j) = false; 
            return
        end
    end
end
end

function proj = ortho_proj(v , M)
proj = M * (M' * v);
end


