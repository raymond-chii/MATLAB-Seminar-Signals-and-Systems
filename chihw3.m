clc
clear
close all
%%%% Lei(Raymond) Chi HW3

%%Space between 

Z = reshape(exp(0:63),[8,8])'; Z(3:6,3:6) = 0;
nonzero_terms = Z(Z ~= 0);
geomean_nonzero = prod(nthroot(nonzero_terms,2));


A = sin(linspace(0,5,100).*linspace(-5,0,100).');[val, ind] = min(abs(A - 1/2)); [min_val, ind2] = min(val); 
value = A(ind(ind2), ind2);
index = [ind(ind2), ind2];


f1 = @(x,y) exp(-(1-x.*y).^2);
f2 = @(x,y) 0.25*(x.^2 + y.^2).^(1/2);
V = integral2(@(x,y) max(0, f1(x,y) - f2(x,y)), -1, 1, -1, 1);
fsurf(@(x,y) f1(x,y), [-1, 1, -1, 1]);
hold on;
fsurf(@(x,y) f2(x,y), [-1, 1, -1, 1]);


%%I need a vacation 

% create matrix A: a circle that is positioned 99,99 with radius 29
A = false(256,256,1);
for i=1:256
    for j=1:256
        if sqrt((i-99)^2+(j-99)^2)<29
            A(i,j) = true;
        end
    end
end

figure
imshow(A)

% create matrix B: a circle that is positioned 62,62 with radius 58
B = false(256,256,1);
for i=1:256
    for j=1:256
        if sqrt((i-62)^2+(j-62)^2)< 58
            B(i,j) = true;
        end
    end
end
figure
imshow(B)

% create matrix C : the area underneath the sin wave
C = false(256,256,1);
for i=1:256
    for j=1:256
        if i-4*sin(j/10)>100
            C(i,j) = true;
        end
    end
end
figure;
imshow(C);

% create matrix S: random signals 
S=rand(256,256,3);
figure
imshow(S)

% matrix M: intersection of a and b
M = A & B; 
figure
imshow(M)
% matrix Z

Z = double(C) .* S + double(M);

figure
imshow(Z)


%%My sinc is broken 

x = linspace(-2*pi, 2*pi, 1000);
y = sinc(x);
dydx = deriv(y, x);
int_y = antideriv(y, x);


[x_extrema, y_extrema] = extrema(y, x);
[x_inflections, y_inflections] = inflections(y, x);
plot(x, y, x, dydx, x, int_y, x_extrema, y_extrema, 'r*', x_inflections, y_inflections, 'bo');
legend('sinc(x)', 'sinc''(x)', 'integral of sinc(x)', 'local extrema', 'inflection points');
xlabel('x');
ylabel('y');

function y = sinc(x)
    y = sin(x)./x;
    y(x==0) = 1;
end

function dydx = deriv(y, x)
    dydx = diff([y(1), y]) ./ diff([x(1), x]); 
end

function y_int = antideriv(y, x)
    y_int = cumtrapz(x, y);
end

function out = switchsign(g)
    out = abs(sign(g) - sign([0, g(1:end-1)])) == 2;
end

function [x_extrema, y_extrema] = extrema(y, x)
    dydx = deriv(y,x); 
    [peaks, locals] = max(abs(diff(sign(dydx)))); 
    x_extrema = x(locals+1-peaks:locals+peaks); y_extrema = y(locals+1-peaks:locals+peaks);
end

function [x_inflections, y_inflections] = inflections(y, x)
    d2ydx2 = deriv(deriv(y, x), x); 
    [x_inflections, y_inflections] = deal(x(switchsign(d2ydx2)), y(switchsign(d2ydx2)));
end
