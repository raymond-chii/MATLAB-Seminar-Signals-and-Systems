clc
clear
close all
%%%% Lei(Raymond) Chi HW3

%%be discrete

num = [12, 23, 37, 0];
den = [2, 13, 0, 12];
zplane(num, den);

title('Pole-Zero Plot of H(z)');
xlabel('Re');
ylabel('Im');

[z, p, k] = tf2zpk(num, den);

impulse_response = impz(num, den, 66);
n = 0:65;
stem(n, impulse_response);
title('Impulse Response of H(z)');
xlabel('Index');
ylabel('Amp');

n = 1:55;
x = sin(n/10);

% Filter using MATLAB's filter function
y1 = filter(num, den, x);

% Filter using convolution with the impulse response
y2 = conv(x, impulse_response);

figure;
subplot(2,1,1);
stem(n, y1);
title('Output of H(z) using filter');
xlabel('Index');
ylabel('Amp');

subplot(2,1,2);
stem(n, y2(1:55));
title('Output of H(z) using convolution');
xlabel('Index');
ylabel('Amp');


%%Poley moley

load handel.mat
soundsc(y, Fs);


k = 0.1;
p = [0.74*exp(0.76j), 0.74*exp(-0.76j), 0.96*exp(1.24j), 0.96*exp(-1.24j)];
z = [exp(2.06j), exp(-2.06j), exp(1.43j), exp(-1.43j)];

b = k * poly(z);
a = poly(p);


[H, w] = freqz(b, a);

figure;
subplot(2,1,1);
plot(w/pi, abs(H));
title('Magnitude');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Mag (dB)');
grid on;

subplot(2,1,2);
plot(w/pi, unwrap(angle(H)));
title('Phase Response');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (rad)');
grid on;

%% bonus 
figure;
zplane(b, a);
title('Pole-Zero Plot');

%question 4
y_filtered = filter(b, a, y);
soundsc(y_filtered, Fs);

%% bonus 
z_new = [exp(2.06j), exp(-2.06j), exp(1.43j), exp(-1.43j), exp(pi/3*j), exp(pi/4*j)];
b_new = k * poly(z_new);
a_new = poly(p);

[H_new, w_new] = freqz(b_new, a_new);

