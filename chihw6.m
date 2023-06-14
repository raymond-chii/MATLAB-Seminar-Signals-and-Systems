clc
clear
close all
%%%% Lei(Raymond) Chi HW6

%% Butterworth
fs = 44100; % Sampling frequency
fpass1 = fs/6;
fpass2 = fs/3; % Passband frequencies
fstop1 = fs/7;
fstop2 = fs/2.5; % Stopband frequencies

Rpass = 1; % Passband ripple
Rstop = 50; % Stopband attenuation

H  = fdesign.bandpass(fstop1, fpass1, fpass2, fstop2, Rstop, Rpass, Rstop, fs);
Hd = design(H, 'butter', 'MatchExactly', 'stopband');



% Plotting frequency response
[H, w] = freqz(Hd); % Frequency response
figure;
plot(w/pi*fs/2,20*log10(abs(H)));
title('Butterworth frequency Response');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
grid on;


% Filter a Gaussian white noise signal
t = 0:1/fs:2; % Time vector
x = randn(size(t)); % Gaussian white noise
y = filter(Hd, x); % Filtered signal


% Plot the filtered signal 
n = length(y); % Fourier transform
n = floor(n);
Y = fft(y)/n;
f = fs*(0:n/2)/n; % Frequency vector

figure;
plot(f, 2*abs(Y(1:n/2+1)));
title('Filtered Signal Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;

%% Elliptic
fpass1 = fs/7;
fpass2 = fs/2.5; % Passband frequencies
fstop1 = fs/6;
fstop2 = fs/3; % Stopband frequencies
Rpass = 1; % Passband ripple
Rstop = 50; % Stopband attenuation
[n, Wn] = ellipord([fpass1 fpass2]/(fs/2), [fstop1 fstop2]/(fs/2), Rpass, Rstop);% Filter order and cutoff frequencies
[a, b] = ellip(n, Rpass, Rstop, Wn, 'stop'); % Filter coefficients

% Plotting frequency response
N = 1024;
[H, w] = freqz(b, a, N, fs);
figure;
plot(w/(2*pi)*fs, 20*log10(abs(H)));
grid on;
title('Elliptic filter frequency response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');


% Gaussian white noise
t = 0:1/fs:2; % Time vector
x = randn(size(t)); % Gaussian white noise
y = filter(a, b, x);

% Fourier transform
figure;
f = fftshift(linspace(-fs/2, fs/2, length(y)));
plot(f, 20*log10(abs(fft(y)/length(y))));
grid on;
title('Filtered Signal with Fourier transform');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%% Chebyshev Type I
fpass = fs/9; % Passband frequency
fstop = fs/8; % Stopband frequency
Rpass = 5; % Passband ripple
Rstop = 40; % Stopband attenuation
[n, Wn] = cheb1ord(fpass/(fs/2), fstop/(fs/2), Rpass, Rstop); % Filter order and cutoff frequency
[a, b] = cheby1(n, Rpass, Wn, 'low'); % Filter coefficients
[H, w] = freqz(b, a, 1024, fs);
% Plotting frequency response
figure;
plot(w/(2*pi)*fs, 20*log10(abs(H)));
grid on;
title('Chebyshev Type I Filter Frequency response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');


% Applying filter to Gaussian white noise
t = 0:1/fs:2; % Time vector
x = randn(size(t)); % Gaussian white noise
y = filter(a, b, x);

% Plotting Fourier transform of filtered signal
figure;
f = fftshift(linspace(-fs/2, fs/2, length(y)));
plot(f, 20*log10(abs(fft(y)/length(y))));
grid on;
title('Filtered Signal with Fourier transform');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');




%% Chebyshev Type II
fpass = fs/3; % Passband frequency
fstop = fs/4; % Stopband frequency
Rpass = 5; % Passband ripple
Rstop = 40; % Stopband attenuation
[n, Wn] = cheb2ord(fpass/(fs/2), fstop/(fs/2), Rpass, Rstop); % Filter order and cutoff frequency
[a, b] = cheby2(n, Rstop, Wn, 'high'); % Filter coefficients
% frequency response
[H, w] = freqz(a, b, 1024, fs);
figure;
plot(w/(2*pi)*fs, 20*log10(abs(H)));
grid on;
title('Chebyshev Type II Filter Frequency response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');


% Applying filter to Gaussian white noise
t = 0:1/fs:2; % Time vector
x = randn(size(t)); % Gaussian white noise
y = filter(b, a, x);

% Plotting frequency response
figure;
f = fftshift(linspace(-fs/2, fs/2, length(y)));
plot(f, 20*log10(abs(fft(y)/length(y))));
grid on;
title('Filtered Signal with Fourier transform');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');



%% Play back the unfiltered and filtered signals
sound(x, fs); % Unfiltered signal
sound(y, fs); % Filtered signal
