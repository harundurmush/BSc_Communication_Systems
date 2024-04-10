clc;
clear all;
close all;
%% 1.1 Sampling
%%
f1 = 2500;
f2 = 6000;
Fs = 80000;
Td = 0.008;
Ts = 1/Fs;
nvec = 0:Ts:(Td-Ts);
x1_n = cos(2*pi*f1.*nvec);
x2_n = sin(2*pi*f2.*nvec);
x_n = x1_n + x2_n;
nvec4 = 0:(Ts*4):(Td-Ts);
nvec12 = 0:(Ts*12):(Td-Ts);

x4_n = x_n(mod((1:length(nvec)), 4) == 1);
x12_n = x_n(mod((1:length(nvec)), 12) == 1);

figure(1)

subplot(311);
stem(nvec, x_n );
title('Original Signal x[n]');
ylabel('Amplitude');
xlabel('Discrete time (n)');
legend('x[n] with ds factor 1 (original)');

subplot(312);
stem(nvec4, x4_n );
title('Downsampled x[n] by 4');
ylabel('Amplitude');
xlabel('Discrete time (n)');
legend('x[n] with ds factor 4');

subplot(313);
stem(nvec12, x12_n);
title('Downsampled x[n] by 12');
ylabel('Amplitude');
xlabel('Discrete time (n)');
legend('x[n] with ds factor 12');

N = length(nvec);
N_4 = length(nvec4);
N_12= length(nvec12);
 
FVec = linspace(-Fs/2, Fs/2, N);
FVec4 = linspace((-Fs/2)/4, (Fs/2)/4, N);
FVec12 = linspace((-Fs/2)/12, (Fs/2)/12, N);
xf = abs(fftshift(fft(x_n,N)))/N;
xf4 = abs(fftshift(fft(x4_n,N)))/N_4;
xf12 = abs(fftshift(fft(x12_n,N)))/N_12;

figure(2)

subplot(311);
plot(FVec, xf );
title('Original Signal x[n]');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('xf');

subplot(312);
plot(FVec4, xf4 );
title('Downsampled x[n] by 4');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('xf4');

subplot(313);
plot(FVec12, xf12);
title('Downsampled x[n] by 12');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('xf12');

x4n_lin = interp1(nvec4,x4_n,nvec,'linear');
x4n_cub = interp1(nvec4,x4_n,nvec,'cubic');
x12n_lin = interp1(nvec12,x12_n,nvec,'linear');
x12n_cub = interp1(nvec12,x12_n,nvec,'cubic');

figure(3)

subplot(211)
plot(nvec,x_n);
hold on;
plot(nvec,x4n_lin);
hold on;
plot(nvec,x4n_cub);
title('Original and reconstructed signals x[n] and linear and cubic interpolated x4[n]');
ylabel('Amplitude');
xlabel('Discrete time (n)');
legend('x[n]','linear x4[n]','cubic x4[n]');

subplot(212)
plot(nvec,x_n);
hold on;
plot(nvec,x12n_lin);
hold on;
plot(nvec,x12n_cub);
title('Original and reconstructed signals x[n] and linear and cubic interpolated x12[n]');
ylabel('Amplitude');
xlabel('Discrete time (n)');
legend('x[n]','linear x12[n]','cubic x12[n]');

%% 1.2 Quantization
%%
f3 = 2000;
f4 = 400;
t = 0:Ts:(Td-Ts);
xt = 3*cos(2*pi*f3.*t) + sin(2*pi*f4.*t);
b = max(xt);
a = min(xt);
N4 = 4;
N6 = 6;
x_quant4 = floor(((xt-a)/(b-a))*(2^N4-1))*((b-a)/(2^N4-1)) + a;
x_quant6 = floor(((xt-a)/(b-a))*(2^N6-1))*((b-a)/(2^N6-1)) + a;

figure(4)

subplot(211);
plot(t, xt);
hold on;
plot(t, x_quant4);
title('Original x2[n] and quantized xq4(t) signal');
ylabel('Amplitude');
xlabel('Time (s)');
legend('x(t)','x_q4(t)');

subplot(212);
plot(t, xt);
hold on;
plot(t, x_quant6);
title('Original x2[n] and quantized xq6(t) signal');
ylabel('Amplitude');
xlabel('Time (s)');
legend('x(t)','x_q6(t)');

e4 = xt - x_quant4;
e6 = xt - x_quant6;

SQNR4_linear = var(xt)/var(e4);
SQNR6_linear = var(xt)/var(e6);
SQNR4_dB = 10*log10(SQNR4_linear)
SQNR6_dB = 10*log10(SQNR6_linear)
