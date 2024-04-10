clc;
clear all;
close all;
%% DSB-SC Modulation and Demodulation
%% a
t = 0.1;
fm = 50;
fc = 1000;
Fs = 50000;
teta = 0; % phase = 0 will be used at future steps
Ts = 1/Fs;
time = 0:Ts:t;
Ac = 1;
Am = 2;
mt = Am*cos(2*pi*fm*time);
ct = Ac*cos(2*pi*fc*time+teta);
%% b
st= mt.*ct;
%% c and d
figure (1)
subplot(211)
plot(time,mt,"blue");
hold on
plot(time,ct,"red");
xlabel("time (s)");
ylabel("m(t) and c(t) (V)");
title("Original message signal m(t) and carrier signal c(t)");
subplot(212)
plot(time,st,"cyan");
xlabel("time (s)");
ylabel("s(t) (V)");
title("Modulated signal s(t)");
%% e
N = length(time);
f = linspace(-Fs/2,Fs/2, N);
Acprime = max(st); 
% Acprime is maximum amplitude of modulated signal.
multiplier = Acprime*cos(2*pi*fc*time+teta);
ut = st.*multiplier;
[order, wcutoff] = buttord(fm/(Fs/2),fc/(Fs/2),1,32);
% stopband attenuation is chosen in order to reduce ripple and phase delay
% to minimum level.
[b,a]=butter(order,wcutoff,"low");
yt = filter(b,a,ut);
uf = fftshift(abs(fft(ut,N)))/N;
mf = fftshift(abs(fft(mt,N)))/N;
sf = fftshift(abs(fft(st,N)))/N;
%% f
figure (2)
subplot(311)
plot(f,mf,"blue");
xlabel("frequency (Hz)");
ylabel("|M(f)| (dB)");
title("Magnitude response of m(t)");
subplot(312)
plot(f,sf,"cyan");
xlabel("frequency (Hz)");
ylabel("|S(f)| (dB)");
title("Magnitude response of s(t)");
subplot(313)
plot(f,uf,"magenta");
xlabel("frequency (Hz)");
ylabel("|U(f)| (dB)");
title("Magnitude response of u(t)");
%% g
figure (3)
plot(time,mt,"blue");
xlabel("time (s)");
ylabel("m(t) and y(t) (V)");
title("Original message signal m(t) vs demodulated signal y(t)");
hold on
plot(time,yt,"red");
legend("m(t)","y(t)");
%% SSB Modulation
%% a
[b1,a1] = butter(5,[1000 2000]/(Fs/2),"bandpass");
%there was an error in this section during lab session, it is fixed at
%report.
%% b
sut = filter(b1,a1,st);
%% c
[h,w]=freqz(b1,a1,Fs/2);
suf = fftshift(abs(fft(sut,N)))/N;
figure (4)
subplot(211)
plot(abs(h),"blue");
xlabel("frequency (Hz)");
ylabel("|BPF| (dB)");
title("Magnitude response of band-pass filter");
subplot(212)
plot(f,suf,"red");
xlabel("frequency (Hz)");
ylabel("|Su(f)| (dB)");
title("Magnitude response of filtered & modulated signal su(t)");






