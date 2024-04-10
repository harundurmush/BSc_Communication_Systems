clc;
clear all;
close all;
%% 5.1 Modulation
%%
P=10;
fm=50;
Am=1;
Fs=2000;
Ts=1/Fs;
T=P/fm;
t=0:Ts:T-Ts;
fc=200;
Ac=1;
%%
mt=Am*sawtooth(2*pi*fm*t);
%%
ct = Ac*cos(2*pi*fc*t);
%%
figure (1)
subplot(211)
plot(t,mt,"blue");
xlabel("time (s)");
ylabel("voltage (V)");
legend("m(t)");
title("Original message signal m(t) and carrier signal c(t)");
subplot(212)
plot(t,ct,"red");
xlabel("time (s)");
ylabel("voltage (V)");
legend("c(t)");
%%
% Threshold value of beta = delta_f/fm and also beta = kf/fm
% to obtain 0.5 and 1.5 times of beta;
% beta1 = delta_f/fm = kf1/fm = 0.5 , kf1 = 25
% beta2 = delta_f/fm = kf2/fm = 1.5 , kf2 = 75
kf1=25;
kf2=75;
%%
kf=[kf1 kf2];
%%
dt=Ts;
st1 = Ac*cos((2*pi*fc*t) + (2*pi*kf(1)*cumsum(mt,2)*dt));
st2 = Ac*cos((2*pi*fc*t) + (2*pi*kf(2)*cumsum(mt,2)*dt));
%%
figure (2)
subplot(211)
plot(t,st1,"blue");
title("FM modulated signals s1(t) and s2(t)");
xlabel("time (s)");
ylabel("voltage (V)");
legend("st1(t)");
subplot(212)
plot(t,st2,"red");
xlabel("time (s)");
ylabel("voltage (V)");
legend("st2(t)");
%% 5.2 Magnitude Frequency Responses
%%
N=length(t);
f=linspace(-Fs/2,Fs/2,N);
%%
mf = abs(fftshift(fft(mt,N)))/N;
cf = abs(fftshift(fft(ct,N)))/N;
%%
figure (3)
subplot(211)
plot(f,mf,"blue");
legend("|M(f)|");
xlabel("frequency (Hz)");
ylabel("magnitude (dB)");
title("Magnitude Frequency Response of m(t) and c(t)");
subplot(212)
plot(f,cf,"red");
legend("|C(f)|");
xlabel("frequency (Hz)");
ylabel("magnitude (dB)");
%%
sf1=fftshift((fft(st1,N)/N));
sf2=fftshift((fft(st2,N)/N));
figure (4)
subplot(211)
plot(f,sf1,"blue");
legend("S1(f)");
xlabel("frequency (Hz)");
ylabel("relative response in (dB)");
title("Frequency Responses of s1(t) and s2(t)");
subplot(212)
plot(f,sf2,"red");
legend("S2(f)");
xlabel("frequency (Hz)");
ylabel("relative response in (dB)");
%% 5.3 Demodulation
%%
mt1_demod = Am*fmdemod(st1,fc,Fs,Am*kf(1));
mt2_demod = Am*fmdemod(st2,fc,Fs,Am*kf(2));
%%
figure (5)
subplot(211)
plot(t,mt,"blue");
title("Original signal m(t) vs. FM demodulated signals");
hold on
plot(t,mt1_demod,"red");
xlabel("time (s)");
ylabel("voltage (V)");
legend("original signal m(t)","demodulated signal with kf1");
subplot(212)
plot(t,mt,"blue");
hold on
plot(t,mt2_demod,"red");
xlabel("time (s)");
ylabel("voltage (V)");
legend("original signal m(t)","demodulated signal with kf2");
%%












