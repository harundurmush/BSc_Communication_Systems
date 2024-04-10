clc;
clear all;
A1=1;
A2=2;
Ac=2;
f1=10;
f2=15;
fc=500;
Fs=4000;
Ts=1/Fs;
Tduration=0.4;
t=0:Ts:Tduration;
mt=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t);
ct=Ac*sin(2*pi*fc*t);
%% Conventional Amplitude Modulation
%% a
figure (1)
plot(t,mt,"black");
grid on
xlabel("time (s)");
ylabel("m(t)");
title("Original message signal m(t)");
%% b
ka1=0.2;
ka2=0.6;
st1=(1+ka1*mt).*ct;
st2=(1+ka2*mt).*ct;

figure (2)
subplot(211)
plot(t,st1,"blue");
xlabel("time (s)");
ylabel("s1(t)");
title("Modulated signal with ka=0.2");
subplot(212)
plot(t,st2,"red");
xlabel("time (s)");
ylabel("s2(t)");
title("Modulated signal with ka=0.6");
%% c
N=length(t);
sf1=fftshift(abs(fft(st1,N)))/N;
sf2=fftshift(abs(fft(st2,N)))/N;
f=linspace(-Fs/2,Fs/2,N);

figure (3)
subplot(211)
plot(f,sf1,"blue");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("ka=0.2");
title("Magnitude response of the modulated signal with ka=0.2");
subplot(212)
plot(f,sf2,"red");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("ka=0.6");
title("Magnitude response of the modulated signal with ka=0.6");

%comment eklenecek
%% Square-Law Envelope Detector
%% a
term1=st1.^2;
term2=st2.^2;
%% b
fcutoff=150; 
% cut-off frequency is choosen according to the equation
% 2*fm < fcutoff < 2*fc also, for minimum phase delay between comparations
% fcutoff=150 Hz is available
[x,y]=butter(5,fcutoff/(Fs/2),"low");
%% c
s1filtered = sqrt(filter(x,y,term1));
s2filtered = sqrt(filter(x,y,term2));
%% d shifted and scaled signals from the step b in conv.AM
mt1prime=(s1filtered-Ac)/ka1;
mt2prime=(s2filtered-Ac)/ka2;
%% e
figure (4)
plot(t,mt,"black");
grid on
hold on
plot(t,mt1prime,"blue");
xlabel("time (s)");
ylabel("Amplitude");
legend("message signal","ka=0.2");
title("Message signal m(t) vs. shifted and scaled signal m''(t)");
figure(5)
plot(t,mt,"black");
grid on
hold on
plot(t,mt2prime,"red");
xlabel("time (s)");
ylabel("Amplitude");
legend("message signal","ka=0.6");
title("Message signal m(t) vs. shifted and scaled signal m''(t)");
%% Modulation Index
%% modulation index is u=ka*Am so,
Am=max(mt); % max value of message signal
u=1;
ka3=u/Am;

st3=(1+ka3*mt).*ct;
term3=st3.^2;
s3filtered = sqrt(filter(x,y,term3));
mt3prime=(s3filtered-Ac)/ka3;

figure(6)
plot(t,mt,"black");
grid on
legend("m(t)");
hold on
plot(t,mt3prime,"green");
xlabel("time (s)");
ylabel("Amplitude");
legend("message signal","ka=0.3461");
title("Message signal m(t) vs. shifted and scaled signal m''(t)");