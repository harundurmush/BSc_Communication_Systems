clc;
clear all;
close all;
%% 2.1 Modulation
%% 2.1.a
Td=0.2;
Fs=2000;
Ts=1/Fs;
time=0:Ts:Td;
%% 2.1.b
Ac=1;
fc=200;
ct=Ac*cos(2*pi*fc*time);
%% 2.1.c
fm1=50;
Am1=1;
m1t=Am1*cos(2*pi*fm1*time);
fm2=10;
Am2=1;
m2t=Am2*cos(2*pi*fm2*time);
%% 2.1.d
% by using the following PM formula: Spm(t)=Ac*cos(2*pi*fc*time + kp*m(t))
kp1=pi/4;
s1t=Ac*cos(2*pi*fc*time+kp1*m1t);
s2t=Ac*cos(2*pi*fc*time+kp1*m2t);
%% 2.1.e
m3t=m1t+m2t;
s3t=Ac*cos(2*pi*fc*time+kp1*m3t);
%% 2.1.f
s12t=s1t+s2t;
figure (1)
subplot(2,1,1)
plot(time,s3t,"color","#7E2F8E");
xlabel("time (s)");
ylabel("voltage (V)");
legend("s3(t)");
title("Comparison of s3(t) vs. sum of s1(t) and s2(t)");
subplot(2,1,2)
plot(time,s12t,"color","#77AC30");
xlabel("time (s)");
ylabel("voltage (V)");
legend("s1(t)+s2(t)");
%% 2.2 The effect of different kp values
%% 2.2.a
spmt=s1t;
kp2=2*pi/4;
spm1t=Ac*cos(2*pi*fc*time+kp2*m1t);
kp3=3*pi/4;
spm2t=Ac*cos(2*pi*fc*time+kp3*m1t);
%% 2.2.b
figure (2)
subplot(4,1,1)
plot(time,m1t,"color","#0072BD");
xlabel("time (s)");
ylabel("voltage (V)");
legend("m1(t)");
title("Effect of different kp values on the phase modulated signals by modulating the message signal m1(t)");
subplot(4,1,2)
plot(time,spmt,"color","#A2142F");
xlabel("time (s)");
ylabel("voltage (V)");
legend("spm(t)");
subplot(4,1,3)
plot(time,spm1t,"color","#D95319");
xlabel("time (s)");
ylabel("voltage (V)");
legend("spm1(t)");
subplot(4,1,4)
plot(time,spm2t,"color","#EDB120");
xlabel("time (s)");
ylabel("voltage (V)");
legend("spm2(t)");
%% 2.3 Demodulation
%% 2.3.a
N=length(time);
% spm(t) demodulation
spmht=hilbert(spmt,N);
phase_ins=angle(spmht);
unwrapped=unwrap(phase_ins);
pm_demod=(unwrapped-2*pi*fc*time)/kp1;
% spm1(t) demodulation
spm1ht=hilbert(spm1t,N);
phase_ins1=angle(spm1ht);
unwrapped1=unwrap(phase_ins1);
pm_demod1=(unwrapped1-2*pi*fc*time)/kp2;
% spm2(t) demodulation
spm2ht=hilbert(spm2t,N);
phase_ins2=angle(spm2ht);
unwrapped2=unwrap(phase_ins2);
pm_demod2=(unwrapped2-2*pi*fc*time)/kp3;
%% 2.3.b
figure (3)
subplot(3,1,1)
plot(time,m1t,"color","#0072BD");
title("Comparison of original m(t) vs demodulated signal spm(t),spm1(t) and spm2(t)");
hold on
plot(time,pm_demod,"color","#A2142F");
xlabel("time (s)");
ylabel("voltage (V)");
legend("m1(t)","demodulated spm(t)");
subplot(3,1,2)
plot(time,m1t,"color","#0072BD");
hold on
plot(time,pm_demod1,"color","#D95319");
xlabel("time (s)");
ylabel("voltage (V)");
legend("m1(t)","demodulated spm1(t)");
subplot(3,1,3)
plot(time,m1t,"color","#0072BD");
hold on
plot(time,pm_demod2,"color","#EDB120");
xlabel("time (s)");
ylabel("voltage (V)");
legend("m1(t)","demodulated spm2(t)");















