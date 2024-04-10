clc;
clear all;
%% step 1.0
Fs=20000;
t=0:1/Fs:0.5;
f1=20;
f2=50;
f3=1000;
f4=5000;
xt=5*sin(2*pi*f1.*t)+cos(2*pi*f2.*t)+15*cos(2*pi*f3.*t)+10*cos(2*pi*f4.*t);
N=length(xt);
f=linspace(-Fs/2,Fs/2,N);
xf=fftshift(abs(fft(xt,N)))/N;
%% step 1.1
figure (1)
subplot(211)
plot(t,xt,"blue");
xlabel("time (s)");
ylabel("x(t)");
legend("x(t)");
title("Original signal x(t)");
subplot(212)
plot(f,xf,"red");
xlabel("frequency (Hz)");
ylabel("Magnitude of x(f) (db)");
legend("x(f)");
title("Frequency Response of the signal x(t), x(f)");
%% step 1.2
figure(2)
subplot(311)
[b,a]= butter(1,40/(Fs/2),"low"); 
% 40 cut-off freq is used not to include 50 Hz in frequency response of the
% signal that might cause an error ( two lines one in 20 Hz other one in 50
% Hz ) that occurs in last step
h1=freqz(b,a,Fs/2);
plot(abs(h1),"blue");
hold on
[bl,al]= butter(7,40/(Fs/2),"low");
h1=freqz(bl,al,Fs/2);
plot(abs(h1),"red");
xlabel("frequency (Hz)");
ylabel("|H(f)|");
legend("1th order","7th order");
title ("LPF (Low-pass) filter");
subplot(312)
[b,a]= butter(1,2000/(Fs/2),"high");
h1=freqz(b,a,Fs/2);
plot(abs(h1),"blue");
hold on
[bh,ah]= butter(7,2000/(Fs/2),"high");
h2=freqz(bh,ah,Fs/2);
plot(abs(h2),"red");
xlabel("frequency (Hz)");
ylabel("|H(f)|");
legend("1th order","7th order");
title ("HPF (High-pass) filter");
subplot(313)
[b,a]= butter(1,[400 2000]/(Fs/2),"bandpass");
h1=freqz(b,a,Fs/2);
plot(abs(h1),"blue");
hold on
[bb,ab]= butter(7,[400 2000]/(Fs/2),"bandpass");
h3=freqz(bb,ab,Fs/2);
plot(abs(h3),"red");
xlabel("frequency (Hz)");
ylabel("|H(f)|");
legend("1th order","7th order");
title ("BPF (Band-pass) filter");
%% step 1.3
xLPF=filter(bl,al,xt);
xHPF=filter(bh,ah,xt);
xBPF=filter(bb,ab,xt);

xLPFf=fftshift(abs(fft(xLPF,N)))/N;
xHPFf=fftshift(abs(fft(xHPF,N)))/N;
xBPFf=fftshift(abs(fft(xBPF,N)))/N;

% for better visualization, expand figure 3 window both vertically and
% horizontally

figure (3)
subplot(221)
plot(f,xf,"blue");
xlabel("frequency (Hz)");
ylabel("Magnitude of x(f) (db)");
legend("x(f)");
title("Frequency Response of the signal x(t), x(f)");
% I thought it makes more sense if I print frequency response of the
% original signal to compare them with other filters. I changed it to
% frequency response graph instead of time domain.
subplot(222)
plot(f,xLPFf,"red");
xlabel("frequency (Hz)");
ylabel("Magnitude of xLPF(f) (dB)");
legend("xLPF(f)");
title("Low-pass Filter Frequency Response of the signal x(t)");
subplot(223)
plot(f,xHPFf,"red");
xlabel("frequency (Hz)");
ylabel("Magnitude of xHPF(f) (dB)");
legend("xHPF(f)");
title("High-pass Filter Frequency Response of the signal x(t)");
subplot(224)
plot(f,xBPFf,"red");
xlabel("frequency (Hz)");
ylabel("Magnitude of xBPF(f) (dB)");
legend("xBPF(f)");
title("Band-pass Filter Frequency Response of the signal x(t)");
%%