clc;
clear all;
close all;
%% 8.1 Image Reading
%%
I = imread('testpat1.png');
y = im2double(I);
%%
Fs = size(y,1)*size(y,2);
%% 8.2 Modulation
%%
mt = reshape(y,1,[]);
%%
N = length(mt);
mf = fftshift(abs(fft(mt,N)))/N;
f = linspace(-Fs/2,Fs/2, N);

figure (1)
plot(f,mf);
title("Magnitude of the Frequency Spectrum of Original Image");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("M(f)");
%%
t = 0:(1/Fs):(numel(mt)-1)/Fs;
%%
Ac = 1;
fc = 10000;
dt=1/Fs;
kf1 = 100;
kf2 = 400;
kf = [kf1 kf2];
st1 = Ac*cos((2*pi*fc*t) + (2*pi*kf(1)*cumsum(mt,2)*dt));
st2 = Ac*cos((2*pi*fc*t) + (2*pi*kf(2)*cumsum(mt,2)*dt));
%%
sf1 = fftshift(abs(fft(st1,N)))/N;
sf2 = fftshift(abs(fft(st2,N)))/N;

figure (2)
subplot (211)
plot(f,sf1);
title("Magnitude of the Frequency Spectrum of Modulated Image with kf1");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("kf1=100");
subplot (212)
plot(f,sf2);
title("Magnitude of the Frequency Spectrum of Modulated Image with kf2");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("kf2=400");
%%
SNR = [0 25 50];
n0_1 = awgn(st1,SNR(1),"measured");
n0_2 = awgn(st2,SNR(1),"measured");
n25_1 = awgn(st1,SNR(2),"measured");
n25_2 = awgn(st2,SNR(2),"measured");
n50_1 = awgn(st1,SNR(3),"measured");
n50_2 = awgn(st2,SNR(3),"measured");
n0_1f = fftshift(abs(fft(n0_1,N)))/N;
n0_2f = fftshift(abs(fft(n0_2,N)))/N;

figure (3)
subplot (211)
plot(f,n0_1f);
title("Magnitude of the Frequency Spectrum of Modulated Noisy Image with kf1");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("R(f) with kf1");
subplot (212)
plot(f,n0_2f);
title("Magnitude of the Frequency Spectrum of Modulated Noisy Image with kf2");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("R(f) with kf2");
%% 8.3 Demodulation
%% 
Am = max(mt);
n25_1_demod = Am*fmdemod(n25_1,fc,Fs,Am*kf(1));
n25_2_demod = Am*fmdemod(n25_2,fc,Fs,Am*kf(2));

n25_1_demodf = fftshift(abs(fft(n25_1_demod,N)))/N;
n25_2_demodf = fftshift(abs(fft(n25_2_demod,N)))/N;

figure (4)
subplot (211)
plot(f,n25_1_demodf);
title("Magnitude of the Frequency Spectrum of Demodulated Noisy Image with kf1");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("kf1=100");
subplot (212)
plot(f,n25_2_demodf);
title("Magnitude of the Frequency Spectrum of Demodulated Noisy Image with kf2");
xlabel("frequency (Hz)");
ylabel("Magnitude (dB)");
legend("kf2=400");
%%
[b,a] = butter(5, 14500/(Fs/2));

n0_1_demod = Am*fmdemod(n0_1,fc,Fs,Am*kf(1));
n0_2_demod = Am*fmdemod(n0_2,fc,Fs,Am*kf(2));
n50_1_demod = Am*fmdemod(n50_1,fc,Fs,Am*kf(1));
n50_2_demod = Am*fmdemod(n50_2,fc,Fs,Am*kf(2));

n0_1_lpf = filter(b,a,n0_1_demod);
n0_2_lpf = filter(b,a,n0_2_demod);
n25_1_lpf = filter(b,a,n25_1_demod);
n25_2_lpf = filter(b,a,n25_2_demod);
n50_1_lpf = filter(b,a,n50_1_demod);
n50_2_lpf = filter(b,a,n50_2_demod);

PSNR0_1 = psnr(n0_1_lpf,mt);
PSNR0_2 = psnr(n0_2_lpf,mt);
PSNR25_1 = psnr(n25_1_lpf,mt);
PSNR25_2 = psnr(n25_2_lpf,mt);
PSNR50_1 = psnr(n50_1_lpf,mt);
PSNR50_2 = psnr(n50_2_lpf,mt);
%%
im0_1=reshape(n0_1_lpf,size(y,1),size(y,2));
im0_2=reshape(n0_2_lpf,size(y,1),size(y,2));
im25_1=reshape(n25_1_lpf,size(y,1),size(y,2));
im25_2=reshape(n25_2_lpf,size(y,1),size(y,2));
im50_1=reshape(n50_1_lpf,size(y,1),size(y,2));
im50_2=reshape(n50_2_lpf,size(y,1),size(y,2));

figure (5)
subplot(2,3,1)
imshow(im0_1);
title("0 dB AWGN image with kf1");
subplot(2,3,4)
imshow(im0_2);
title("0 dB AWGN image with kf2");
subplot(2,3,2)
imshow(im25_1);
title("25 dB AWGN image with kf1");
subplot(2,3,5)
imshow(im25_2);
title("25 dB AWGN image with kf2");
subplot(2,3,3)
imshow(im50_1);
title("50 dB AWGN image with kf1");
subplot(2,3,6)
imshow(im50_2);
title("50 dB AWGN image with kf2");
%% 8.4 Mean Square Error (MSE) and PSNR Comparison
%%
Pwithkf1 = [PSNR0_1 PSNR25_1 PSNR50_1];
Pwithkf2 = [PSNR0_2 PSNR25_2 PSNR50_2];
figure (6)
plot(SNR,Pwithkf1);
title("Power of the Noisy Demodulated & Filtered Images with kf1");
xlabel("SNR");
ylabel("POWER");
hold on
plot(SNR,Pwithkf2);
title("Power of the Noisy Demodulated & Filtered Images with kf2");
xlabel("SNR");
ylabel("POWER");
%%
MSE0_1 = immse(n0_1_lpf,mt);
MSE0_2 = immse(n0_2_lpf,mt);
MSE25_1 = immse(n25_1_lpf,mt);
MSE25_2 = immse(n25_2_lpf,mt);
MSE50_1 = immse(n50_1_lpf,mt);
MSE50_2 = immse(n50_2_lpf,mt);

MSEwithkf1 = [MSE0_1 MSE25_1 MSE50_1];
MSEwithkf2 = [MSE0_2 MSE25_2 MSE50_2];

figure (7)
plot(SNR,MSEwithkf1);
xlabel("SNR");
ylabel("MSE");
title("Mean Square Error (MSE) of the Noisy Demodulated & Filtered Images with kf1");
hold on
plot(SNR,MSEwithkf2);
title("Mean Square Error (MSE) of the Noisy Demodulated & Filtered Images with kf2");
xlabel("SNR");
ylabel("MSE");