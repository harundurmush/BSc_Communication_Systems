clc;
clear all;
close all;
%% 7.1 Image Reading
%% 
I = imread('cameraman.tif');
y = im2double(I);
%%
Fs = size(y,1)*size(y,2);
%% 7.2 Modulation
%% 
y_imgvec = reshape(y,1,[]);
t = 0:(1/Fs):(numel(y_imgvec)-1)/Fs;
%%
fc=20000;
ct = cos(2*pi*fc*t);
%%
ka=0.9;
m = (1+ka.*y_imgvec).*ct;
%%
SNR = [0 10 30];
%%
P = sum(abs(m).^2) / length(m);
%%
SNRlin0 = 10.^(0.1*SNR(1));
SNRlin10 = 10.^(0.1*SNR(2));
SNRlin30 = 10.^(0.1*SNR(3));
%%
var0 = P/SNRlin0;
var10 = P/SNRlin10;
var30 = P/SNRlin30;
var = [ var0 var10 var30];
%%
n0 = sqrt(var(1))*randn(1,length(m));
n10 = sqrt(var(2))*randn(1,length(m));
n30 = sqrt(var(3))*randn(1,length(m));
%%
r0 = m + n0;
r10 = m + n10;
r30 = m + n30;
%% 7.3 Demodulation
%% 
term0 = r0.^2;
term10 = r10.^2;
term30 = r30.^2;
[b, a] = butter(2, fc/2/(Fs/2));
filtered0 = sqrt(filter(b,a,term0));
filtered10 = sqrt(filter(b,a,term10));
filtered30 = sqrt(filter(b,a,term30));
rprime0=(filtered0-max(ct))/ka;
rprime10=(filtered10-max(ct))/ka;
rprime30=(filtered30-max(ct))/ka;
%%
reshaped0 = reshape(rprime0,256,256);
reshaped10 = reshape(rprime10,256,256);
reshaped30 = reshape(rprime30,256,256);
%%
figure (1)
subplot(221);
imshow(y);
title("Original Image");
subplot(222);
imshow(reshaped0);
title("Demodulated Image with SNR 0 dB");
subplot(223);
imshow(reshaped10);
title("Demodulated Image with SNR 10 dB");
subplot(224);
imshow(reshaped30);
title("Demodulated Image with SNR 30 dB");
%% 7.4 Mean Square Error (MSE) and comparison
%%
mse0 = sum(sum((y-reshaped0).^2),2)/numel(y);
mse10 = sum(sum((y-reshaped10).^2),2)/numel(y);
mse30 = sum(sum((y-reshaped30).^2),2)/numel(y);
MSE = [mse0 mse10 mse30];
figure (2)
plot(SNR,MSE,"color","#A2142F");
title("Mean Square Error (MSE) / SNR");
ylabel('MSE');
xlabel('SNR');
%%