clc;
clear all;
close all;
%% reading audio file
[x,fs] = audioread('handel.wav');
N = length(x);
% mx = max(x)/2; % threshold for encoding
% mx = 0.2;
%% quantization
b = max(x);
a = min(x);
Nq = 3; % quantization number
quantized = floor(((x-a)/(b-a))*(2^Nq-1))*((b-a)/(2^Nq-1)) + a;
mx = max(quantized)/2;
%% encoding (sanırım Tx Filter)
x_enc = [];
for i = 1:N
    if(quantized(i) < mx)
        x_enc(i) = 0;
    end 
    if (quantized(i) >= mx)
        x_enc(i) = 1;
    end
end    
encoded = transpose(x_enc);
p_encoded = 2*encoded-1;
M = 1; % modulation index

figure (1)
subplot(411)
plot(x);
subtitle("original audio signal");

subplot(412)
plot(quantized);
subtitle("quantized audio signal");

subplot(413)
plot(encoded);
subtitle("encoded audio signal");

subplot(414)
plot(p_encoded);
subtitle("modulated audio signal");
%% defining carrier signal
fc = floor(fs/10);
ts = 1/fs;
n_enc = length(encoded);
n = M*N; % bunun ne olduğunu anlamadım
tc = 1/fc;
t = 0:ts:n*ts-1;
carrier = sin(2*pi*fc*t);
%% converting impulse data into square data
% tp = 0:ts:tc*M;
% pulse_data_1 = [];
% for i=1:n_enc
%     for j=1:length(tp)/10-1
%         pulse_data_1 = [pulse_data p_encoded(i)];
%     end
% end
%% modulation
% pulse_data = transpose(pulse_data);
% modsig = pulse_data.*carrier;
modsig = p_encoded.*carrier;

figure(2)
grid on;
plot(carrier, 'g-');
hold on;
plot(modsig(1,:));
%%
demod = modsig(1,:).*carrier;
figure (3)
plot(demod);
