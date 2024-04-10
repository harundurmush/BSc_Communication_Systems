clc;
clear all;
close all;
%% reading audio file
[x,fs] = audioread('handel.wav');
N = length(x);
% mx = max(x)/2; % threshold for encoding
% mx = 0.2;
%% quantization
b = max(audio);
a = min(audio);
%% encoding (sanırım Tx Filter)
% x_enc = [];
% for i = 1:N
%     if(x(i) < mx)
%         x_enc(i) = 0;
%     end 
%     if (x(i) >= mx)
%         x_enc(i) = 1;
%     end
% end    
% encoded = transpose(x_enc);
% p_encoded = 2*encoded-1;
% M = 0.1; % modulation index

figure (1)
subplot(311)
plot(x);
subtitle("original audio signal");

subplot(312)
plot(encoded);
subtitle("encoded audio signal");

subplot(313)
plot(p_encoded);
subtitle("modulated audio signal");
%% defining carrier signal
fc = floor(fs/10);
ts = 1/fs;
n_enc = length(encoded);
n = M*N; % bunun ne olduğunu anlamadım
tc = 1/fc;
t = 0:ts:n*tc;
carrier = sin(2*pi*fc*t);
%% converting impulse data into square data
tp = 0:ts:tc*M;
pulse_data = [];
for i=1:n_enc
    for j=1:length(tp)-1
        pulse_data = [pulse_data p_encoded(i)];
    end
end
%% modulation
pulse_data = transpose(pulse_data);
modsig = pulse_data.*carrier;

figure(2)
% plot(pulse_data, 'r-', 'LineWidth', 4);
grid on;
hold on;
plot(carrier, 'g-');
hold on;
% plot(modsig(1,:), 'b-', 'LineWidth', 2);

