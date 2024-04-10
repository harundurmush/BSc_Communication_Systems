clc;
clear all;
close all;
%% signal generation
%%
N_bits = 10000;
bit_duration = 1e-3;
samples_per_bit = 10;
bits = randi([0 1], 1, N_bits);
%% reshaping
%%
signal = repelem(bits, samples_per_bit);
%% creating time vector
%%
fs = 1 / (bit_duration/samples_per_bit);
time = 0:1/fs:(length(bits)*bit_duration) - 1/fs;
%% figures
%%
figure;
subplot(2,1,1);
stem(bits, 'Marker', 'none');
title('Original Bits');
subplot(2,1,2);
plot(time, signal);
title('Signal with 1 ms Bit Duration');
%% modulation (B-PSK)
%%
fc = 10000; % Adjust the carrier frequency based on your requirements
mxsig = zeros(1, length(signal));
for i=1:length(signal)
    if signal(1,i) == 0
        mxsig(1,i) = cos(2*pi*fc*time(i)) ;
    else
        mxsig(1,i) = cos(2*pi*fc*time(i)+pi);
    end
end 
%% figures
%%
figure;
subplot(2,1,1);
stem(bits, 'Marker', 'none');
title('Original Bits');
subplot(2,1,2);
plot(time, mxsig);
title('BPSK Modulated Signal');
scatterplot(mxsig, 1, 0, 'b.');
title('BPSK Constellation Diagram');
%% channel model
%%
channel_response = [0.05 0.1 0.2 0.9 0.3 0.1 0.05]; % channel 1. for channel 1, 5-tap filter is ideal
% channel_response = [0.2 0.9 0.3]; % channel 1
% channel_response = [0.1 0.2 0.1]; % channel 2
% channel_response = [0.1 0.16 0.27 0.15 0.1]; % channel 3
chan_out = conv(mxsig,channel_response);
Eb_N0_dB = 10;
%% noise addition
%%
pavg = sum(abs(mxsig).^2)/length(mxsig);
snr_lin = 10^(0.1*Eb_N0_dB);
var_noise = pavg/snr_lin;
noise = sqrt(var_noise)*randn(1,length(chan_out));
noisy_out = chan_out + noise; % additive white gaussian noise
%% zero-forcing
%%
z_channel_zf = ztrans(channel_response);
















