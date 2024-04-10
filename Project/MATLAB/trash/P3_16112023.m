clc;
clear all;
close all;
%%
% Parameters
numSamples = 1000;  % Number of samples in the signal
SNR_dB = 20;       % Signal-to-Noise Ratio in decibels
t = linspace(0, 1, numSamples);
signal = sin(2 * pi * 5 * t);  % Example signal (5 Hz sinusoidal)
SNR_linear = 10^(SNR_dB / 10);
noisePower = var(signal) / SNR_linear;
noise = sqrt(noisePower) * randn(1, numSamples);
receivedSignal = signal + noise;
noiseCoefficients = sqrt(var(signal) / 10^(SNR_dB / 10)) * randn(1, numSamples);
figure;
plot(t, signal, 'LineWidth', 2, 'DisplayName', 'Transmitted Signal');
hold on;
plot(t, receivedSignal, 'LineWidth', 2, 'DisplayName', 'Received Signal with AWGN');
title('AWGN Channel Simulation');
xlabel('Time');
ylabel('Amplitude');
legend('show');
grid on;
hold off;





