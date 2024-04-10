clc;
clear all;
close all;
%% PPM (Pulse Position Modulation) on audio signal
%% Reading Audio File
[audio,fs] = audioread('tired.wav');
audio = audio / max(abs(audio));
threshold = 0.5;
binary_sequence = (audio > threshold);
%% Applying PPM
time_interval = 0.01;
symbol_duration = 0.005;

ppm_signal = zeros(1, round(length(binary_sequence) * symbol_duration * fs));
symbol_index = 1;

for i = 1:length(binary_sequence)
    if binary_sequence(i) == 1
        pulse_position = round((symbol_index - 1) * symbol_duration * fs + 1);
        ppm_signal(pulse_position : pulse_position + round(symbol_duration * fs) - 1) = 1;
    end

    symbol_index = symbol_index + 1;
end
%% Demodulation
Pavg= sum(abs(ppm_signal).^2)/length(ppm_signal);
snr_db=[30 15 0];
snr_lin(1)=10^(0.1*snr_db(1));
snr_lin(2)=10^(0.1*snr_db(2));
snr_lin(3)=10^(0.1*snr_db(3));
var(1)=Pavg/snr_lin(1);
var(2)=Pavg/snr_lin(2);
var(3)=Pavg/snr_lin(3);
awgn30=sqrt(var(1))*randn(1,length(ppm_signal));
awgn15=sqrt(var(2))*randn(1,length(ppm_signal));
awgn0=sqrt(var(3))*randn(1,length(ppm_signal));
rt30 = ppm_signal + awgn30;
rt15 = ppm_signal + awgn15;
rt0 = ppm_signal + awgn0;
channel_coefficients = [1, 0.5, 0.2];  % Example channel coefficients
%% Zero Forcing Equalizer
equalizer_taps = fliplr(pinv(toeplitz(channel_coefficients)));
rt30 = conv(transpose(rt30), channel_coefficients);
rt30 = conv(rt30, equalizer_taps(:));

rt15 = conv(transpose(rt15), channel_coefficients);
rt15 = conv(rt15, equalizer_taps(:));

rt0 = conv(transpose(rt0), channel_coefficients);
rt0 = conv(rt0, equalizer_taps(:));

demodulated_signal30 = zeros(1, length(rt30));
demodulated_signal15 = zeros(1, length(rt15));
demodulated_signal0 = zeros(1, length(rt0));

for i = 1:length(demodulated_signal30)
    if rt30(i) > threshold
        demodulated_signal30(i) = 1;
    else
        demodulated_signal30(i) = 0;
    end
end

for i = 1:length(demodulated_signal15)
    if rt15(i) > threshold
        demodulated_signal15(i) = 1;
    else
        demodulated_signal15(i) = 0;
    end
end

for i = 1:length(demodulated_signal0)
    if rt0(i) > threshold
        demodulated_signal0(i) = 1;
    else
        demodulated_signal0(i) = 0;
    end
end
figure(1)
subplot(311)
plot(ppm_signal);
hold on
plot(demodulated_signal30);
title("PPM signal vs. PPM Demodulated signal SNR:30dB");
legend("ppm","demodulated snr=30");
xlabel("sample time (Ts)");
ylabel("bits");

subplot(312)
plot(ppm_signal);
hold on
plot(demodulated_signal15);
title("PPM signal vs. PPM Demodulated signal SNR:15dB");
legend("ppm","demodulated snr=15");
xlabel("sample time (Ts)");
ylabel("bits");

subplot(313)
plot(ppm_signal);
hold on
plot(demodulated_signal0);
title("PPM signal vs. PPM Demodulated signal SNR:0dB");
legend("ppm","demodulated snr=0");
xlabel("sample time (Ts)");
ylabel("bits");
%% MMSE Equalizer
mmse_equalizer = inv(conv(channel_coefficients, conj(flip(channel_coefficients))) + var(3));
equalized_signal0 = zeros(size(rt0));
equalized_signal0 = conv(rt0, mmse_equalizer, 'full');
demodulated_signal0 = zeros(size(equalized_signal0));
demodulated_signal0 = equalized_signal0(i, :) > threshold;
figure(1)
plot(demodulated_signal0);