clc;
clear all;
close all;
%% BPSK (Binary Phase Shift Keying) on audio signal
%% Modulation
[audio,fs] = audioread('tired.wav');
audio = audio / max(abs(audio));
% sound(audio,fs);
fc = 1e6;
t = (0:1/fs:(length(audio)/fs)-1/fs)';
c_signal = cos(2*pi*fc*t);
pskmod_signal = audio .* c_signal;
rt = awgn(pskmod_signal,20);
% sound(rt,fs);
demodulated_signal = rt .* c_signal;
filtered_signal = lowpass(demodulated_signal, fc/2, fs);
decoded_audio = sign(filtered_signal);
sound(decoded_audio, fs);




