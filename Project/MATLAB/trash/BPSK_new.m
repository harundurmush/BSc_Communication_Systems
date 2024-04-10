clc;
clear all;
close all;

% Generate 100 bits with a bit duration of 1 ms and 10 samples per bit:
N_bits = 1000;
bit_duration = 1;
samples_per_bit = 10;
bits = randi([0 1], 1, N_bits);

% Create the signal by repeating each bit 10 times:
signal = repelem(bits, samples_per_bit);

% Time vector and sampling frequency:
fs = 1 / (bit_duration / samples_per_bit);
time = 0:1/fs:(length(bits) * bit_duration) - 1/fs;

figure;
subplot(2,1,1);
stem(bits, 'Marker', 'none');
title('Original Bits');
subplot(2,1,2);
plot(time, signal);
title('Signal with 1 ms Bit Duration');

% BPSK modulation
fc = 100*fs/samples_per_bit; % Adjust the carrier frequency based on your requirements
mxsig = zeros(1, length(signal));
for i =1: length(signal)
    if signal( 1 , i ) == 0
        mxsig( 1 , i ) = cos(2*pi*fc*time(i)) ;
    else
        mxsig( 1 , i ) = sin(2*pi*fc*time(i)) ;
    end
end 

figure;
subplot(2,1,1);
stem(bits, 'Marker', 'none');
title('Original Bits');
subplot(2,1,2);
plot(time, mxsig);
title('BPSK Modulated Signal');

scatterplot(mxsig, 1, 0, 'yx');
title('BPSK Constellation Diagram');
