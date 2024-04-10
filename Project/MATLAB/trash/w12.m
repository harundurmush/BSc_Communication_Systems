clc;
clear all;
close all;
[signal, Fs] = audioread('handel.wav'); % Fs is the sampling frequency
% Assuming signal is mono-channel
signal = signal(:,1);

% Normalize the signal
signal = (signal - min(signal)) / (max(signal) - min(signal));

% Convert to binary
binarySignal = de2bi(round(signal * 7)); % 7 because we have 8 levels (0 to 7)

% Reshape to groups of 3 bits
numSymbols = floor(length(binarySignal)/3);
binarySymbols = binarySignal(1:numSymbols*3);
binarySymbols = reshape(binarySymbols, [], 3);

theta = 0:pi/4:2*pi-pi/4; % 8 phase shifts for 8-PSK
modulatedSignal = zeros(size(binarySymbols, 1), 1);

for i = 1:size(binarySymbols, 1)
    symbol = binarySymbols(i, :);
    index = bi2de(symbol, 'left-msb') + 1;
    modulatedSignal(i) = exp(1j*theta(index));
end

demodulatedSymbols = zeros(size(binarySymbols));

for i = 1:length(modulatedSignal)
    [~, index] = min(abs(modulatedSignal(i) - exp(1j*theta)));
    demodulatedSymbols(i, :) = de2bi(index-1, 3, 'left-msb');
end

% Flatten the array back to binary
demodulatedBinary = reshape(demodulatedSymbols', 1, []);

% Truncate or pad the demodulated binary signal
if length(demodulatedBinary) > length(binarySignal)
    demodulatedBinary = demodulatedBinary(1:length(binarySignal));
elseif length(demodulatedBinary) < length(binarySignal)
    % Pad with zeros if necessary
    demodulatedBinary = [demodulatedBinary zeros(1, length(binarySignal) - length(demodulatedBinary))];
end

% Convert back to decimal
demodulatedSignal = bi2de(demodulatedBinary, 'left-msb') / 7;

% Scale back to original signal range
demodulatedSignal = demodulatedSignal * (max(signal) - min(signal)) + min(signal);

% Ensure the demodulated signal has the same number of samples as the original
if length(demodulatedSignal) > length(signal)
    demodulatedSignal = demodulatedSignal(1:length(signal));
elseif length(demodulatedSignal) < length(signal)
    % Pad with the mean value (or another suitable value) if necessary
    demodulatedSignal = [demodulatedSignal; ones(length(signal) - length(demodulatedSignal), 1) * mean(demodulatedSignal)];
end

% Plotting the original and recovered signals
figure
subplot(2,1,1)
plot(signal)
title('Original Signal')

subplot(2,1,2)
plot(demodulatedSignal)
title('Recovered Signal')



