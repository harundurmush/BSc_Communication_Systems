clc;
clear all;
%%
% Define the time vector
t = 0:0.001:4; % Time from 0 to 4 seconds with a sampling interval of 1 ms

% Initialize the signal as zeros
signal = zeros(size(t));

% Set the signal to 0 for t between 0 and 1 second
signal(t < 1) = 0;

% Set the signal to cos(2*pi*2000*t) for t between 1 and 3 seconds
signal = cos(2*pi*2000*t);

% Set the signal to 0 for t between 3 and 4 seconds
signal(t > 3) = 0;

% Plot the signal
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal: 0 to 1 second, cos(2*pi*2000*t) from 1 to 3 seconds, 0 from 3 to 4 seconds');
