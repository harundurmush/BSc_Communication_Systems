clc;
clear all;
close all;
%%
Fs = 1000;
Eb = (1/sqrt(2)) * sqrt(1/Fs);
sz = 250;
figure(1)
scatter(0, 0, sz, 'x', 'k', 'MarkerEdgeColor', 'black', 'LineWidth', 2); % Bold 'X'
hold on 
scatter(Eb, 0, sz, 'x', 'k', 'MarkerEdgeColor', 'black', 'LineWidth', 2); % Bold 'X'
grid on
xlabel('Ø_0(t)');
ylabel('Ø_1(t)');
title('Signal Constellation for B-ASK');
