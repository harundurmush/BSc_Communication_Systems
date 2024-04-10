clear all;
clc;
axis = -5:5;
channel_response1 = [0 0 0.05 0.1 0.2 0.9 0.3 0.1 0.05 0 0]; % channel 1. for channel 1, 7-tap filter is ideal
channel_response2 = [0 0 0 0 0.2 0.9 0.3 0 0 0 0]; % channel 2
figure;
subplot(211)
stem(axis, channel_response1,'LineWidth', 1.5);
xlabel("sampled time (T)");
ylabel("amplitude");
title("CHANNEL 1");
subplot(212)
stem(axis, channel_response2,'LineWidth', 1.5);
xlabel("sampled time (T)");
ylabel("amplitude");
title("CHANNEL 2");