clc;
clear all;
close all;
%% 5.1      B-ASK Modulation and Demodulation
%% 5.1.a    (defining parameters and the data to be transmitted)
%%
bits = [0 1 0 1 1];
N_bits = length(bits);
tb = 0.02;
fs = 10000;
ts = 1/fs;
time = 0:ts:(N_bits*tb)-ts;
%% 5.1.b    (reshaping)
%%
bits_reshaped = reshape(bits, N_bits, 1);
spb = tb*fs; % sample per bit
message = repmat(bits_reshaped, 1, spb);
message = reshape(message', 1, []);
%% 5.1.c    (modulation)
%%
fc = 2500;
carrier_bask = cos(2*pi*fc*time);
mod_bask = 5 * message .* carrier_bask;
% modulation for bit 0 equals to 0, since it multiplies with 0.
%% 5.1.d    (noise addition)
%%
pavg_bask = sum(abs(mod_bask).^2)/length(mod_bask);
snr_db = 10;
snr_lin = 10^(0.1*snr_db);
var_noise_bask = pavg_bask/snr_lin;
awgn_bask = sqrt(var_noise_bask)*randn(1,length(mod_bask));
noisy_mod_bask = mod_bask + awgn_bask;
%% 5.1.e    (demodulation)
%%
demod_bask = zeros(size(message));
max_lag = 0; % defined maximum lag
for i=1:N_bits
    segment_bask = noisy_mod_bask((i-1)*spb + (1:spb));
    correlation_bask = xcorr(segment_bask, carrier_bask(1:spb), max_lag);
    demod_bask(i) = sum(abs(correlation_bask));
end
threshold_bask = max(demod_bask)/2;
demodulation_bask = zeros(1,N_bits);
for i = 1:N_bits
    if demod_bask(i) > threshold_bask 
        demodulation_bask(i) = 1; 
    else 
        demodulation_bask(i) = 0; 
    end
end
demodulation_bask = reshape(demodulation_bask, length(demodulation_bask), 1);
demodulation_bask = repmat(demodulation_bask, 1, spb);
demodulation_bask = reshape(demodulation_bask', 1, []);
%% 5.1.f    (results)
%%
figure(1)
subplot(411);
plot(time, message,"color","#0072BD","LineWidth",1.5);
title("MESSAGE SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("message signal [0 1 0 1 1]");
grid on;
subplot(412);
plot(time, mod_bask,"color","#A2142F");
title("B-ASK MODULATED SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("b-ask modulated signal");
grid on;
subplot(413);
plot(time, noisy_mod_bask,"color","#D95319");
title("NOISY B-ASK MODULATED SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("noisy b-ask modulated signal");
grid on;
subplot(414);
plot(time, demodulation_bask,"color","#EDB120","LineWidth",1.5);
title("B-ASK DEMODULATED SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("b-ask demodulated signal");
grid on;
%% 5.2      B-FSK Modulation and Demodulation
%% 5.2.a    (carrier and modulation)
%%
f0=500;
f1=750;
carrier_fsk_0=cos(2*pi*f0*time);
carrier_fsk_1=cos(2*pi*f1*time);
mod_fsk = zeros(1, length(message));
for i=1:N_bits
    start_index = (i-1)*spb + 1;
    end_index = i*spb;
    if bits(i) == 0
        mod_fsk(start_index:end_index) = cos(2*pi*f0*time(start_index:end_index));
    else
        mod_fsk(start_index:end_index) = cos(2*pi*f1*time(start_index:end_index));
    end
end
%% 5.2.d    (noise addition)
%%
pavg_fsk = sum(abs(mod_fsk).^2)/length(mod_fsk);
var_noise_fsk = pavg_fsk/snr_lin;
awgn_fsk = sqrt(var_noise_fsk)*randn(1,length(mod_fsk));
noisy_mod_fsk = mod_fsk + awgn_fsk;
%% 5.2.e    (demodulation and detection)
%%
demod_fsk = zeros(size(message));
for i = 1:N_bits
    segment_fsk = noisy_mod_fsk((i-1)*spb + (1:spb));
    correlation_fsk_0 = xcorr(segment_fsk, carrier_fsk_0(1:spb), max_lag);
    correlation_fsk_1 = xcorr(segment_fsk, carrier_fsk_1(1:spb), max_lag);
    demodulation_fsk_0(i) = sum(abs(correlation_fsk_0));
    demodulation_fsk_1(i) = sum(abs(correlation_fsk_1));
end
demodulation_fsk_2 = demodulation_fsk_1 - demodulation_fsk_0;
threshold_fsk = max(demodulation_fsk_2)/2;
demodulation_fsk = zeros(1,N_bits);
for i = 1:N_bits
    if demodulation_fsk_2(i) > threshold_fsk 
        demodulation_fsk(i) = 1; 
    else 
        demodulation_fsk(i) = 0; 
    end
end
demodulation_fsk = reshape(demodulation_fsk, length(demodulation_fsk), 1);
demodulation_fsk = repmat(demodulation_fsk, 1, spb);
demodulation_fsk = reshape(demodulation_fsk', 1, []);
%% 5.2.f    (results)
%%
figure(2)
subplot(411);
plot(time, message,"color","#0072BD","LineWidth",1.5);
title("MESSAGE SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("message signal [0 1 0 1 1]");
grid on;
subplot(412);
plot(time, mod_fsk,"color","#A2142F");
title("B-FSK MODULATED SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("b-fsk modulated signal");
grid on;
subplot(413);
plot(time, noisy_mod_fsk,"color","#D95319");
title("NOISY B-FSK MODULATED SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("noisy b-fsk modulated signal");
grid on;
subplot(414);
plot(time, demodulation_fsk,"color","#EDB120","LineWidth",1.5);
title("B-FSK DEMODULATED SIGNAL");
xlabel("message duration (s)");
ylabel("amplitude (V)");
legend("b-fsk demodulated signal");
grid on;