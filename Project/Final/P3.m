clc;
clear all;
close all;
%% signal generation
%%
N_bits = 100;
bit_duration = 1e-3;
samples_per_bit = 10;
bits = randi([0 1], 1, N_bits);
%% reshaping
%%
signal = repelem(bits, samples_per_bit);
%% creating time vector
%%
fs = 1 / (bit_duration/samples_per_bit);
ts = 1/fs;
time = 0:ts:(length(bits)*bit_duration) - ts;
%% figures
%%
figure;
subplot(211);
stem(bits, 'Marker', 'none','LineWidth',1.5);
title('Original Bits');
subplot(212);
plot(time, signal,'LineWidth',1.5);
title('Signal with 1 ms Bit Duration');
%% modulation (B-PSK)
%%
fc = fs*10; % Adjust the carrier frequency based on our requirements
mxsig = zeros(1, length(signal));
for i=1:length(signal)
    if signal(1,i) == 1
        mxsig(1,i) = cos(2*pi*fc*time(i)) ;
    elseif signal(1,i) == 0
        mxsig(1,i) = cos(2*pi*fc*time(i)+pi);
    end
end 
%% figures
%%
figure;
subplot(2,1,1);
stem(bits, 'Marker', 'none','LineWidth',1.5);
title('Original Bits');
subplot(2,1,2);
plot(time, mxsig,'LineWidth',1.5);
title('BPSK Modulated Signal');
scatterplot(mxsig, 1, 250, 'yx');
title('BPSK Constellation Diagram');
%% channel model
%%
channel_response = [0 0.2 0.3 0.9 0.3 0 0]; % channel 1. for channel 1, 7-tap filter is ideal
% channel_response = [0.2 0.9 0.3]; % channel 2
chan_out = conv(mxsig,channel_response);
%% noise addition
%%
Eb_N0_dB = 0:1:10;
for i = 1:length(Eb_N0_dB)
    pavg = sum(abs(mxsig).^2)/length(mxsig);
    snr_lin = 10^(0.1*Eb_N0_dB(i));
    var_noise = pavg/snr_lin;
    noise = sqrt(var_noise)*randn(1,length(chan_out));
    noisy_out = chan_out + noise; % additive white gaussian noise
    %% ZERO-FORCING EQUALIZATION
    %%
    L = length(channel_response);
    k = (L-1)/2; % (2*k+1)-tap equalizer design
    channel_matrix_zf = toeplitz([channel_response((L+1)/2:end) zeros(1,2*k+1-L+(L-1)/2)], [channel_response((L+1)/2:-1:1) zeros(1,2*k+1-L+(L-1)/2)]);
    d_zf = zeros(1,2*k+1);
    d_zf(k+1) = 1;
    channel_zf  = transpose(inv(channel_matrix_zf)*transpose(d_zf));
    %% matched filter
    N_msg = length(signal);
    yFilt_zf = conv(noisy_out,channel_zf);
    yFilt_zf = yFilt_zf(k+(L+1)/2:end); 
    yFilt_zf = conv(yFilt_zf,ones(1,1)); % convolution
    ySamp_zf = yFilt_zf(1:1:N_msg);  % sampling at time T
    %% receiver - hard decision decoding
    yDecoded_zf = real(ySamp_zf)>0;
    if(Eb_N0_dB(i) == 0)
        yout1_zf = ySamp_zf;
    elseif(Eb_N0_dB(i) == 10)
        yout2_zf = ySamp_zf;
    end
    %% counting the errors
    nErr_zf(i) = size(find(signal - yDecoded_zf),2);
    %% MMSE EQUALIZATION
    %%
    h_autocorr = conv(channel_response,fliplr(channel_response));
    channel_matrix_mmse = toeplitz([h_autocorr(L:end) zeros(1,2*k+1-L)], [h_autocorr(L:end) zeros(1,2*k+1-L)]);
    channel_matrix_mmse = channel_matrix_mmse + (1/2)*(10^(-Eb_N0_dB(i)/10))*eye(2*k+1);
    d_mmse = zeros(1,2*k+1);
    d_mmse([-(L-1)/2:(L-1)/2]+k+1) = fliplr(channel_response);
    channel_mmse  = transpose(inv(channel_matrix_mmse)*transpose(d_mmse));
    %%  matched filter
    yFilt_mmse = conv(noisy_out,channel_mmse);
    yFilt_mmse = yFilt_mmse(k+(L+1)/2:end); 
    yFilt_mmse = conv(yFilt_mmse,ones(1,1)); % convolution
    ySamp_mmse = yFilt_mmse(1:1:N_msg);  % sampling at time T
    %% receiver - hard decision decoding
    yDecoded_mmse = real(ySamp_mmse)>0;
    if(Eb_N0_dB(i) == 0)
        yout1_mmse = ySamp_mmse;
    elseif(Eb_N0_dB(i) == 10)
        yout2_mmse = ySamp_mmse;
    end
    %% counting the errors
    nErr_mmse(i) = size(find(signal - yDecoded_mmse),2);
end
simBer_zf = nErr_zf/N_msg; % simulated ber
simBer_mmse = nErr_mmse/N_msg; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(0.1*Eb_N0_dB))); % theoretical ber
%% figures
%%
figure;
semilogy(Eb_N0_dB,theoryBer,'-.','LineWidth',1.5);
hold on
semilogy(Eb_N0_dB,simBer_zf,'-.','LineWidth',1.5);
% hold on
semilogy(Eb_N0_dB,simBer_mmse,'-.','LineWidth',1.5);
hold off
grid on
legend('theoritical response','7-tap zf','7-tap mmse');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BIT ERROR PROBABILITY CURVE OVER DIFFERENT SNR VALUES');

eyediagram(yout1_zf, samples_per_bit, ts, 0);
title("Eye Diagram for SNR = 0 dB (zero-forcing)");

eyediagram(yout1_mmse, samples_per_bit, ts, 0);
title("Eye Diagram for SNR = 0 dB (MMSE)");

eyediagram(yout2_zf, samples_per_bit, ts, 0);
title("Eye Diagram for SNR = 10 dB (zero-forcing)");

eyediagram(yout2_mmse, samples_per_bit, ts, 0);
title("Eye Diagram for SNR = 10 dB (MMSE)");

scatterplot(yout1_zf, 1, 1, 'y.');
title("zero-forcing 0 dB");
scatterplot(yout1_mmse, 1, 1, 'y.');
title("MMSE 0 dB");
scatterplot(yout2_zf, 1, 1, 'y.');
title("zero-forcing 10 dB");
scatterplot(yout2_mmse, 1, 1, 'y.');
title("MMSE 10 dB");