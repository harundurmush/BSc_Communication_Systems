clc;
clear all;
close all;
%% zero forcing equalization
%%
Eb_N0_dB = -10:15; % multiple Eb/N0 values
nTAP = 4;
for i = 1:length(Eb_N0_dB)
    bits = randi([0 1],1,100000);
    N_bits = length(bits);
    R = 1000; % Bit rate in bits per second
    tb = 1/R; % Bit time
    fc = 10 * R; % Carrier frequency, for example, 10 times the bit rate
    fs = 4 * fc; % Sampling frequency, 4 times the carrier frequency
    ts = 1/fs;
    % time = 0:ts:(N_bits*tb)-ts;
    %% modulation
    %%
    mxsig = 2*bits-1;
    % %% pulse shaping for carrier multiplication
    % %%
    % bits_reshaped = reshape(bits, N_bits, 1);
    % spb = tb*fs; % sample per bit
    % message = repmat(bits_reshaped, 1, spb);
    % message = reshape(message', 1, []);
    % %% defining carrier
    % %%
    % carrier = sin(2*pi*fc*time);
    % %% multiplication
    % %% 
    % mxsig = message.*carrier;
    % %% figures
    % %%
    % figure(1)
    % subplot(311)
    % plot(time,message);
    % subplot(312)
    % plot(time,carrier);
    % subplot(313)
    % plot(time,mxsig);
    %% channel model (multipath channel)
    %%
    channel_response = [0.2 0.9 0.3]; 
    chan_out = conv(mxsig,channel_response);
    %% noise addition
    %%
    pavg_channel = sum(abs(chan_out).^2)/length(chan_out);
    snr_lin = 10^(0.1*Eb_N0_dB(i));
    var_noise = pavg_channel/snr_lin;
    noise = sqrt(var_noise)*randn(1,length(chan_out));
    noisy_out = chan_out + noise; % additive white gaussian noise
    %% zero forcing equalizer
    %%
    for k = 1:nTAP
         L  = length(channel_response);
         channel_matrix = toeplitz([channel_response(2:end) zeros(1,2*k+1-L+1)], [ channel_response(2:-1:1) zeros(1,2*k+1-L+1) ]);
         d  = zeros(1,2*k+1);
         d(k+1) = 1;
         channel  = (inv(channel_matrix)*d.').';
         %%  matched filter
         yFilt = conv(noisy_out,channel);
         yFilt = yFilt(k+2:end); 
         yFilt = conv(yFilt,ones(1,1)); % convolution
         ySamp = yFilt(1:1:N_bits);  % sampling at time T
         %% receiver - hard decision decoding
         ipHat = real(ySamp)>0;
         %% counting the errors
         nErr(k,i) = size(find(bits- ipHat),2);
     end
end
simBer = nErr/N_bits; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber
figure
semilogy(Eb_N0_dB,simBer(1,:),'bs-','Linewidth',1);
hold on
semilogy(Eb_N0_dB,simBer(2,:),'gd-','Linewidth',1);
semilogy(Eb_N0_dB,simBer(3,:),'ks-','Linewidth',1);
semilogy(Eb_N0_dB,simBer(4,:),'mx-','Linewidth',1);
% axis([0 10 10^-3 0.5])
grid on
legend('sim-3tap', 'sim-5tap','sim-7tap','sim-9tap');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK in ISI with ZF equalizer');