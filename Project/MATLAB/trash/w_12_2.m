clc;
clear all;
close all;
%%
Eb_N0_dB = -30:20; % multiple Eb/N0 values
nTAP = 4;
for i = 1:length(Eb_N0_dB)
    % bits = randi([0 1],1,1000000);
    % N_bits = length(bits);
    % Rb = 1000; % Bit rate in bits per second
    % tb = 1/Rb; % Bit time
    % fc = 10 * Rb; % Carrier frequency, for example, 10 times the bit rate
    % fs = fc; % Sampling frequency, 4 times the carrier frequency
    % ts = 1/fs;
    % k = 1; % binary
    % Rs = Rb/k;
    % sps = fs/Rs;
    % time = 0:ts:(N_bits*tb)-ts;
    bits = [0 1 0 0 1 1 0 1 0 1 1 0 0 0 1];
    N_bits = length(bits);
    tb = 0.1;
    fs = 300000;
    ts = 1/fs;
    %% pulse shaping for carrier multiplication
    %%
    bits_reshaped = reshape(bits, N_bits, 1);
    spb = tb*fs; % sample per bit
    message = repmat(bits_reshaped, 1, spb);
    message = reshape(message', 1, []);
    N_msg = length(message);
    % time = 0:ts:(N_msg*tb)-ts;
    %% modulation
    %%
    mxsig = 2*message-1;
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
    channel_response = [0.2 0.9 0.3]; % channel 1
    % channel_response = [0.1 0.2 0.1]; % channel 2
    % channel_response = [0.1 0.27 0.1]; % channel 3
    chan_out = conv(mxsig,channel_response);
    %% noise addition
    %%
    pavg_channel = sum(abs(mxsig).^2)/length(mxsig);
    snr_lin = 10^(0.1*Eb_N0_dB(i));
    var_noise = pavg_channel/snr_lin;
    noise = sqrt(var_noise)*randn(1,length(chan_out));
    noisy_out = chan_out + noise; % additive white gaussian noise
    %% equalization
    %%
    for k = 1:nTAP
        L  = length(channel_response);
        %% zero forcing equalizer
        %%
         channel_matrix = toeplitz([channel_response(2:end) zeros(1,2*k+1-L+1)], [ channel_response(2:-1:1) zeros(1,2*k+1-L+1) ]);
         d  = zeros(1,2*k+1);
         d(k+1) = 1;
         channel_zf  = (inv(channel_matrix)*d.').';
         %%  matched filter
         yFilt_zf = conv(noisy_out,channel_zf);
         yFilt_zf = yFilt_zf(k+2:end); 
         yFilt_zf = conv(yFilt_zf,ones(1,1)); % convolution
         ySamp_zf = yFilt_zf(1:1:N_msg);  % sampling at time T
         %% receiver - hard decision decoding
         ipHat_zf = real(ySamp_zf)>0;
         %% counting the errors
         nErr_zf(k,i) = size(find(message- ipHat_zf),2);
        %% mmse equalizer
        %%
         hAutoCorr = conv(channel_response,fliplr(channel_response));
         channel_matrix = toeplitz([hAutoCorr([3:end]) zeros(1,2*k+1-L)], [ hAutoCorr([3:end]) zeros(1,2*k+1-L) ]);
         channel_matrix = channel_matrix + 1/2*10^(-Eb_N0_dB(i)/10)*eye(2*k+1);
         d  = zeros(1,2*k+1);
         d([-1:1]+k+1) = fliplr(channel_response);
         channel_mmse  = [inv(channel_matrix)*d.'].';
         %%  matched filter
         yFilt_mmse = conv(noisy_out,channel_mmse);
         yFilt_mmse = yFilt_mmse(k+2:end); 
         yFilt_mmse = conv(yFilt_mmse,ones(1,1)); % convolution
         ySamp_mmse = yFilt_mmse(1:1:N_msg);  % sampling at time T
         %% receiver - hard decision decoding
         ipHat_mmse = real(ySamp_mmse)>0;
         %% counting the errors
         nErr_mmse(k,i) = size(find(message- ipHat_mmse),2);
     end
end
simBer_zf = nErr_zf/N_msg; % simulated ber
simBer_mmse = nErr_mmse/N_msg; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber
%% plotting
figure (1)
% semilogy(Eb_N0_dB,simBer_zf(1,:),'-.');
hold on
semilogy(Eb_N0_dB,simBer_zf(2,:),'-.','LineWidth',1.5);
% semilogy(Eb_N0_dB,simBer_zf(3,:),'-.');
% semilogy(Eb_N0_dB,simBer_zf(4,:),'-.');
% semilogy(Eb_N0_dB,simBer_mmse(1,:),'s-');
semilogy(Eb_N0_dB,simBer_mmse(2,:),'s-','LineWidth',1.5);
% semilogy(Eb_N0_dB,simBer_mmse(3,:),'s-');
% semilogy(Eb_N0_dB,simBer_mmse(4,:),'s-');
% axis([0 10 10^-3 0.5])
grid on
legend('sim-5tap zf','sim-5tap mmse');
% legend('sim-3tap zf', 'sim-5tap zf','sim-7tap zf','sim-9tap zf','sim-3tap mmse', 'sim-5tap mmse','sim-7tap mmse','sim-9tap mmse');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK in ISI with MMSE equalizer');
hold off

% eyediagram(simBer_zf(1,:),spb,ts,0);
% eyediagram(simBer_zf(2,:),spb,ts,0);
% eyediagram(simBer_zf(3,:),spb,ts,0);
% eyediagram(simBer_zf(4,:),spb,ts,0);
% eyediagram(simBer_mmse(1,:),spb,ts,0);
% eyediagram(simBer_mmse(2,:),spb,ts,0);
% eyediagram(simBer_mmse(3,:),spb,ts,0);
% eyediagram(simBer_mmse(4,:),spb,ts,0);
