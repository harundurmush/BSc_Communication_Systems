clear all;
close all;
clc;

[y,Fs] = audioread('handel.wav');
N = length(y);
mx = max(y)./2;
t = 0:1/Fs:((N-1)/Fs);
% sound(y,Fs);
Y = ['The sample rate of the signal is ',num2str(Fs)];   
disp(Y);

figure (1)
plot(t,y);
title('Input Audio Signal in Time domain');
xlabel('Time Axis');
ylabel('Amplitude of the Audio signal');

FFT = fft(y);

figure (2)
stem(abs(FFT));
title('Input Audio Signal in Frequency domain');
xlabel('Number of samples');
ylabel('Frequency Amplitude');
%% Encoding
 for i = 1 : N
     if(y(i) < mx)
         y(i) = 0;
     end
     
     if (y(i) >= mx)
         y(i) = 1;
     end
 end    
input = transpose(y);
txsig = 2*input-1;
M = 2;
%% Modulation
txSig = pskmod(y,M,pi/M);
%s = pskmod(input,M);
figure (3)
plot(t,txsig);
title('BPSK Modulated Signal in Time domain');
xlabel('Time Axis');
ylabel('Amplitude of the wave');
%% AWGN noise n generated
noise = 1/sqrt(2) * [randn(1,N) + j*randn(1,N)];
Eb_No = [-5:10]; % multiple Eb/N0 values

n = length(Eb_No);
%% r = x + n part
for i = 1 : n
   y = txsig + 10 ^ (-Eb_No(i) / 20) * noise;
   noisy = real(y) > 0;
   %BER Calculation
   BER(i) = size(find([input- noisy]),2);

end

BerT = 0.5 * erfc( sqrt(10 .^ (Eb_No / 10)) );
BerS = BER/N;

figure (4)
semilogy(Eb_No, BerT, 'b.-');
hold on
semilogy(Eb_No, BerS, 'mx-');
axis([-5 11 10^-5 0.5])
grid on
legend('Theoritical Values', 'Simulated Values');
xlabel('Eb/No in dB');
ylabel('BIT ERROR RATE');
title('BIT ERROR RATE for BPSK Modulated Audio Signal with AWGN');