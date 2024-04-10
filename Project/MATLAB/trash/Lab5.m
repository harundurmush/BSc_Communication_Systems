clear all; clc; close all;

%% 5.1 Modulation 
%a
[b,Fs] = audioread('handel.wav');
% Fs = 45000;
Ts = 1/Fs;
fc = 3000;
Tb = 0.002;
bitSample = Tb/Ts;
% b = randi([0 1],10,1);
% b2 = repmat(b,1,bitSample)';
% b2 = b2(:)';

duration = length(b)*Tb;
time = 0:Ts:(duration-Ts);
tim = 0:Ts:(Tb-Ts);

bpsk_mod = [];
for i = 1:length(b)
    if (b(i) == 0) 
        bpsk_mod = [bpsk_mod cos(2*pi*fc*(((i-1)*Tb)+tim))];
    else 
        bpsk_mod = [bpsk_mod cos(2*pi*fc*(((i-1)*Tb)+tim)+pi)];
    end
end

bpsk_dmd = [];
for i = 1:length(b)
    L0 = sum(bpsk_mod((1:bitSample)+bitSample*(i-1)).*(cos(2*pi*fc*(((i-1)*Tb)+tim))));
    L1 = sum(bpsk_mod((1:bitSample)+bitSample*(i-1)).*(cos(2*pi*fc*(((i-1)*Tb)+tim)+pi)));
    L = L1-L0;
    if (L > 0)
        bpsk_dmd = [bpsk_dmd 1];
    else 
        bpsk_dmd = [bpsk_dmd 0];
    end
end

bpsk_dmd = repmat(bpsk_dmd',1,bitSample)';
bpsk_dmd = bpsk_dmd(:)';

figure(3) 
subplot(311);
plot(time, b2);
title('Bit Stream');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(312);
plot(time,bpsk_mod);
title('Binary PSK Modulation');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(313);
plot(time,bpsk_dmd);
title('Demodulated Bit Stream');
xlabel('Time (s)'); ylabel('Amplitude');
