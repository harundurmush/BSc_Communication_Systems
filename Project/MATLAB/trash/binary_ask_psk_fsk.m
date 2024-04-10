clc; clear all; close all;

[b,fs] = audioread('handel.wav');
% fs=40e3; 
fc=1000; 
tb=0.005;
t=0:1/fs:(length(b)*tb+(length(b)-1)/fs);
%% EXPRESSING THE GENERATED BIT STREAM
b1=[];
for i = 1:length(b)
    t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
    if (b(i)==1)
        b1 = [b1,ones(1,length(t1))];
    else 
        b1 = [b1,zeros(1,length(t1))];
    end
end 
% %% BASK MODULATION 
% ASK_modulated=[];
% for i = 1:length(b)
%     t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
%     if (b(i)==1)
%         ASK_modulated = [ASK_modulated,cos(2*pi*fc*t1)];
%     else 
%         ASK_modulated = [ASK_modulated,zeros(1,length(t1))];
%     end
% end
% % BASK DEMODULATION
% for i = 1:length(b)
%     x1=0;
%     x = ASK_modulated(((i-1)*tb*fs+i):(i*tb*fs+i));
%     t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
%     x1 = sum(x.*cos(2*pi*fc*t1));
%     if (x1>0)
%         demodulated_stream(i) = 1;
%     else
%         demodulated_stream(i) = 0;
%     end
% end
% ASK_demodulated=[];
% for i = 1:length(b)
%     t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
%     if (demodulated_stream(i)==1)
%         ASK_demodulated = [ASK_demodulated,ones(1,length(t1))];
%     else 
%         ASK_demodulated = [ASK_demodulated,zeros(1,length(t1))];
%     end
% end
%% BPSK MODULATION
PSK_modulated=[];
phase=pi;
for i = 1:length(b)
    t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
    if (b(i)==1)
        PSK_modulated = [PSK_modulated,cos(2*pi*fc*t1)];
    else 
        PSK_modulated = [PSK_modulated,(cos(2*pi*fc*t1+phase))];
    end
end
% BPSK DEMODULATION
for i = 1:length(b)
    x1=0;
    x = PSK_modulated(((i-1)*tb*fs+i):(i*tb*fs+i));
    t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
    x1 = sum(x.*cos(2*pi*fc*t1));
    if (x1>0)
        demodulated_stream(i) = 1;
    else
        demodulated_stream(i) = 0;
    end
end
PSK_demodulated=[];
for i = 1:length(b)
    t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
    if (demodulated_stream(i)==1)
        PSK_demodulated = [PSK_demodulated,ones(1,length(t1))];
    else 
        PSK_demodulated = [PSK_demodulated,zeros(1,length(t1))];
    end
end
% %% BFSK MODULATION 
% f1=2e3; 
% f2=1e3;
% FSK_modulated=[];
% for i = 1:length(b)
%     t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
%     if (b(i)==1)
%         FSK_modulated = [FSK_modulated,cos(2*pi*f1*t1)];
%     else 
%         FSK_modulated = [FSK_modulated,cos(2*pi*f2*t1)];
%     end
% end
% % BFSK DEMODULATION
% for i = 1:length(b)
%     x = FSK_modulated(((i-1)*tb*fs+i):(i*tb*fs+i));
%     t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
%     x1 = sum(x.*cos(2*pi*f1*t1));
%     x2 = sum(x.*cos(2*pi*f2*t1));
%     y=x1-x2;
%     if (y>0)
%         demodulated_stream(i) = 1;
%     else
%         demodulated_stream(i) = 0;
%     end
% end
% FSK_demodulated=[];
% for i = 1:length(demodulated_stream)
%     t1=t(((i-1)*tb*fs+i):(i*tb*fs+i));
%     if (demodulated_stream(i)==1)
%         FSK_demodulated = [FSK_demodulated,ones(1,length(t1))];
%     else 
%         FSK_demodulated = [FSK_demodulated,zeros(1,length(t1))];
%     end
% end
%% FIGURES
% figure;
% subplot(3,1,1);plot(t,b1);
% xlabel('t'); ylabel('amplitude'); title('Generated Bitstream');
% subplot(3,1,2);plot(t,ASK_modulated);
% xlabel('t'); ylabel('amplitude'); title('BASK Modulated Signal');
% subplot(3,1,3);plot(t,ASK_demodulated);
% xlabel('t'); ylabel('amplitude'); title('BASK Demodulated Bitstream');
figure;
subplot(3,1,1);plot(t,b1);;
xlabel('t'); ylabel('amplitude'); title('Generated Bitstream')
subplot(3,1,2);plot(t,PSK_modulated);
xlabel('t'); ylabel('amplitude'); title('BPSK Modulated Signal');
subplot(3,1,3);plot(t,PSK_demodulated); 
xlabel('t'); ylabel('amplitude'); title('BPSK Demodulated Bitstream');
% figure;
% subplot(3,1,1);plot(t,b1);
% xlabel('t'); ylabel('amplitude'); title('Generated Bitstream');
% subplot(3,1,2);plot(t,FSK_modulated); 
% xlabel('t'); ylabel('amplitude'); title('BFSK Modulated Signal'); 
% subplot(3,1,3);plot(t,FSK_demodulated);
% xlabel('t'); ylabel('amplitude'); title('BFSK Demodulated Bitstream');