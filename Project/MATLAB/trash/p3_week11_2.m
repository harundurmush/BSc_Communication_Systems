clc;
clear all;
close all;

bit_level = 8;
SNRinDb = 10;
%Load_the_wav_file
[y,Fs]=audioread('handel.wav');  
Length=length(y); 
sound(y,Fs);
t = [1/Fs:1/Fs:length(y)/Fs]; %kaç_saniye_sürüyor ?
mono_y = y(:,1); %stereo to mono

 


figure(1)
subplot(3,1,1);
plot(t,mono_y);
title('Original Signal');
xlabel('Time');
ylabel('Amplitude');

 

%Quantization
quant = (max(mono_y)+abs(min(mono_y)))/(2^(bit_level)); %calculating step interval. dividing interval between max value and min value by interval number.
y_quantized = round(mono_y/quant); %representing signal amplitude as stepsize.

 

 

%Plot quantized
subplot(3,1,2);
hold on
plot(t,mono_y);
stairs(t,y_quantized*(0.0068)); %multiplying quantized values with stepsize to scale quantized values with original signal.

 

title('Quantized Signal and Original Signal');
xlabel('Time');
ylabel('Amplitude');

 

 

%dec to binary
signe = floor((sign(y_quantized)+1)/2); %floor used for rounding 0.5 value to 0
b_out = [dec2bin(signe) dec2bin(abs(y_quantized),8)];

 

%bit_level+1 x Length matrix to 1 x 9*68867 array
%binary_array = zeros(1,(bit_level+1)*Length);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%matrix conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%long way
% for i = 1 : Length
%     for j = 1 : (bit_level+1)
%         binary_array(((i-1)*bit_level+1)+(j-1)) = b_out(i,j)-48; %b_out-48 converts char to number
%     end
% end

 

%short way
binary_array = reshape(b_out-48,[1,(bit_level+1)*Length]); %b_out-48 converts char to number

 


check = zeros(1,length(binary_array));
for i=1 : length(binary_array)
    if(binary_array(i) == 0)
        check(i) = (-1);
    else
        check(i) = 1;
    end
end

 


%BPSK modulation
nt = [1/Fs:1/Fs:length(binary_array)/Fs];
Fc = Fs*10000;
modulated_bpsk = check .* cos(2*pi*Fc*nt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR_arry = -10:50;%0 : 2 : 10;
SNR=10.^(SNR_arry/10);
error_arry = zeros(1,length(SNR_arry));
out_awgn_channel = awgn(modulated_bpsk, SNRinDb,'measured');
%plotting modulated_bpsk and out_awgn_channel
subplot(3,1,3);
hold on
plot(nt,modulated_bpsk);