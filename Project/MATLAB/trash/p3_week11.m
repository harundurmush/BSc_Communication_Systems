clc;
clear all;
close all;
%% RUN ETMESİ UZUN SÜRÜYOR. mod_psk 1x1843200 uzunluğunda olduğu için (line 142)
%% reading audio file
[x,fs] = audioread('handel.wav');
N = length(x);
%% quantization
b = max(x);
a = min(x);
Nq = 3; % quantization number
quantized = floor(((x-a)/(b-a))*(2^Nq-1))*((b-a)/(2^Nq-1)) + a;
mx = max(quantized)/2;

figure (1)
subplot(211)
plot(x);
subplot(212)
plot(quantized);
% sound(quantized,fs);
% sound(x,fs);
%% huffman source-coding theorem
inputs = unique(quantized);
occurances = zeros(1,8);
for i=1:length(inputs)
    for j=1:length(quantized)
        if(inputs(i)==quantized(j))
            occurances(i) = occurances(i) +1;
        end
    end
end
probabilities = transpose(occurances./N);
% Input -0.800018310546875: Huffman Code - 0101001
% Input -0.571446010044643: Huffman Code - 01011
% Input -0.342873709542411: Huffman Code - 00
% Input -0.114301409040179: Huffman Code - 1
% Input 0.114270891462054:  Huffman Code - 011
% Input 0.342843191964286:  Huffman Code - 0100
% Input 0.571415492466518:  Huffman Code - 010101
% Input 0.799987792968750:  Huffman Code - 0101000
%% açıklama
% Huffman codelarını chat gpt ile hesaplattık. probabilityleri bulduk
% leveller belli. bunları gptye verdik hesapla dedik.
%% encoding audio signal
encoded = [];
for j=1:length(quantized)
    if(quantized(j)==inputs(1))
        encoded = [encoded 0 1 0 1 0 0 1];
    end
    if(quantized(j)==inputs(2))
        encoded = [encoded 0 1 0 1 1];
    end
    if(quantized(j)==inputs(3))
        encoded = [encoded 0 0];
    end
    if(quantized(j)==inputs(4))
        encoded = [encoded 1];
    end
    if(quantized(j)==inputs(5))
        encoded = [encoded 0 1 1];
    end
    if(quantized(j)==inputs(6))
        encoded = [encoded 0 1 0 0];
    end
    if(quantized(j)==inputs(7))
        encoded = [encoded 0 1 0 1 0 1];
    end
    if(quantized(j)==inputs(8))
        encoded = [encoded 0 1 0 1 0 0 0];
    end
end
figure(2)
plot(encoded);
%% defining carrier
fc = fs*2; % fc fs'ten büyük olması lazım ama mod_psk çok büyük olmadığı için x2 yaptık sadece
t = linspace(0,2*pi);
c1 = cos(2*pi*fc*t);
c2 = cos(2*pi*fc*t + pi/4);
c3 = cos(2*pi*fc*t + 2*pi/4);
c4 = cos(2*pi*fc*t + 3*pi/4);
c5 = cos(2*pi*fc*t + 4*pi/4);
c6 = cos(2*pi*fc*t + 5*pi/4);
c7 = cos(2*pi*fc*t + 6*pi/4);
c8 = cos(2*pi*fc*t + 7*pi/4);

% carrier = [c1, c2, c3, c4, c5, c6, c7, c8];

figure (3)
subplot(421)
plot(c1);
grid on;
subtitle("carrier c1(t)");

subplot(422)
plot(c2);
grid on;
subtitle("carrier c2(t)");

subplot(423)
plot(c3);
grid on;
subtitle("carrier c3(t)");

subplot(424)
plot(c4);
grid on;
subtitle("carrier c4(t)");

subplot(425)
plot(c5);
grid on;
subtitle("carrier c5(t)");

subplot(426)
grid on;
plot(c6);
subtitle("carrier c6(t)");

subplot(427)
plot(c7);
grid on;
subtitle("carrier c7(t)");

subplot(428)
plot(c8);
grid on;
subtitle("carrier c8(t)");
%% modulation
% Input -0.800018310546875: Huffman Code - 0101001
% Input -0.571446010044643: Huffman Code - 01011
% Input -0.342873709542411: Huffman Code - 00
% Input -0.114301409040179: Huffman Code - 1
% Input 0.114270891462054:  Huffman Code - 011
% Input 0.342843191964286:  Huffman Code - 0100
% Input 0.571415492466518:  Huffman Code - 010101
% Input 0.799987792968750:  Huffman Code - 0101000
%% açıklama
% aşağıda yaptığımız şey, her bir encoded sequence (huffman codeları) aynı uzunlukta olmadığı için quantized sinyalde ilgili 
% leveli gördüğü zaman ilgili sequence ile ilgili carrierı çarparak modüle
% etmek ve uç uca eklemek, aslında serial to parallel conversion yapıp
% sonra carrierlarla çarpıp tekrar birleştirmiş olduk
%%
i=1;
mod_psk = [];
for k=1:N
    if(quantized(k)==inputs(1))
        temp_modulation = transpose(encoded(i:i+6)).*c1;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 7;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(2))
        temp_modulation = transpose(encoded(i:i+4)).*c2;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 5;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(3))
        temp_modulation = transpose(encoded(i:i+1)).*c3;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 2;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(4))
        temp_modulation = transpose(encoded(i)).*c4;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 1;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(5))
        temp_modulation = transpose(encoded(i:i+2)).*c5;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 3;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(6))
        temp_modulation = transpose(encoded(i:i+3)).*c6;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 4;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(7))
        temp_modulation = transpose(encoded(i:i+5)).*c7;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 6;
        temp_modulation = 0;
    end
    if(quantized(k)==inputs(8))
        temp_modulation = transpose(encoded(i:i+6)).*c8;
        mod_psk = [mod_psk temp_modulation(1,:)];
        i = i + 6;
        temp_modulation = 0;
    end
end
figure (4)
plot(mod_psk);

