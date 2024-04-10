% Generate random bits
    bits = randi([0 3], 1, 100);
    N_bits = length(bits);

    % Modulation parameters
    fs = 200000; % Sampling frequency
    ts = 1/fs;   % Sampling time
    fc = 1000;   % Carrier frequency for BPSK

    % Time vector
    time = 0:ts:(N_bits*ts)-ts;

    % BPSK modulation using pskmod
    M = 4; % BPSK modulation
    mxsig = pskmod(bits, 4, pi/8); % Modulate using BPSK

    % Visualization (optional)
    figure;
    subplot(2,1,1);
    plot(time, bits);
    title('Original Bits');

    subplot(2,1,2);
    plot(time, mxsig);
    title('BPSK Modulated Signal');

    scatterplot(mxsig)