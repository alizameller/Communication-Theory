%% Project 3
% Jacob Khalili, Aliza Meller, Brian Khaimov
clear; clc; close all;

% Generate a sequence of +/-1 symbols
N = 20000;
msg1 = randi([0,1], 1, N); % message with 1's and 0's
msg = msg1;
msg(~msg) = -1; % message with 1's and -1's

% Passthrough Filter
lf = [3, 5, 7];
step = [0.05, 0.08, 0.095];

%MSE Results
MSE1 = zeros(1,N);
MSE2 = zeros(1,N);
MSE3 = zeros(1,N);

SNR = 0; % variable to keep track of whether or not we are using AWGN because we use the same function for both
[MSE1, w1] = LMS(N, msg, lf(1), step(1), MSE1, SNR);
[MSE2, w2] = LMS(N, msg, lf(2), step(2), MSE2, SNR);
[MSE3, w3] = LMS(N, msg, lf(3), step(3), MSE3, SNR);

figure;
MSE_plot = semilogy(1:N,MSE1,1:N,MSE2,1:N,MSE3);
grid on;
MSE_plot(1).LineWidth = 2;
MSE_plot(2).LineWidth = 2;
MSE_plot(3).LineWidth = 2;
title('Mean Squared Error Performance');
xlabel('Number of Symbols');
ylabel('Mean Squared Error');
legend({'Length = 3','Length = 5','Length = 7'},'Location', 'Northeast');

% MSE with AWGN Results
SNR = 1:2:11;
MSE_AWGN1 = zeros(length(SNR),N);
MSE_AWGN2 = zeros(length(SNR),N);
MSE_AWGN3 = zeros(length(SNR),N);
for i = 1:length(SNR)
    MSE_AWGN1(i,:) = LMS(N,msg,lf(1), 0.1, MSE_AWGN1(i,:), SNR(i));
    MSE_AWGN2(i,:) = LMS(N,msg,lf(2), 0.1, MSE_AWGN2(i,:), SNR(i));
    MSE_AWGN3(i,:) = LMS(N,msg,lf(3), 0.1, MSE_AWGN3(i,:), SNR(i));
end

% Plot the MSE Results
figure;
s1 = semilogy(MSE_AWGN1.');
s1(1).LineWidth = 2;s1(2).LineWidth = 2;s1(3).LineWidth = 2;s1(4).LineWidth = 2;s1(5).LineWidth = 2;s1(6).LineWidth = 2;
grid on;
title('Mean Squared Error Performance (Filter Length = 3)');
xlabel('Number of Symbols');
ylabel('Mean Squared Error');
legend({'SNR = 1dB','SNR = 3dB','SNR = 5dB','SNR = 7dB','SNR = 9dB','SNR = 11dB'},'Location', 'Northeast');

% Plot the MSE Results
figure;
s2 = semilogy(MSE_AWGN2.');
s2(1).LineWidth = 2; s2(2).LineWidth = 2; s2(3).LineWidth = 2;s2(4).LineWidth = 2; s2(5).LineWidth = 2;s2(6).LineWidth = 2;
grid on;
title('Mean Squared Error Performance (Filter Length = 5)');
xlabel('Number of Symbols');
ylabel('Mean Squared Error');
legend({'SNR = 1dB','SNR = 3dB','SNR = 5dB','SNR = 7dB','SNR = 9dB','SNR = 11dB'},'Location', 'Northeast');

% Plot the MSE Results
figure;
s3 = semilogy(MSE_AWGN3.');
s3(1).LineWidth = 2; s3(2).LineWidth = 2; s3(3).LineWidth = 2;s3(4).LineWidth = 2; s3(5).LineWidth = 2;s3(6).LineWidth = 2;
grid on;
title('Mean Squared Error Performance (Filter Length = 7)');
xlabel('Number of Symbols');
ylabel('Mean Squared Error');
legend({'SNR = 1dB','SNR = 3dB','SNR = 5dB','SNR = 7dB','SNR = 9dB','SNR = 11dB'},'Location', 'Northeast');

% BER
SNR = 1:2:20;
[t1, y1_ber] = BER(w1, msg.', SNR);
received1 = y1_ber > 0;
berw1 = sum(received1 ~= msg1.')/N;

[t2, y2_ber] = BER(w2, msg.', SNR);
received2 = y2_ber > 0;
berw2 = sum(received2 ~= msg1.')/N;

[t3, y3_ber] = BER(w3, msg.', SNR);
received3 = y3_ber > 0;
berw3 = sum(received3 ~= msg1.')/N;

figure();
s = semilogy(SNR, berw1, SNR, berw2, SNR, berw3);
grid on;
s(1).LineWidth = 2;s(2).LineWidth = 2;s(3).LineWidth = 2;
title("BER of LMS for different lengths");
xlabel("SNR (dB)");
ylabel("Bit Error Rate (ratio)");
legend({'Length = 3','Length = 5','Length = 7'},'Location', 'Northeast');

% function defs for BER and LMS calculations
function [MSE,w] = LMS(N,msg,l,step, MSE, SNR)
            z = conv(msg, [1 0.2 0.4]);
            z = z(1:length(msg));

       if SNR ~= 0
            noise = awgn(z + eps*1i, SNR);
            z = real(noise);
       end
        
        % Initialize the algorithm by setting tap weights, w = 0
        w = zeros(1,l);
        
        % Initialize the error signal variable 
        e = zeros(1,N);
        
        % Step size parameter of the adaptive filter
        mu = step;
        
        for i = l:N
            % y is the ouput of inner product of d_n and w
            y = w*z(i:-1:i-l+1).';
            % e is the error between the desired signal and the actual response 
            e(i) = z(i) - y;
            % Adaptations of the tap weights of the equalizer
            w = w + mu*e(i)*z(i:-1:i-l+1);
        end

        % Calculate Mean Squared Error
        for i = 1:length(e)
            MSE(i) = sum(e(1:i).^2)/numel(e(1:i));
        end
end

function [r , y] = BER(w, t, SNR)
    y = zeros(length(t), 10);

    % using deterimined filter coefficient attempt to retrieve the signal
    for i = 1:length(SNR)
        r = awgn(filter([1 0.2 0.4], 1, t), SNR(i));
        y(:,i) = filter(w, 1, r);

    end
end

% Given enough symbols, the LMS algorithm is able to produce a low mean
% sqaured error of about 5e-3. When the legnth = 5, the mse is lowest. When
% we add white gaussian noise, however, as the SNR increases, the MSE
% decreases, but the MSE is lowest (i.e. better) without the noise. In
% addition, the BER graph shows that overall as the legnth increases, the 
% BER decreases. 

