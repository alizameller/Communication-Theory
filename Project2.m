% Jacob Khalili, Aliza Meller, Brian Khaimov
clear;
close all;
clc;

SNRs = 0:10; %SNR test values
num_bits = 1e6;

%Generating input bits
data = randi([0,1],num_bits,length(SNRs));
data_temp = data;
data_temp(data==0) = -1;

% Mapping data to a square pulse
rect_pulse = rectpulse(data_temp, 8) + eps*1j;

% Adding white noise
noisy_pulse = zeros(num_bits*8,length(SNRs));
parfor i = 1:length(SNRs)
    noisy_pulse(:,i) = awgn(rect_pulse(:,i),SNRs(i) + 10*log10(1/8));
end
% Putting the new noisy function into a matched filter
received1 = intdump(noisy_pulse, 8);
actually1 = output(received1, data, num_bits);
% Creating the theoretical data
theorical = berawgn(SNRs, 'psk', 2, 'nondiff');
% Ploting data against theoretical data
semilogy(SNRs, theorical, "r")
hold on
semilogy(SNRs, actually1, "g")

grid on
legend

hold off

% Function to calculate the vit error rate given a threshold
function actually = output(received, data, num_bits)
    lambda = 0;
    received = (received > lambda);
    errors = sum(data ~= received);
    actually = sum(errors,1)./num_bits;
end

% After plotting the results of our data vs. the theoretical data, and
% adding the proper scaling factor, they appear to map almost identicaly.
% We also see that the higher the SNR is the less likely we are to produce
% errors.