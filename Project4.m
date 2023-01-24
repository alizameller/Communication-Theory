% Jacob Khalili, Aliza Meller, Brian Khaimov
clear; clc; close all;


N = 1e7;
bits = randi(1,N,1);

SNRs = 1:1:50;
results = zeros(length(SNRs),4);

h1 = (1/sqrt(2))*(randn(N,1)+1i*randn(N,1));
h2 = (1/sqrt(2))*(randn(N,1)+1i*randn(N,1));
h3 = (1/sqrt(2))*(randn(N,1)+1i*randn(N,1));
h4 = (1/sqrt(2))*(randn(N,1)+1i*randn(N,1));

mod_input1 = pskmod(bits,2);
mod_input2 = pskmod(bits,2);

arrayMRRC = zeros(length(SNRs),N);
berMRRC2X = zeros(1, length(SNRs));
berMRRC4X = zeros(1, length(SNRs));
berMRRC1X = zeros(1, length(SNRs));
berNew2T = zeros(1,length(SNRs));
berNew2R = zeros(1,length(SNRs));

for SNR = SNRs

    ro1 = awgn(h1.*mod_input1, SNR);
    ro2 = awgn(h2.*mod_input1, SNR);
    ro3 = awgn(h3.*mod_input1, SNR);
    ro4 = awgn(h4.*mod_input1, SNR);
    
    so1 = conj(h1) .* ro1;
    so2 = conj(h1) .* ro1 + conj(h2) .* ro2;
    so4 = conj(h1) .* ro1 + conj(h2) .* ro2 + conj(h2) .* ro3 + conj(h3) .* ro3 + conj(h4) .* ro4;
    
    rc1 = 0.0 + real(so1) < 0;
    rc2 = 0.0 + real(so2) < 0;
    rc4 = 0.0 + real(so4) < 0;

    [c1, berMRRC1X(SNR)] = biterr(bits, rc1);
    [c2, berMRRC2X(SNR)] = biterr(bits, rc2);
    [c3, berMRRC4X(SNR)] = biterr(bits, rc4);


    rc21 = awgn(h1 .* (mod_input1/(sqrt(2))) + h2 .* (mod_input2/(sqrt(2))), SNR);
    rc22 = awgn(- h1 .* conj(mod_input1/(sqrt(2))) + h2 .* conj(mod_input2/(sqrt(2))), SNR);
    rc23 = awgn(h3 .* (mod_input1/(sqrt(2))) + h4 .* (mod_input2/(sqrt(2))), SNR);
    rc24 = awgn(- h3 .* conj(mod_input1/(sqrt(2))) + h4 .* conj(mod_input2/(sqrt(2))), SNR);

    rcn1 = conj(h1) .* rc21 + h2 .* conj(rc22); 
    rcn2 = conj(h1) .* rc21 + h2 .* conj(rc22) + conj(h3) .* rc23 + h4 .*conj(rc24);

    rc2T = 0.0 + real(rcn1) < 0;
    rc2R = 0.0 + real(rcn2) < 0;

    [c4, berNew2T(SNR)] = biterr(bits, rc2T);
    [c5, berNew2R(SNR)] = biterr(bits, rc2R);
end



semilogy(SNRs,berMRRC1X, SNRs, berMRRC2X, SNRs,  berMRRC4X, SNRs, berNew2T, SNRs, berNew2R);
legend ('no diversity', 'MRCC 1T 2X', 'MRCC 1T 4X', 'New Scheme 2T 1X', 'New Scheme 2T 2X')
