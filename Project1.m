% Aliza Meller, Brian Khaimov and Jacob Khalili

clc;
close all;
clear;

% rate 7-4 BCH code
m = 3;
k = 4;
n = 2^m - 1;
%nwords = 250000;
nwords = 25000;

msgTx = gf(randi([0 1], nwords, k)); %Generating the trancieved message

probs = 0:0.05:0.5; %Testing for a range of probabilities
encode = bchenc(msgTx, n, k);
errors_74 = zeros(1,length(probs));
i = 1;

for prob = probs %Testing decoding efectiveness for a range of probabilities
    noisycode = bsc(encode, prob); %Flipping bits with given probability
    decoded = bchdec(noisycode, n, k);
    error_count = accuracy(decoded, msgTx); %Calls accuracy function defined below
    errors_74(i) =  error_count/(k*nwords); %Calculates the bit error rate
    i = i + 1;
end

% BCH 15/7
%Runs similar to BCH code above
k = 7;
m = 4;
n = 2^m - 1;
%nwords = 142858;
nwords = 14286;


msgTx = gf(randi([0 1], nwords, k));

probs = 0:0.05:0.5;
encode = bchenc(msgTx, n, k);
errors_157 = zeros(1,length(probs));
i = 1;

for prob = probs
    noisycode = bsc(encode, prob);
    decoded = bchdec(noisycode, n, k);
    error_count = accuracy(decoded, msgTx);
    errors_157(i) =  error_count/(k*nwords);
    i = i + 1;
end



% Convolutional Code 1/2

trellis = poly2trellis(3,[6 7]); % Generates a trellis of length 3
% and describes the code generator using octal values
probs = 0:0.05:0.5;
errors_poly12 = zeros(1,length(probs));
i = 1;

for prob = probs
    %data = randi([0 1],1000000,1); % Generates random input data
    data = randi([0 1],70,1);
    codedData = convenc(data,trellis);
    noisyData = bsc(codedData, prob); % Corrupts data with specified probability
    tbdepth = 34; % Traceback depth
    decodedData = vitdec(noisyData,trellis,tbdepth,'trunc','hard');% Decodes data
    errors_poly12(i) = biterr(data,decodedData)/70; % Computes bit error rate
    i = i + 1;
end

% Convolutional Code 2/3
% Works similar to convolutional code above
trellis = poly2trellis([5 4],[23 35 0; 0 5 13]); % Generates a trellis with upper and lower contraints of 5 and 4 respectively

probs = 0:0.05:0.5;
errors_poly23 = zeros(1,length(probs));
i = 1;

for prob = probs
     %data = randi([0 1],1000000,1);
    data = randi([0 1],70,1);
    codedData = convenc(data,trellis);
    noisyData = bsc(codedData, prob);
    tbdepth = 34;
    decodedData = vitdec(noisyData,trellis,tbdepth,'trunc','hard');
    errors_poly23(i) = biterr(data,decodedData)/70;
    i = i + 1;
end


figure; hold on; %Plot results of code effectiveness on graph
a1 = plot(probs, errors_74, "r"); M1 = "BCH (7-4)";
a2 = plot(probs, errors_157, "b"); M2 = "BCH (15-7)";
a3 = plot(probs, errors_poly12, "g"); M3 = "CONV 1/2";
a4 = plot(probs, errors_poly23, "m"); M4 = "CONV 2/3";

legend([a1,a2,a3,a4],[M1,M2,M3,M4]);



%Function to count ammount of bits flipped in recieved message
function error_count = accuracy(decoded, msgTx)
    change = decoded.x ~= msgTx.x;
    error_count = sum(sum(change));
end

%Summmary of results: We learned how to implement BCH code and
%conviolutional code in matlab, which helped us understand each step of
%the process. All codes were pretty unreliable in decoding the recieved
%message when give large error porbabilities(<0.2). The convoltional 2/3
%code behaved strangely at eroor rates of about 0.2, the bit error rate
%usually shot up around then. The rest of the codes trended upward mostly
%relative to the error probability.
