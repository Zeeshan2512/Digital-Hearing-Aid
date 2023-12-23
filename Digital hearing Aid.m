clear all;
close all;
clc;
% Input sound
input = 'Sample/late_noise.wav';
% Medical Parameters
g = 50;
transitionV = [1000, 1500, 2550, 5000];
Psat=100;
% Sound Output
y = hearingAidF(input,g,Psat,transitionV);
function y = denoise(x,fs)
% y = denoise(x);
% method to denoise a given signal using wavelets
% x is the input Matlab sound file
%THR is the threshold, SORH is for soft or hard thresholding, KEEPAPP allows you to keep
% approximation coefficients
[thr,sorh,keepapp]=ddencmp('den','wv',x);
% returns a de-noised version xc of input signal x (our one-dimensional speech signal)
[y, ~, ~, ~, ~]=wdencmp('gbl',x,'db3',2,thr,sorh,keepapp);
% comparing noisy and denoise signal
x_length = length(x);
t=(0:1/fs:(x_length-1)/fs);
figure;
subplot(2,1,1);
plot(t,x,'r');
title('Signal with Noise');
subplot(2,1,2);
plot(t,y);
title('Signal after Denoising');
function y = freqshape(x,g,transitionV,fs)
% y = applySkiSlope(x,g,transitionV,fs)
% Creates the gain filter for a patient with ski slope hearing loss.
% The maximum gain will be g and the minimum gain will be one. The magnitude
% of the gain function will be the concatenation of preset piecewise functions.
% However the time of the transitions from one piecewise function to another can
% be set by the user in the elements of the transitionV. The final frequency used
% will be fs/2 since that's the highest frequency that the input signal will contain.
% The output will be the filtered signal
first = transitionV(1);
second = transitionV(2);
third = transitionV(3);
fourth = transitionV(4);
x_length = length(x);
n = nextpow2(x_length);
N = 2^n;
T = 1/fs;
X = fft(x,N);
gain = zeros(N,1);
% Sets the gain for the first stage of frequencies
firstC = (.3*(g-1))/first;
k=0;
while(k/N <= first/fs)
gain(k+1) = firstC*k/(N*T) + 1;
gain(N-k) = gain(k+1);
k=k+1;
end;
% Sets the gain for the second stage of frequencies
secondC = firstC*first +1;
secondC2 = (second-first)/5;
while(k/N <= second/fs)
gain(k+1) = 1 + (secondC-1)*exp(-((k/(N*T))-first)/secondC2);
gain(N-k) = gain(k+1);
k=k+1;
end;
% Sets the gain for the third stage of frequencies
thirdC = 1 + (secondC-1)*exp(-second/secondC2);
thirdC2 = (third-second)/5;
while(k/N <= third/fs)
gain(k+1) = g + (thirdC-g)*exp(-((k/(N*T)-second))/thirdC2);
gain(N-k) = gain(k+1);
k=k+1;
end;
% Sets the gain for the fourth stage of frequencies
while(k/N <= fourth/fs)
gain(k+1) = g;
gain(N-k) = gain(k+1);
k=k+1;
end;
% Sets the gain for the fifth stage of frequencies
fifthC = g;
fifthC2 = (fs/2-fourth)/5;
while(k/N <= .5)
gain(k+1) = 1 + (fifthC-1)*exp(-((k/(N*T))-fourth)/fifthC2);
gain(N-k) = gain(k+1);
k=k+1;
end;
k_v = (0:N-1)/N;
figure; %non-redundant filter transfer function
k_v = k_v*fs;
k_v = k_v(1:N/2+1);
plot(k_v,gain(1:N/2+1));
title('Frequency Shaper Transfer Function');
xlabel('Frequency (Hertz)');
ylabel('Gain');
xlim([0 10000]);
Y = X+gain; % for X refer line no.27
y = real(ifft(Y,N));
y = y(1:x_length);
function y = powerCompress(input, Psat,Fs)
% y = powerCompress(input, Psat,Fs)
% Takes in a a signal makes sure that the maximum power in any frequency
% is less than or equal to Psat. Also had some denoising capabilities, by
% zeroing out very low power frequencies.
% input - input Matlab sound file
% Psat - Saturation power
% FS - Sampling frequency of the input signal
x=input;
len=Fs*0.1;
iter=floor(length(x)/len);
Plow=0.008;
for rg=0:1:iter;
start=rg*len+1;
en= rg*len+len;
if rg*len+len>length(x)
en=length(x);
end
clear signal X X_pow Y_pow Y y z;
signal=x(start:en);
n = nextpow2(len);
N = 2^n;
X = fft(signal,N);
X_phase=angle(X); % Save the old phase information
X_pow = abs(X)/N;
Y_pow = X_pow;
Y=zeros(N,1);
for k=0:N/2
if Y_pow(k+1)<Plow % Take out noise
Y_pow(k+1)=0;
Y_pow(N-k)=0;
elseif Y_pow(k+1)>Psat % Clip amplitudes higher than Psat
Y_pow(k+1)=Psat;
Y_pow(N-k)=Psat;
end;
Y(k+1) = Y_pow(k+1)*(cos(X_phase(k+1))+1i*sin(X_phase(k+1)));
Y(N-k) = Y_pow(N-k)*(cos(X_phase(N-k))+1i*sin(X_phase(N-k)));
end;
y = real(ifft(Y,N));
z = y(1:en-start+1);
sig_out(start:en)=z;
end;
y = sig_out*5000; % Multiplying 5000 is just increasing the intensity of the o/p signal
function y = hearingAidF(input,g,Psat,transitionV)
% input - the input signal to the system. Should be a wave file.
% g - the maximum gain that will be applied to the signal
% Psat - the cut-off power. The output power will not be higher than this
% transitionV - 4 element vector that has the values of where the gain changes
% to the next piecewise function
% Inputing and reading audio file
[x,fs] = audioread(input);
x = x(:, 1);
% Add Noise to the signal (Step 1)
x = awgn(x,20);
% Denoising filter (Step 2)
xc = denoise(x,fs);
% Frequency shaping filter (Step 3)
xf = freqshape(xc,g,transitionV,fs);
% Amplitude shaping filter (Step 4)
y = powerCompress(xf, Psat,fs);
% Comparing Spectrograms
figure;
subplot(2,1,1);
specgram(x);
title('Spectrogram of Original Signal 2');
subplot(2,1,2);
specgram(y);
title('Spectrogram of Adjusted Signal 2');
% Output sound (Step 5)
disp('Adjusted Sound');
sound(y,fs);