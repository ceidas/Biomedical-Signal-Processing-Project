%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018 - Biomedical Signal Processing        %
% MSc Biomedical Engineering - upatras       %
% Elissaios Petrai                           %
% petrai AT ceid.upatras.gr                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
clc;
clear;
close all

% reading the signal
[x, fs] = audioread('F3n.wav');

figure('NumberTitle', 'off', 'name', 'Initial Signal')
N = length(x);
t1 = [0:N - 1] / fs;
hold on;
plot(t1, x, 'b');
title('Initial Signal')

%DC removal
x1 = x - mean(x);

%%Powerline removal - implement Notch filter to remove 50Hz
w = 50 / (44100/2);
bw = w;
[num, den] = iirnotch(w, bw); % notch filter implementation
ecg_notch = filter(num, den, x);

N1 = length(ecg_notch);
t1 = [0:N1 - 1] / fs;
plot(t1, ecg_notch, 'r');
title('Signal after IIR Notch ')
legend('ORIGINAL SIGNAL', 'IIR NOTCH FILTERED SIGNAL')
xlabel('time')
ylabel('amplitude')
hold off;

%%WAVELET DENOISING
xd = wden(ecg_notch, 'sqtwolog', 's', 'mln', 100, 'sym5');
figure('NumberTitle', 'off', 'name', 'Final filtered Signal')
subplot(2, 1, 1)
plot(xd);
title('Final filtered signal')
subplot(2, 1, 2)
plot(xd(570000:650000));
title('Final filtered signal plotting for 570000 to 650000 samples')

% %%CD removal on frequency domain
% figure
% fx=fft(xd);
% plot(abs(fx));
% fx(2:50)=0;
% figure
% final=real(ifft(fx));
% plot(final);

% deb = x(1);
% xd = wden(x-ecg_notch,'sqtwolog','s','mln',3,'db6')+deb;
% plot(xd);
