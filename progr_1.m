%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018 - Biomedical Signal Processing        %
% MSc Biomedical Engineering - upatras       %
% Elissaios Petrai                           %
% petrai AT ceid.upatras.gr                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
[ecg, f_s] = audioread('118e00m.wav');
N = length(ecg);
ti = [0:N - 1] / f_s; %time period(total sample/Fs )

%%IIR Notch filter to remove powerline (50Hz)
w = 50 / (360/2);
bw = w;
[num, den] = iirnotch(w, bw); % notch filter implementation
ecg_notch = filter(num, den, ecg);

%plot ecg after Notch filter
figure
N1 = length(ecg_notch);

%obtain each signal seperately
x1 = ecg_notch(:, 1);
x2 = ecg_notch(:, 2);
L1 = length(x1);
L2 = length(x2);
x = ecg(:, 1);
hold on
plot(x); title('Channel 1')
xlabel('time')
ylabel('amplitude')

plot(x1, 'r');
legend('Initial Signal of Channel 1', 'Signal of Channel 1 after IIR Notch Filter')
xlabel('time')
ylabel('amplitude')

%wavelet implementation
for i = 1:2

    if i == 1
        j = x1;
    elseif i == 2
        j = x2;
    end

    %%wavelet decomposition of the signal
    [e, l] = wavedec(j, 10, 'sym4'); % Wavelet implementation, e=wavelet decomposition vector, l=bookkeeping vector
    approx = appcoef(e, l, 'sym4');
    [cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10] = detcoef(e, l, [1 2 3 4 5 6 7 8 9 10]);

    %%reduce baseline wandering cd10=a10=0
    a10 = appcoef(e, l, 'sym4', 10);
    a10 = wrcoef('a', e, l, 'sym4', 10);
    cd10 = wrcoef('d', e, l, 'sym4', 10);
    e(1:1290) = 0;

    %%Reduction of EMG noise

    %cd4
    elem4 = sum(l(1:7) + 1);
    elem41 = sum(l(1:8));
    [px1, ax1] = ecdf(e(elem4:elem41)); %elements of cd4
    lo_th_1 = floor(length(ax1) * 15/100);
    up_th_1 = floor(length(ax1) * 85/100);
    pinakas = e(elem4:elem41);
    pinakas(pinakas > ax1(lo_th_1) & pinakas < ax1(up_th_1)) = 0; %hard thresholding
    cd4 = pinakas;
    e(elem4:elem41) = pinakas;

    %cd3
    elem3 = sum(l(1:8) + 1);
    elem31 = sum(l(1:9));
    [px2, ax2] = ecdf(e(elem3:elem31));
    lo_th_2 = floor(length(ax2) * 10/100);
    up_th_2 = floor(length(ax2) * 90/100);
    pinakas2 = e(elem3:elem31);
    pinakas2(pinakas2 > ax2(lo_th_2) & pinakas2 < ax2(up_th_2)) = 0; %hard thresholding
    cd3 = pinakas2;
    e(elem3:elem31) = pinakas2;

    %%reduce of motion artifacts. Remove the very large magnitude cD8 , cD9 coefficients above threshold
    %cd9
    s9 = median(abs(cd9)) / 0.6457;
    t9 = s9 * sqrt(2 * log10(l(3)));
    cd9 = wthresh(cd9, 's', t9); %soft thresholding

    %cd8
    s8 = median(abs(cd8)) / 0.6457;
    t8 = s8 * sqrt(2 * log10(l(4)));
    cd8 = wthresh(cd8, 's', t8); %soft thresholding

    %%high frequencies elimination - cd1=0; & cd2=0;
    a = sum(l(1:9)) + 1; %162,594
    t = sum(l(1:11)); %650,106
    e(a:t) = 0;

    if i == 1
        figure;
        plot(ax1, px1);
        line([ax1(lo_th_1) ax1(lo_th_1)], [0 1], 'Color', 'blue');
        line([ax1(up_th_1) ax1(up_th_1)], [0 1], 'Color', 'blue');
        title('First Channel CDF - cd4')

        figure;
        plot(ax2, px2);
        line([ax2(lo_th_2) ax2(lo_th_2)], [0 1], 'Color', 'blue');
        line([ax2(up_th_2) ax2(up_th_2)], [0 1], 'Color', 'blue');
        title('First Channel CDF - cd3')

        %reconstruct the signal
        signal_rec1 = waverec(e, l, 'sym4');
        signal_rec1 = smooth(signal_rec1);
    else
        figure;
        plot(ax1, px1);
        line([ax1(lo_th_1) ax1(lo_th_1)], [0 1], 'Color', 'blue');
        line([ax1(up_th_1) ax1(up_th_1)], [0 1], 'Color', 'blue');
        title('Second Channel CDF - cd4')

        figure;
        plot(ax2, px2);
        line([ax2(lo_th_2) ax2(lo_th_2)], [0 1], 'Color', 'blue');
        line([ax2(up_th_2) ax2(up_th_2)], [0 1], 'Color', 'blue');
        title('Second Channel CDF - cd3')

        %reconstruct the signal
        signal_rec2 = waverec(e, l, 'sym4');
        signal_rec2 = smooth(signal_rec2); % using average filter to remove glitches
        %to increase the performance of peak detection

    end

end

%Concatenate
final_signal = [signal_rec1 signal_rec2];

%plot final and initial signal
figure
subplot(2, 1, 1)
plot(ecg(1:2000, :)); title('Raw ECG Data plotting ')
xlabel('samples')
ylabel('amplitude')

subplot(2, 1, 2)
plot(final_signal(1:2000, :)); title('Final Signal Plotting ')
xlabel('samples')
ylabel('amplitude')

figure
subplot 211
plot(final_signal)
title('Baseline Wander Removal')
subplot 212
plot(final_signal(1:2000))
title('Baseline Wander Removal')
legend('Samples 1:2000')
