load sampleEEGdata.mat
times = EEG.times;

% simulation parameters
sfreq = 7;
ntrials = 10;
srate = 256;
time  = -1:1/srate:1.5;
gwidth = 0.2;
noiseamp = 0.1; % noise standard deviation

signal = zeros(ntrials, length(time));
signal2 = zeros(ntrials, length(time));

% wavelet... more on this in Chapters 13 and 14
wavetime   = -1:1/EEG.srate:1;
n_conv     = length(wavetime)+length(signal)-1;
waveletfft = fft(exp(2*1i*pi*7.*wavetime) .* exp(-wavetime.^2./(2*(5/(2*pi*10))^2))/7,n_conv);
data10hz   = zeros(ntrials,length(signal));
data10hz2   = zeros(ntrials,length(signal));


phasedata10hz   = zeros(ntrials,length(signal));
phasedata10hz2   = zeros(ntrials,length(signal));


for ti = 1 :ntrials
    
    swave  = 3*sin(2*pi*sfreq*time);
    gausw  = exp( -4*log(2)*(time).^2 / gwidth^2 );
    signal(ti,:) = swave .* gausw;
    
    
    for i=101:length(signal)
        
        signal2(ti,i) = signal(ti,i-100);
        
    end
    
    
    convolution_result_fft = ifft(waveletfft.*fft(signal(ti,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft = convolution_result_fft(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz(ti,:) = abs(convolution_result_fft)*2;
    phasedata10hz(ti,:) = angle(convolution_result_fft);

    
    convolution_result_fft2 = ifft(waveletfft.*fft(signal2(ti,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft2 = convolution_result_fft2(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz2(ti,:) = abs(convolution_result_fft2)*2;
    phasedata10hz2(ti,:) = angle(convolution_result_fft2);

    
    
end

% create data and multiple channels plus noise
data1 = signal ;
% + noiseamp * randn(ntrials,length(time));
data2 = signal2 ;
% + noiseamp * randn(ntrials,length(time));


xx = squeeze(mean(data10hz(:,1:end-1),1));
yy = squeeze(mean(data10hz2(:,1:end-1),1));


h = figure;
set(gcf,'color','w');

subplot(221)
plot(times,squeeze(mean(data1(:,1:end-1),1)),'color','k','linew',2); hold on; plot(times, xx,'r:','linew',2)
hold on
plot(times,squeeze(mean(data2(:,1:end-1),1))-10,'k','linew',2);hold on; plot(times, yy-10,'r:','linew',2)
legend({'Amplitude', 'Power'})

% xlim([-600 700])
xlim([-600 600])

set(gca, 'ytick',[], 'xtick',[], 'fontsize',20)
title({'Sequential','Directed- A->B ≠ B->A'}) % axis square
set(gca,'XColor','none', 'YColor','none')
legend boxoff                   % Hides the legend's axes (legend border and background)

angles1 = angle(hilbert( data1' ));
angles2 = angle(hilbert( data2' ));


% show phase angle time series
subplot(223)
plot(EEG.times,angles1(1:end-1,1),'k','linew',2)
hold on
plot(EEG.times,angles2(1:end-1,1)-10,'k','linew',2)
legend({'Phase'})
xlim([-800 800])

set(gca, 'ytick', [],'xtick', [], 'fontsize', 20)
set(gca,'XColor','none', 'YColor','none')

% h = findobj('type', 'axes');  % Find all sets of axes
% set(h(1), 'visible', 'off') 
legend boxoff                   % Hides the legend's axes (legend border and background)


%% 


load sampleEEGdata.mat
times = EEG.times;

% simulation parameters
sfreq = 7;
ntrials = 10;
srate = 256;
time  = -1:1/srate:1.5;
gwidth = 0.2;
noiseamp = 0.1; % noise standard deviation
data10hz3   = zeros(ntrials,length(signal));


signal = zeros(ntrials, length(time));


for ti = 1 :ntrials
    
    swave  = 3*sin(2*pi*sfreq*time);
    gausw  = exp( -4*log(2)*(time).^2 / gwidth^2 );
    signal(ti,:) = swave .* gausw;
    
    
    convolution_result_fft2 = ifft(waveletfft.*fft(signal(ti,:),n_conv)) * sqrt(5/(2*pi*10));
    convolution_result_fft2 = convolution_result_fft2(floor(length(wavetime)/2)+1:end-floor(length(wavetime)/2));
    data10hz3(ti,:) = abs(convolution_result_fft2)*2;
    
end

% create data and multiple channels plus noise
data1 = signal ;
% + noiseamp * randn(ntrials,length(time));
data2 = signal;
% + noiseamp * randn(ntrials,length(time));

xx = squeeze(mean(data10hz3(:,1:end-1),1));
yy = xx;


subplot(222)
plot(times,squeeze(mean(data1(:,1:end-1),1)),'color','k','linew',2); hold on; plot(times, xx,'r:','linew',2)
hold on
plot(times,squeeze(mean(data2(:,1:end-1),1))-10,'k','linew',2);hold on; plot(times, yy-10,'r:','linew',2)
% xlim([-500 500])
xlim([-600 600])

set(gca, 'ytick',[], 'xtick',[], 'fontsize',20)
title({'Concurrent','Non-directed- A->B = B->A'}) % axis square
set(gca,'XColor','none', 'YColor','none')



angles1 = angle(hilbert( data1' ));
angles2 = angle(hilbert( data2' ));

% show phase angle time series
subplot(224)
plot(EEG.times,angles1(1:end-1,1),'k','linew',2)
hold on
plot(EEG.times,angles2(1:end-1,1)-10,'k','linew',2)
% ylabel('Phase angle')
% title('Phase angle over time')
xlim([-800 800])
set(gca, 'ytick', [],'xtick', [], 'fontsize', 20)
set(gca,'XColor','none', 'YColor','none')

