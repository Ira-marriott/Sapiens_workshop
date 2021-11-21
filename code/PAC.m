


fs = 2e3;
t = 0:1/fs:0.3-1/fs;

l = [0 130.81 146.83 164.81 174.61 196.00 220 246.94];
m = [0 261.63 293.66 329.63 349.23 392.00 440 493.88];
h = [0 523.25 587.33 659.25 698.46 783.99 880 987.77];
note = @(f,g) [1 1 1]*sin(2*pi*[l(g) m(g) h(f)]'.*t);

mel = [3 2 1 2 3 3 3 0 2 2 2 0 3 5 5 0 3 2 1 2 3 3 3 3 2 2 3 2 1]+1;
acc = [3 0 5 0 3 0 3 3 2 0 2 2 3 0 5 5 3 0 5 0 3 3 3 0 2 2 3 0 1]+1;

song = [];
for kj = 1:length(mel)
    song = [song note(mel(kj),acc(kj)) zeros(1,0.01*fs)];
end
song = song/(max(abs(song))+0.1);
% sound(song,fs)

figure;plot(song)

figure; pspectrum(song,fs,'spectrogram','TimeResolution',0.31, ...
    'OverlapPercent',0,'MinThreshold',-60)


%% 

fs = 2000;

freqs = [270 425];

order = round(7*fs/freqs(2) );

slope = 0.2;

nyquist = fs/2;

frequencies = [0, [freqs(1)-freqs(1)*slope, freqs(1), freqs(2), freqs(2)+freqs(2)*slope]/nyquist 1]; 

amplitude = [0 0 1 1 0 0];

filter = firls(order, frequencies, amplitude);

fvtool(filter,1)

nbsong = filtfilt(filter,1,song);

figure; pspectrum(nbsong,fs,'spectrogram','TimeResolution',0.31, ...
    'OverlapPercent',0,'MinThreshold',-60)

% sound(nbsong,fs)

figure;plot(nbsong)




%% 


h = hilbert(nbsong);
figure;plot(abs(h).^2)

dur = length(song)/fs;


wavetime = -4.4:1/fs:4.6;

a = 0.1 * sin(2*pi*1.6*wavetime);

% figure;plot(a)

figure;plot(abs(h).^2, 'k','linew' ,2); hold on; plot(a, 'linew' ,2)
legend({'~350Hz power', '2Hz phase'})
xlim([1 10000])

h2 = hilbert(a);
phase = angle(h2(1:end-21));
pwr = abs(h).^2;

pac = abs(mean(pwr.*exp(1i*phase)));

% plot(phase);hold on; plot(pwr)



%% 



% observed cross-frequency-coupling (note the similarity to Euler's formula)
obsPAC = abs(mean(pwr.*exp(1i*phase)));
% obsPAC_bias = abs(mean(power_bias.*exp(1i*phase_bias)));

num_iter = 1000;

permutedPAC = zeros(2,num_iter);

for i=1:num_iter
    
    % select random time point
    random_timepoint = randsample(round(length(eeg)*.8),1)+round(length(eeg)*.1);
    
    % shuffle power
    timeshiftedpwr      = [ pwr(random_timepoint:end) pwr(1:random_timepoint-1) ];
    
    % compute PAC
    permutedPAC(1,i) = abs(mean(timeshiftedpwr.*exp(1i*phase)));
end

% compute PACz
pacz(1) = (obsPAC-mean(permutedPAC(1,:)))/std(permutedPAC(1,:));
% pacz(2) = (obsPAC_bias-mean(permutedPAC(2,:)))/std(permutedPAC(2,:));
zval = norminv(1-(.05/length(pacz)));


figure
% subplot(121)
% plot(abs(h).^2); hold on; plot(a)


% subplot(122)
hist(permutedPAC(1,:),50);
hold on
plot([obsPAC obsPAC],get(gca,'ylim')/2,'m','linew',3)
legend({'Permuted values';'Observed value'})
xlabel('Modulation strength'), ylabel('Number of observations')
title([ 'PAC_z = ' num2str(pacz(1)) ])