% Workshop part three


% Task 1: In groups work together to extract and plot alpha phase and power
% for supplied EEG data



%% Generate song 

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
sound(song,fs)

figure; pspectrum(song,fs,'spectrogram','TimeResolution',0.31, ...
    'OverlapPercent',0,'MinThreshold',-60)


%% low pass filtering

% most of this is the code same code that is already written within the supplied Filthilb
% function that everyone can use but could in short go over it here

% general filtering paramters 
order = round(7*fs/freqs(2) );

slope = 0.2;

nyquist = fs/2;


%% 


freqs = [0 450];

frequencies = [freqs(1), [freqs(2)-freqs(2)*slope, freqs(2)+freqs(2)*slope]/nyquist 1]; 
    
amplitude = [1 1 0 0];

filter = firls(order, frequencies, amplitude);

long = filtfilt(filter,1,song);

figure; pspectrum(long,fs,'spectrogram','TimeResolution',0.31, ...
    'OverlapPercent',0,'MinThreshold',-60)

% sound(long,fs)



%% high pass filtering

freqs = [0 250];

frequencies = [freqs(1), [freqs(2)-freqs(2)*slope, freqs(2)+freqs(2)*slope]/nyquist 1]; 

amplitude = [0 0 1 1];

filter = firls(order, frequencies, amplitude);

lpsong = filtfilt(filter,1,song);

figure; pspectrum(lpsong,fs,'spectrogram','TimeResolution',0.31, ...
    'OverlapPercent',0,'MinThreshold',-60)

% sound(lpsong,fs)


%% bandpass filter

freqs = [270 425];

frequencies = [0, [freqs(1)-freqs(1)*slope, freqs(1), freqs(2), freqs(2)+freqs(2)*slope]/nyquist 1]; 

amplitude = [0 0 1 1 0 0];

filter = firls(order, frequencies, amplitude);


nbsong = filtfilt(filter,1,song);

figure; pspectrum(nbsong,fs,'spectrogram','TimeResolution',0.31, ...
    'OverlapPercent',0,'MinThreshold',-60)




%% 










