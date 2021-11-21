% Workshop part three


% Task 1: In groups work together to compute two measures of non-directed/ concurrent entrainment

%% Part one: pick two electrodes and compute correlations in power and in phase
% (note phase is circular data)/ PLV – compute this over all time frequency points and produce time frequency plot of result.



% load in data here
load sampleEEGdata.mat


% intialise varibles for time frequency decomposition

data = EEG.data;                             % take data from eeglab structure
freqs = linspace(2, 20, 18);                 % frequencies used in time frequency decomposition
srate = EEG.srate;                            % sampling rate of EEG data

%  extract power and phase using wavelet convolution


% pick two electrodes
channel1 = 'o1';
channel2 = 'o2';

chan1 = find(strcmpi({EEG.chanlocs.labels}, channel1)==1);
chan2 = find(strcmpi({EEG.chanlocs.labels}, channel2)==1);


% get data for selected channels only 
data1 = squeeze(EEG.data(chan1,:,:));
data2 = squeeze(EEG.data(chan2,:,:));


% define filter parameters 
freqbloom = 1.5; % defines the width of the filter in the frequency domain


% Filter the data over specified frequency range
tf_res = FiltHilb(data1, freqs, freqbloom, srate);
tf_res2 = FiltHilb(data2, freqs, freqbloom, srate);


% permute data as hilbert function computes hilbert over first dimension
temphildata =  hilbert(permute(tf_res,[2,1,3]));
hildata =  permute(temphildata, [2,1,3]);

temphildata2 =  hilbert(permute(tf_res2,[2,1,3]));
hildata2 =  permute(temphildata2, [2,1,3]);


% re permtue data back into original format
pwr = abs(hildata).^2;
pwr2 = abs(hildata2).^2;


% define baseline times for power normalisation
bt = dsearchn(EEG.times',[-700 -400]');

% decibel conversion/ baseline normalisation of power
baseline_power1 = mean(pwr(:,bt(1):bt(2),:),2);  dbpwr = 10*log10( bsxfun(@rdivide,pwr,baseline_power1));
baseline_power2 = mean(pwr2(:,bt(1):bt(2),:),2);  dbpwr2 = 10*log10( bsxfun(@rdivide,pwr2,baseline_power2));


% get phase angle time series
phasedata = angle(hildata);
phasedata2 = angle(hildata2);



% plot time frequency decibel normalised power
figure
contourf(EEG.times, freqs, squeeze(mean(dbpwr,3)), 40, 'linecolor', 'non');
xlim([-500 800])
colormap jet
title('Power')


figure;
contourf(EEG.times, freqs, abs(mean(exp(1i*phasedata),3)), 40, 'linecolor', 'non');
xlim([-500 800])
colormap jet
title('ITC')


%% entrainment part one power correlations


corr_ts = tfPow_corr(dbpwr,dbpwr2);

figure
subplot(221)
contourf(EEG.times, freqs, squeeze(mean(dbpwr,3)), 40, 'linecolor', 'non');
xlim([-500 800])
colormap jet
colorbar

subplot(222)
contourf(EEG.times, freqs, squeeze(mean(dbpwr,3)), 40, 'linecolor', 'non');
xlim([-500 800])
colormap jet
colorbar

subplot(2,2,[3 4])
contourf(EEG.times, freqs, squeeze(mean(corr_ts,3)), 40, 'linecolor', 'non');
xlim([-500 800])
% ylim([2 30])
colormap jet


%% part two phase correlations - circular

pnts = EEG.pnts;
trials = EEG.trials;

data2corr1 = phasedata;
data2corr2 = phasedata2;

% tf_plv = abs(mean(exp(1i* (data2corr1 - data2corr2)), 3));

times = -800:10:1000;
times2save = dsearchn(EEG.times', times');

circcorr = zeros(length(freqs), length(times2save));


for i = 1:length(freqs)
    
    for ti = 1:length(times2save)
        
        timeindx = times2save(ti)-32:times2save(ti)+32;
        
        tempdat = reshape(squeeze(data2corr1(i,timeindx,:)), 1, length(timeindx)*trials);
        tempdat2 = reshape(squeeze(data2corr2(i,timeindx,:)), 1, length(timeindx)*trials);
        
        circcorr(i, ti) = circular_corr(tempdat,tempdat2);
        
    end
    
end
        
       

figure
contourf(times, freqs, circcorr, 40, 'linecolor', 'non');
% xlim([-500 800])
colormap jet
colorbar
ylim([2 20])

figure
contourf(EEG.times, freqs, abs(mean(exp(1i* (data2corr1 - data2corr2)), 3)), 40, 'linecolor', 'non');
% xlim([-500 800])
colormap jet
colorbar
ylim([2 20])


        
        
        
 %% maybe? Part three: compute connectivity between all pairs of channels – this will create a channel x channel matrix of connectivity that you can also then visualize
        
 nbchans = EEG.nbchan;
 
 corr_ts = zeros(nbchans, nbchans, length(freqs));
 
 for i = 1:nbchans
     
     for j = 1:nbchans
         
         
         for fi = 1:length(freqs)
         
        
data2corr1 = dbpwr(i,fi,:,:);
data2corr2 = dbpwr(j,fi,:,:);

corr_ts(i,j,fi) = squeeze(mean(tfPow_corr(data2corr1,data2corr2)));


         end
     end
 end
 
 
 
 
 figure; imagesc(squeeze(corr_ts(:,:,8))); colormap jet
 

 
 
 
 
 %% 
 
 % Remove phase locked part of signal (ERP) can recompute connectivity- how has the PLV changed? What does this tell you about how increasing in PLV might be generated?
