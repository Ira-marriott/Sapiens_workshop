% Workshop part three


% Task 2:

% part one pick two electrodes and compute granger causality for one time
% segement - choose whichever time window length you want

% part two repeat this process for all time segement and visualise the
% resulting time varying granger causality

% compute time varing granger causality for first and last 40 trials.
% visualise the results in time and space - imagine these are two
% experimental conditions - dicuss in your groups any differences between
% these two groups in time or space and how you would interpret these
% observations.



%% Example one  - here we are going to fit AR models to the infant and mum arousal datasets that we used in the previous section.
% then we are going to calculate granger causality between the infant and
% mum arousal data

% your task for this example is just to dicuss what the interpretation of
% the visualisation is.



FromDir = 'C:\Users\Ira Marriott\Desktop\Sapiens workshop\home awake synched_orig\';
cd(FromDir)
files = dir('*.mat');

load(strcat(FromDir, files(2).name))

morder = 5; % model order determines the number of past samples included in the AR model fit


% univariate AR model fit 
[Ax,Ex] = armorf(ret(:,1)', 1, length(ret), morder);

% univariate AR model fit 
[Ay,Ey] = armorf(ret(:,2)', 1, length(ret), morder);

% bivariate AR model fit
[Axy,E] = armorf(ret', 1, length(ret), morder);


y2x = log(Ex/E(1,1));
x2y = log(Ey/E(2,2));

 var(Ex)/var(E)



figure; bar([x2y, y2x])
set(gca, 'XTickLabel', {'Infant to mum', 'Mum to infant'}, 'fontsize', 15)
ylabel('Granger causal estimate')



%%  part one pick two electrodes and compute granger causality for one time
% segement - choose whichever time window length you want

% load in data here
load sampleEEGdata.mat


% intialise varibles for time frequency decomposition

data = EEG.data;                             % take data from eeglab structure
srate = EEG.srate;                            % sampling rate of EEG data
ntrials = EEG.trials;
pnts = EEG.pnts;

morder = 5;

% pick two electrodes
channel1 = 'o1';
channel2 = 'o2';

chan1 = find(strcmpi({EEG.chanlocs.labels}, channel1)==1);
chan2 = find(strcmpi({EEG.chanlocs.labels}, channel2)==1);


% get data for selected channels only
data1 = squeeze(EEG.data(chan1,:,:));
data2 = squeeze(EEG.data(chan2,:,:));


twin = 500;
twin_pnts = round(twin/(1000/srate));

timepnt = dsearchn(EEG.times', 200');


X = [data1(timepnt-floor(twin_pnts/2):timepnt+floor(twin_pnts/2),:);...
    data2(timepnt-floor(twin_pnts/2):timepnt+floor(twin_pnts/2),:)];

tempdata = reshape(X,2,(twin_pnts+1)*ntrials);



[Ax,Ex] = armorf(tempdata(1,:),ntrials,twin_pnts,morder);
[Ay,Ey] = armorf(tempdata(2,:),ntrials,twin_pnts,morder);
[Axy,E] = armorf(tempdata     ,ntrials,twin_pnts,morder);


y2x = log(Ex/E(1,1));
x2y = log(Ey/E(2,2));

figure; bar([x2y, y2x])
set(gca, 'XTickLabel', {'Channel 1 to channel 2', 'Channel 2 to channel 1'}, 'fontsize', 15)
ylabel('Granger causal estimate')



%% part two repeat this process for all time segement and visualise the
% resulting time varying granger causality

fs = 256;
twin  = 128;
times2save = -800:50:1000;
times2saveidx = dsearchn(EEG.times',times2save');
twin_pnts = round(twin/(1000/fs));

chan1 = 1;
chan2 = 18;

y2x = zeros(1, length(times2save));
x2y = zeros(1, length(times2save));


for i = 1:length(times2save)
    
    
    X = [EEG.data(chan1,times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),:); ...
        EEG.data(chan2, times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),:)];
    
    %    X = permute(X,[3 1 2]);
    
    
    tempdata = reshape(X,2,(twin_pnts)*ntrials);
    
    %             % fit AR models (model estimation from bsmart toolbox)
    [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,twin_pnts,morder);
    [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,twin_pnts,morder);
    [Axy,E] = armorf(tempdata     ,EEG.trials,twin_pnts,morder);
    
    
    %             % time-domain causal estimate
    y2x(i)=log(Ex/E(1,1));
    x2y(i)=log(Ey/E(2,2));
    
end


figure; plot(times2save,x2y)
hold on
plot(times2save,y2x)
xlabel('Time ms')
ylabel('Granger causal estimate')
set(gca,  'fontsize', 15)
legend({'Channel 1- channel 2', 'Channel 2- channel 1'})

% discuss how to interpret these results



%% compute time varing granger causality for first and last 40 trials.
% visualise the results in time and space - imagine these are two
% experimental conditions - dicuss in your groups any differences between
% these two groups in time or space and how you would interpret these
% observations.


condition1 = 1:40;
condition2 = 60:ntrials;


fs = 256;
twin  = 128;
times2save = -800:50:1000;
times2saveidx = dsearchn(EEG.times',times2save');
twin_pnts = round(twin/(1000/fs));

chan1 = 1;
chan2 = 18;

y2x = zeros(2, length(times2save));
x2y = zeros(2, length(times2save));


for i = 1:length(times2save)
    
    
    X = [EEG.data(chan1,times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),condition1); ...
        EEG.data(chan2, times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),condition1)];
    
    
    X2 = [EEG.data(chan1,times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),condition2); ...
        EEG.data(chan2, times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),condition2)];
    
    
    
    tempdata = reshape(X,2,(twin_pnts)*length(condition1));
    tempdata2 = reshape(X2,2,(twin_pnts)*length(condition1));
    
    
    %             % fit AR models (model estimation from bsmart toolbox)
    [~,Ex] = armorf(tempdata(1,:),length(condition1),twin_pnts,morder);
    [~,Ey] = armorf(tempdata(2,:),length(condition1),twin_pnts,morder);
    [~,E] = armorf(tempdata     ,length(condition1),twin_pnts,morder);
    
    
    %             % fit AR models (model estimation from bsmart toolbox)
    [~,Ex2] = armorf(tempdata2(1,:),length(condition1),twin_pnts,morder);
    [~,Ey2] = armorf(tempdata2(2,:),length(condition1),twin_pnts,morder);
    [~,E2] = armorf(tempdata2     ,length(condition1),twin_pnts,morder);
    
    
    
    %             % time-domain causal estimate
    y2x(1,i)=log(Ex/E(1,1));
    x2y(1,i)=log(Ey/E(2,2));
    
    y2x(2,i)=log(Ex2/E2(1,1));
    x2y(2,i)=log(Ey2/E2(2,2));
    
end


figure;

subplot(211);plot(times2save,x2y(1,:))
hold on
plot(times2save,y2x(1,:))
xlabel('Time ms')
ylabel('Granger causal estimate')
set(gca,  'fontsize', 15)
legend({'Channel 1- channel 2', 'Channel 2- channel 1'})



subplot(212); plot(times2save,x2y(2,:))
hold on
plot(times2save,y2x(2,:))
xlabel('Time ms')
ylabel('Granger causal estimate')
set(gca,  'fontsize', 15)
legend({'Channel 1- channel 2', 'Channel 2- channel 1'})







