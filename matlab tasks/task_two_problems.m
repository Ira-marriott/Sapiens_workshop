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




%% example one compute granger causality between infant and mum arousal data


FromDir = 'C:\Users\Ira Marriott\Desktop\Sapiens workshop\home awake synched_orig\';
cd(FromDir)
files = dir('*.mat');

load(strcat(FromDir, files(2).name))

morder = 5;
% model order determines the number of past samples included in the AR model fit


% univariate AR model fit 
[~,Ex] = armorf(ret(:,1)', 1, length(ret), morder);

% univariate AR model fit 
[~,Ey] = armorf(ret(:,2)', 1, length(ret), morder);

% bivariate AR model fit
[~,E] = armorf(ret', 1, length(ret), morder);


y2x = log(Ex/E(1,1));
x2y = log(Ey/E(2,2));


figure; bar([x2y, y2x])
set(gca, 'XTickLabel', {'Infant to mum', 'Mum to infant'}, 'fontsize', 15)
ylabel('Granger causal estimate')





%%  part one pick two electrodes and compute granger causality for one time
% segement - choose whichever time window length you want

% load in data here
load( )

% intialise varibles for time frequency decomposition
data = EEG.data;                             % take data from eeglab structure
srate = EEG.srate;                            % sampling rate of EEG data
ntrials = EEG.trials;
pnts = EEG.pnts;
morder = 5;


% pick two electrodes...
channel1 = ;
channel2 = ;

% find channel numbers
chan1 = find(strcmpi({EEG.chanlocs.labels}, channel1)==1);
chan2 = find(strcmpi({EEG.chanlocs.labels}, channel2)==1);


% get data for selected channels only..
data1 = ;
data2 = ;



% define time window
twin =  ;
twin_pnts = round(twin/(1000/srate));


% pick a time point to look at
timepoint = ;
timepnt = dsearchn(EEG.times', timepoint');


% take data for both channels for the time window and all trials
% X should be a 3d matrix channels (2) by time points by trials
X = ;

% reshape data
tempdata = reshape(X,2,(twin_pnts+1)*ntrials);


% fit univariate AR models for data1 (x) and data2 (y)
[Ax,Ex] = armorf(tempdata(1,:),ntrials,twin_pnts,morder);
[Ay,Ey] = armorf(tempdata(2,:),ntrials,twin_pnts,morder);

% fit bivariate AR models for data1 (x) and data2 (y)
[Axy,E] = armorf(tempdata     ,ntrials,twin_pnts,morder);

% define granger causality
y2x = log(Ex/E(1,1));
x2y = log(Ey/E(2,2));


%% part two repeat this process for all time segement and visualise the
% resulting time varying granger causality


% define time points to loop over
timepointstoloop = -800:50:1000;
times2saveidx = dsearchn(EEG.times',timepointstoloop');


% define time window
twin =  ;
twin_pnts = round(twin/(1000/srate));


% intialise output stuctures
y2x = zeros(1, length(timepointstoloop));
x2y = zeros(1, length(timepointstoloop));


for i = ; %     need to loop over times2saveidx
    
    
    
    % take data for both channels for the time window and all trials
    % X should be a 3d matrix channels (2) by time points by trials
    % take X for only time points in current iteration of loop
    
    X = ;
    
    
    tempdata = reshape(X,2,(twin_pnts)*ntrials);
    
    %             % fit AR models (model estimation from bsmart toolbox)
    [Ax,Ex] = ;
    [Ay,Ey] = ;
    [Axy,E] = ;
    
    
    %             % time-domain causal estimate
    y2x(i) = ;
    x2y(i) = ;
    
end


figure; plot(timepointstoloop,x2y)
figure; plot(timepointstoloop,y2x)



%% compute time varing granger causality for first and last 40 trials.
% visualise the results in time and space - imagine these are two
% experimental conditions - dicuss in your groups any differences between
% these two groups in time or space and how you would interpret these
% observations.



% initalise variable same as before
fs = 256;
twin  = 128;
timepointstoloop = -800:50:1000;
times2saveidx = dsearchn(EEG.times',timepointstoloop');
twin_pnts = round(twin/(1000/fs));



condition1 = ;
condition2 = ;

chan1 = ;
chan2 = ;


% find channel numbers
chan1 = find(strcmpi({EEG.chanlocs.labels}, channel1)==1);
chan2 = find(strcmpi({EEG.chanlocs.labels}, channel2)==1);


% initalise output varibles - not 2 here meaning two conditions
y2x = zeros(2, length(timepointstoloop));
x2y = zeros(2, length(timepointstoloop));


for i = ; %     loop over time points
    
    % take data for both channels for the time window and all trials
    % X should be a 3d matrix channels (2) by time points by trials
    % take X for only time points in current iteration of loop
    
    
    
    % take data only for condition 1 here
    
    X = ;
    
    
    % and then take data for condition 2 here
    X2 = ;
    
    
    %     reshape data
    tempdata = reshape(X,2,(twin_pnts)*length(condition1));
    tempdata2 = reshape(X2,2,(twin_pnts)*length(condition1));
    
    
    %             % fit AR models for condition 1 data
    [~,Ex] = ;
    [~,Ey] = ;
    [~,E] =  ;
    
    
    %             % repeat for condition 2 data
    [~,Ex2] = ;
    [~,Ey2] = ;
    [~,E2] =  ;
    
    
    
    %             % time-domain causal estimate
    y2x(1,i)= ;
    x2y(1,i)= ;
    
    y2x(2,i)= ;
    x2y(2,i)= ;
    
end


figure; plot(timepointstoloop,x2y)
figure; plot(timepointstoloop,y2x)

