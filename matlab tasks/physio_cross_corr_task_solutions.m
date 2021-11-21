

% for this task we are going to look at cross correlation of
% data and interpret these results.

% Part one - using matlabs task one compute cross correlation between physio data 

% Part two - using matlabs medfilt function - filter the arousal data.
% Visualise the filtered arousal data. Compute the cross correlation on the
% filtered data. Visualise this result. Compare the cross correlation
% between the filtered and unfiltered data


% Part two - using matlabs medfilt function - filter the arousal data.
% Visualise the filtered arousal data. Compute the cross correlation on the
% filtered data. Visualise this result. Compare the cross correlation
% between the filtered and unfiltered data

% Part three add random noise to the singals using rand. Vary the amplitude
% of the noise and see how this affects both how the signals look and how
% this effects the resulting cross correlation estimate. Disucss what
% features of the data are chaning here?


% Part four - using matlabs detrend function - detrend the arousal data.
% Visualise the filtered arousal data. Compute the cross correlation on the
% filtered data. Visualise this result. Compare the cross correlation
% between the filtered and unfiltered and detrended data



FromDir = 'C:\Users\Ira Marriott\Desktop\Sapiens workshop\home awake synched_orig\';
cd(FromDir)
files = dir('*.mat');


%% Compute cross correlation as dot product between two vectors 

% intialise two vectors 
dp1 = zeros(3,1);
dp2 = zeros(3,1);

% create two basic signals that will peak in the middle
v1 = [1 1 4 4 1 1];
v2 = [1 1 4 4 1 1];

% visualise these signals
figure; subplot(121); plot(v1);
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)
title('signal X')
subplot(122); plot(v2)
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)
title('signal Y')


% compute the dot product at lag zero between the two signals 
dp1(2) = sum(v1.*v2);
dp2(2) = dot(v1,v2);



% shuffle signal two back one sample and compute the dot product between the two signals at negative 1 lag -
% notice zero padding 
v1 = [0 1 1 4 4 1];
v2 = [1 1 4 4 1 1];

dp1(1) = sum(v1.*v2);
dp2(1) = dot(v1,v2);


% shuffle signal forward one sample and compute the dot product between the two signals at negative 1 lag -
% again notice how zero padding had changed here

v1 = [1 4 4 1 1 0];
v2 = [1 1 4 4 1 1];

dp1(3) = sum(v1.*v2);
dp2(3) = dot(v1,v2);


% visualise cross correlation
figure; plot(dp1, 'o-')
set(gca, 'xticklabels', -1:0.2:1)
xlabel('Lag in samples')
ylabel('Cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation between x and y')


%% example here using matlabs cross correlation function 

v1 = [1 1 4 4 1 1];
v2 = [1 1 4 4 1 1];

[r, lags] = xcorr(v1,v2,1);
% note third function input is number of lags
figure;stem(lags,r)
xlabel('Lag in samples')
ylabel('Cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation using xcorr')



%% % Part one - using matlabs task one compute cross correlation between physio data 

% load in physio data
load(strcat(FromDir, files(2).name))

% take a second to look at the data
whos ret


% remove zeroed out epochs 
x = isnan(ret(:,1));
toremove = find(x==1);
ret(toremove,:) = [];


% visualise mum and infant arousal data 
figure; subplot(211);
plot(ret(1:500,1))
title('Infant arousal')
set(gca, 'fontsize', 15)
set(gca, 'XTick',[], 'YTick',[])
subplot(212); plot(ret(1:500,2))
title('Mum arousal')
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)

% visualise cross correlation between mum and infant arousal data
[r, lags] = xcorr(ret(:,1), ret(:,2),500, 'normalized');
figure;stem(lags,r) 
xlabel('Lag in samples')
ylabel('Normalised cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation in arousal')


%% mean filter the data - reduce variability and see how this effects cross correlation value

% Part two - using matlabs medfilt function - filter the arousal data.
% Visualise the filtered arousal data. Compute the cross correlation on the
% filtered data. Visualise this result. Compare the cross correlation
% between the filtered and unfiltered data



% mean filter the data
ret2 = medfilt1(ret,20);


% visualise the mean filtered data
figure; subplot(211);
plot(ret2(1:500,1))
title('Infant arousal filtered')
set(gca, 'fontsize', 15)
subplot(212); 
plot(ret2(1:500,2))
title('Mum arousal filtered')
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)


% visualise cross correlation between filtered mum and infant arousal data
[r, lags] = xcorr(ret2(:,1), ret2(:,2),500, 'normalized');
figure;stem(lags,r) 
xlabel('Lag in samples')
ylabel('Normalised cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation in filtered arousal')

%% %% explanation on impact of variability 

% define amplitde of noise 
noiseamp = 0.5;

% add noise to infant and mum data
ret2(:,1) = ret(:,1)+ noiseamp*randn(length(ret),1);
ret2(:,2) = ret(:,2)+ noiseamp*randn(length(ret),1);


% visualise the mean filtered data
figure; subplot(211);
plot(ret2(1:500,1))
title('Infant arousal noisey')
set(gca, 'fontsize', 15)
subplot(212); 
plot(ret2(1:500,2))
title('Mum arousal noisey')
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)


[r, lags] = xcorr(ret2(:,1), ret2(:,2),500, 'normalized');
figure;stem(lags,r) 
xlabel('Lag in samples')
ylabel('Normalised cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation in noisey arousal')




%% remove linear trends from the data and see how this effects cross correlation value

% Part four - using matlabs detrend function - detrend the arousal data.
% Visualise the filtered arousal data. Compute the cross correlation on the
% filtered data. Visualise this result. Compare the cross correlation
% between the filtered and unfiltered and detrended data


% detrend the data
ret3 = detrend(ret,10);


% visualise the detrended data
figure; subplot(211);
plot(ret3(1:500,1))
title('Infant arousal detrended')
set(gca, 'fontsize', 15)
subplot(212); 
plot(ret3(1:500,2))
title('Mum arousal detrended')
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)


% visualise the cross correlation for detrended data
[r, lags] = xcorr(ret3(:,1), ret3(:,2),500, 'normalized');
figure;stem(lags,r) 
xlabel('Lag in samples')
ylabel('Normalised cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation in detrended arousal')




%% %% explanation on impact of linear trend



% t = 0:20;
t = 0:0.1:2;

v1 = 3*sin(t) + t + randn(length(t),1)';

v2 = v1+randn(length(v1),1)';
figure;plot(t,v1,t,v2-5,':k')
% legend('Input Data','Detrended Data','Trend','Location','northwest') 
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)
title('Raw signals')




[r, lags] = xcorr(v1,v2, 'normalized');
% note third function input is number of lags
figure;stem(lags,r)
xlabel('Lag in samples')
ylabel('Cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation using xcorr')

%% detrend


v1 = detrend(v1);
v2 = detrend(v2);
figure;plot(t,v1,t,v2-5,'k')
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)
title('Detrended signals')
set(gca,'Visible','off')


[r, lags] = xcorr(v1,v2, 'normalized');
% note third function input is number of lags
figure;stem(lags,r)
xlabel('Lag in samples')
ylabel('Cross correlation')
set(gca, 'fontsize', 15)
title('Cross correlation using xcorr')



