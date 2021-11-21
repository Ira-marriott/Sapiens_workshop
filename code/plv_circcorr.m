

load sampleEEGdata.mat

data  = EEG.data;
freqs = linspace(2,30,40);
fs =256;

[~, angles1,~]=conv_mf2(data, fs, freqs);


n = length(data);


r = sum(exp(1i*angles1),4);

% obtain mean by
temp = angle(r);


%% 

a = squeeze(angles1(1,8,:));
b = squeeze(angles1(2,8,:));


figure;plot(a); hold on; plot(b); xlim([1 500])


ji = 1:fs:length(a);

r = zeros(1, length(ji));

for i = 1:length(ji)-1
    
    indx = ji(i);
    
    r(i) = circular_corr( a(indx:indx+fs/2), b(indx:indx+fs/2) );
    
end

figure;plot(r)


%% 
figure;plot(a(1:500)); hold on; plot(b(1:500)); 

% 
a = squeeze(angles1(1,8,:));
b = squeeze(angles1(2,8,:));


ji = 1:fs:length(a);

plv = zeros(1, length(ji));

for i = 1:length(ji)-1
    
    indx = ji(i);
    
    plv(i) = abs(mean(exp(1i* (a(indx:indx+fs/2) - b(indx:indx+fs/2)) ),1));
    
end


figure;plot(plv)
hold on
plot(r)
legend({'plv' , 'circcorr'})















