%% 

load sampleEEGdata.mat;

% simulation parameters 
f = 7;
srate = 256;
t  = -1:1/srate:1.5;
ntrials = 100;
d = zeros( length(t), ntrials);
s = zeros(length(t), ntrials);
noiseamp=0.2;

pnts = length(t);


for i = 1:ntrials

r = 2*pi*rand(1);

s(1:100,i) = sin(2*pi*f*t(1:100)+r); % constant phase

s(101:340,i) =  sin(2*pi*f*t(101:340)  );
% 
s(341:641,i) = sin(2*pi*f*t(341:641)+r); % constant phase



d(1:100,i) = sin(2*pi*f*t(1:100)+r); % constant phase

d(101:340,i) = sin(2*pi*f*t(101:340 ) ); % 

d(341:641,i) = sin(2*pi*f*t(341:641)+r); % constant phase

end


signal1 = s + noiseamp * randn(length(t),ntrials);

signal2 = d + noiseamp * randn(length(t),ntrials);

% figure(3)
% set(gcf,'color','w');


% extract angles from Hilbert transform
angles1 = angle(hilbert( signal1 ));
angles2 = angle(hilbert( signal2 ));

% show phase angle time series
subplot(223)
plot(EEG.times,angles1(1:end-1,1:10),'k')
hold on
plot(EEG.times,angles2(1:end-1,1:10)-10,'k')
% ylabel('Phase angle')
% title('Phase angle over time')
xlim([-1000 800])
set(gca, 'ytick', [],'xtick', [], 'fontsize', 30)


%% 


% simulation parameters
ntrials = 100;
f = 7;
srate = 256;
t  = -1:1/srate:1.5;
d = zeros( length(t), ntrials);
s = zeros(length(t), ntrials);
noiseamp=0.2;
pnts = length(t);


for i = 1:ntrials

r = 2*pi*rand(1);

s(1:100,i) = sin(2*pi*f*t(1:100)+r); % constant phase

s(101:340,i) =  sin(2*pi*f*t(101:340)  );
% 
s(341:640,i) = sin(2*pi*f*t(341:640)+r); % constant phase



d(1:200,i) = sin(2*pi*f*t(1:200)+r); % constant phase

d(201:440,i) = sin(2*pi*f*t(201:440 ) ); % 

d(441:640,i) = sin(2*pi*f*t(441:640)+r); % constant phase


end


data1 = s + noiseamp * randn(length(t),ntrials);
data2 = d + noiseamp * randn(length(t),ntrials);

angles1 = angle(hilbert( data1 ));
angles2 = angle(hilbert( data2 ));

% show phase angle time series
subplot(224)
% figure
plot(EEG.times,angles1(1:end-1,1:10),'k')
hold on
plot(EEG.times,angles2(1:end-1,1:10)-10,'k')
% ylabel('Phase angle')
% title('Phase angle over time')
xlim([-1000 1000])
set(gca, 'ytick', [],'xtick', [], 'fontsize', 30)











