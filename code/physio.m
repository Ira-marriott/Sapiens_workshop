% FromDir = '/Users/iramarriottharesign/Desktop/home awake synched_orig/';
FromDir = 'C:\Users\Ira Marriott\Desktop\Sapiens workshop\home awake synched_orig\';



cd(FromDir)
files = dir('*.mat');

%% 
load(strcat(FromDir, files(2).name))

% use dataset two for example here.


% 1 is good example of negative correlation
% 2 has high correlation 


% plot(ret(:,1))
% hold on 
% plot(ret(:,2))

% 
% [r, lags] = xcorr(ret(:,1), 2000,'normalized');
% figure;stem(lags,r) 




x = isnan(ret(:,1));
toremove = find(x==1);


ret(toremove,:) = [];



figure; subplot(211);plot(ret(1:500,1))
set(gca, 'XTick',[], 'YTick',[])

% hold on 
subplot(212); plot(ret(1:500,2))
% set(gca, 'XTick',[], 'YTick',[])
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)

[r, lags] = xcorr(ret(:,1), ret(:,2),500, 'normalized');
figure;stem(lags,r) 


%% 


ret2 = medfilt1(ret,20);

figure; subplot(211);plot(ret2(1:500,1))
set(gca, 'XTick',[], 'YTick',[])

% hold on 
subplot(212); plot(ret2(1:500,2) )
% set(gca, 'XTick',[], 'YTick',[])
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)

[r, lags] = xcorr(ret2(:,1), ret2(:,2),500, 'normalized');
figure;stem(lags,r) 


%% 

ret3 = detrend(ret,10);

figure; subplot(211);plot(ret3(1:500,1))
set(gca, 'XTick',[], 'YTick',[])

% hold on 
subplot(212); plot(ret3(1:500,2) )
% set(gca, 'XTick',[], 'YTick',[])
xlabel('Samples')
ylabel('Amplitude')
set(gca, 'fontsize', 15)

[r, lags] = xcorr(ret3(:,1), ret3(:,2),500, 'normalized');
figure;stem(lags,r) 


%% 



% figure; plot(ret(1:300,1), 'linew', 2)
% set(gca, 'XTick',[], 'YTick',[])
% 
% % hold on 
% figure; plot(ret(1:300,2),  'linew', 2)
% set(gca, 'XTick',[], 'YTick',[])
% xlabel('Samples')
% ylabel('Amplitude')
% set(gca, 'fontsize', 20)

% 


r2 = corrcoef(ret(:,1), ret(:,2))


%% 


x=1;
for i=2:30
    x(i) = 1.1*x(i-1) + randn;
end

subplot(121)
plot(x,'kp-','linewid',1,'markerface','g','markersize',10)
set(gca,'xlim',[0 31])
title('Univariate autoregression')
legend('x_t = 1.1*x_t_-_1 + randn')
% title('Non-stationary autoregressive process')
set(gca, 'fontsize',20)
set(gca, 'XTick',[], 'YTick',[])



x=[1 1.5];
for i=3:30
    x(i) = 1.2*x(i-2) + -.3*x(i-1) + randn;
end

subplot(122)
plot(x,'kp-','linewid',1,'markerface','g','markersize',10)
set(gca,'xlim',[0 31])
legend('x_t = -0.3*x_t_-_1 + 1.2*x_t_-_2 + randn')
% title('Non-stationary autoregressive process')

set(gca, 'XTick',[], 'YTick',[])
set(gca, 'fontsize',20)



%% 


x = randn;
for i=2:30
    x(i) = exp(cos(pi*x(i-1))) + randn;
end

subplot(211)
plot(x,'mo-','linewid',1,'markerface','k','markersize',6)
set(gca,'xlim',[0 31])
% legend('x_t = e^c^o^s^(^\pi^xt-1^) + randn')
% title('Stationary autoregressive process')
set(gca, 'XTick',[], 'YTick',[])


x=randn(2,1);
for i=3:30
    x(i) = .2*x(i-1) - .4*x(i-2) + randn;
end

subplot(212)
plot(x,'mo-','linewid',1,'markerface','k','markersize',6)
set(gca,'xlim',[0 31])
% legend('x_t = .2x_t_-_1 - .4x_t_-_2 + randn')
% title('Stationary autoregressive process')

set(gca, 'XTick',[], 'YTick',[])








