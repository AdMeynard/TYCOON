close all; clear all;

load('ExampleAf.mat');
time = t(401:100:end-100) ;
Fs = 10 ;

am1 = smooth(cumsum(randn(length(time),1)) ./ 32, 300, 'loess') ;
am1 = 1 + (am1 + max(abs(am1))) ./ (3*max(abs(am1))) ;
phi1 = cumsum(rri(401:100:end-100))/Fs ;
if1 = resample(rri(401:100:end-100),1,1);
cf1 = (if1(2:end) - if1(1:end-1)) * Fs ;
cf1 = [cf1; cf1(end)] ;

s = am1 .* cos(2*pi*phi1);

time0 = time - time(1);
figure;
subplot('Position',[0.075 0.6 0.9 0.37]);
plot(time0,s,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$y$','interpreter','latex') ; axis tight ;

subplot('Position',[0.075 0.1 0.41 0.37]);
plot(time0,am1,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$A_1$','interpreter','latex') ; axis tight ;

subplot('Position',[0.565 0.1 0.41 0.37]);
plot(time0,if1,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$\phi''_1$','interpreter','latex') ; axis tight ;

% save('SingleComp','s','if1','cf1','am1','Fs');