close all; clear all;
addpath(genpath('../OtherMethods/'));
load('ExampleAf.mat');
time = t(401:100:end-100) ;
Fs = 10 ;

for qqq=1:100
    am1 = smooth(cumsum(randn(length(time),1)) ./ 32, 300, 'loess') ;
    am1 = 1 + (am1 + max(abs(am1))) ./ (3*max(abs(am1))) ;
    phi1 = cumsum(rri(401:100:end-100))/Fs ;
    if1 = resample(rri(401:100:end-100),1,1);
    cf1 = (if1(2:end) - if1(1:end-1)) * Fs ;
    cf1 = [cf1; cf1(end)] ;

    am2 = smooth(cumsum(randn(length(am1),1)) ./ Fs, 300, 'loess') ;
    am2 = 1 + (am2 + max(abs(am2))) ./ (3*max(abs(am2))) ;
    if2c = smooth(cumsum(randn(length(am1),1)) ./ Fs, 20, 'loess') ;
    if2c = if2c - mean(if2c) + 1.8956 ;
    if2 = 1 + 2*(if2c + 0.5*max(abs(if2c))) ./ ( 1*max(abs(if2c)) ) ;
    phi2 = cumsum(if2) / Fs ;

    mask2 = ones(length(am1), 1) ; mask2(1:202) = 0 ;
    if2 = if2 - sin(time) ;
    if2(1:202) = nan ;

    s1 = am1 .* cos(2*pi*phi1) ;
    s2 = am2 .* cos(2*pi*(phi2+cos(time))) .* mask2 ;
    s_nonoise = s1 + s2 ;

    s = s1 + s2;
    time0 = time - time(1);
    
    ss{qqq} = s;
    aa1{qqq} = am1;
    aa2{qqq} = am2;
    ff1{qqq} = if1;
    ff2{qqq} = if2;
end

figure;
subplot('Position',[0.075 0.6 0.9 0.37]);
plot(time0,s1,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$y$','interpreter','latex') ; axis tight ;

subplot('Position',[0.075 0.1 0.41 0.37]);
plot(time0,am1,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$A_1$','interpreter','latex') ; axis tight ;

subplot('Position',[0.565 0.1 0.41 0.37]);
plot(time0,if1,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$\phi''_1$','interpreter','latex') ; axis tight ;

figure;
subplot('Position',[0.075 0.6 0.9 0.37]);
plot(time0,s2,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$y$','interpreter','latex') ; axis tight ;

subplot('Position',[0.075 0.1 0.41 0.37]);
plot(time0,am2,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$A_2$','interpreter','latex') ; axis tight ;

subplot('Position',[0.565 0.1 0.41 0.37]);
plot(time0,if2,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$\phi''_2$','interpreter','latex') ; axis tight ;


figure;
subplot('Position',[0.075 0.6 0.9 0.37]);
plot(time0,s,'linewidth',2); set(gca, 'fontsize', 24) ;
xlabel('Time (s)') ; ylabel('$y$','interpreter','latex') ; axis tight ;

supp = 12;
N = length(time0);
[tfr, ~, ~, ~,~] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 101, 1, supp, 1, 0, 0, 0) ;
subplot('Position',[0.075 0.1 0.9 0.37]);
M = ceil(N/2);
omega = linspace(0,Fs/2,M);
imagesc(time0,omega,log1p(abs(tfr))); set(gca, 'fontsize', 18) ;
xlabel('Time (s)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 Fs/2]); axis xy ; colormap(1-gray) ;