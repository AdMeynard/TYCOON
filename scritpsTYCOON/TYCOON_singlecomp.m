clear all; close all; clc;

addpath('../SynthSig');
addpath('../algorithmTYCOON/');

load('SingleComp');
N = length(s);

M = ceil(N/2);
F = zeros(M,N); % TF matrix

alpha = zeros(1,N); % Alpha vector

nbtmu = 10;
tMuVect = logspace(1,-10,nbtmu); % on teste differents hyperparam lambda

BigAlpha = zeros(N,nbtmu);
BigF = zeros(M,N,nbtmu);

lambda = 0.99;
gamma = 5e-4;

nbxp = 0;
stop_eps = 5e-4;

DEBUG = 0;
for tmu = tMuVect
    nbxp = nbxp+1;
    Fold = F;
    alphaold = alpha;
    mu = tmu*lambda;
    nu = tmu*(1-lambda);
    alpha = zeros(1,N); % Alpha vector

    if DEBUG
        fprintf('\n\nXP %d\n',nbxp);
    end

    [F,alpha] = tycoon(s,Fs,mu,nu,gamma,F,alpha,stop_eps,DEBUG);

    if DEBUG
        t = linspace(0,N/Fs,N);
        figure(1)
        plot(t,s,'b',t,Fs/2/(M-1)*real(sum(F)),'r');

        figure(2)
        plot(t,alpha) ; set(gca, 'fontsize', 24) ;
        xlabel('Time (sec)') ; ylabel('alpha') ; axis tight ;

        figure(3);
        omega = linspace(0,Fs/2,M);
        imagesc(t,omega,log1p(abs(F))); set(gca, 'fontsize', 18) ;
        xlabel('Time (sec)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 Fs/2]); axis xy ; colormap(1-gray) ;
    end
    drawnow;

    BigF(:,:,nbxp) = F;
    BigAlpha(:,nbxp) = alpha(:);
end


%% Analysis

t = linspace(0,N/Fs,N);
fmax = 3;
figure(1)
plot(t,s,'b',t,Fs/2/(M-1)*real(sum(F)),'r');

figure;
subplot('Position',[0.085 0.1 0.9 0.32]);
plot(t,cf1/norm(cf1,2),'b',t,alpha/norm(alpha,2),'r','linewidth',2) ; set(gca, 'fontsize', 24) ; 
legend('$\phi''''_1$','Estimation','interpreter','latex')
xlabel('Time (s)') ; ylabel('Chirp factor') ; axis tight ;

subplot('Position',[0.085 0.52 0.41 0.45]);
omega = linspace(0,Fs/2,M);
imagesc(t,omega,log1p(abs(F))); set(gca, 'fontsize', 18) ;
xlabel('Time (s)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 fmax]); axis xy ; colormap(1-gray) ;

subplot('Position',[0.575 0.52 0.41 0.45]);
omega = linspace(0,Fs/2,M);
imagesc(t,omega,log1p(abs(F)));
hold on; plot(t,if1,'r','linewidth',2) ;
set(gca, 'fontsize', 18) ;
xlabel('Time (s)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 fmax]); axis xy ; colormap(1-gray) ;

%% Comparison
addpath(genpath('OtherMethods'))

[tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 101, 1, 6, 1, 0, 0, 0) ;

figure;
subplot('Position',[0.085 0.15 0.41 0.75]);
omega = linspace(0,Fs/2,M);
imagesc(t,omega,log1p(abs(tfr))); set(gca, 'fontsize', 18) ;
xlabel('Time (s)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 fmax]); axis xy ; colormap(1-gray) ;

subplot('Position',[0.575 0.15 0.41 0.75]);
omega = linspace(0,Fs/2,M);
imagesc(t,omega,log1p(abs(tfrsq)));
set(gca, 'fontsize', 18) ;
xlabel('Time (s)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 fmax]); axis xy ; colormap(1-gray) ;