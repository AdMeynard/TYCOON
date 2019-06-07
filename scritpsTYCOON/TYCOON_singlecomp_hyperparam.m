clear all; close all; clc;

addpath('../SynthSig');
addpath('../AlgorithmTYCOON/');

load('SingleComp');
N = length(s);

M = ceil(N/2);
F0 = zeros(M,N); % TF matrix

alpha0 = zeros(1,N); % Alpha vector

lambda = 0.99;
nbtmu = 3;
tMuVect = logspace(0,-6,nbtmu);
nblambda = 4;
lambdaVect = [0.95 0.99 0.999 0.9999];

gamma = 5e-4;
stop_eps = 5e-4;

BigAlpha = zeros(N,nbtmu*nblambda);
BigF = zeros(M,N,nbtmu*nblambda);

DEBUG = 0;
nbxp = 0;
for tmu = tMuVect
    for lambda = lambdaVect
        nbxp = nbxp+1;
        
       mu = tmu*lambda;
       nu = tmu*(1-lambda);

        if DEBUG
            fprintf('\n\nXP %d\n',nbxp);
        end

        [F,alpha] = tycoon(s,Fs,mu,nu,gamma,F0,alpha0,stop_eps,DEBUG);

        BigF(:,:,nbxp) = F;
        BigAlpha(:,nbxp) = alpha(:);
    end
end

%% Hyperparameters influence
t = linspace(0,N/Fs,N);
omega = linspace(0,Fs/2,M);
fmax = 3;

figure;
for nbxp = 1:nbtmu*nblambda
    subplot(nbtmu,nblambda,nbxp);
    Fn = BigF(:,:,nbxp);
    imagesc(t,omega,log1p(abs(Fn))); set(gca, 'fontsize', 14) ;
    xlabel('Time (s)');ylabel('Frequency (Hz)'); axis([0 N/Fs 0 fmax]); axis xy ; colormap(1-gray) ;
end