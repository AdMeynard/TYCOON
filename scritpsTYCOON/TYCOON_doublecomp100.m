clear variables; close all; clc;

addpath('../SynthSig');
addpath('../AlgorithmTYCOON/');
addpath(genpath('../OtherMethods'));
addpath('../PerfEvaluation/');

load('DoubleComp100');

lambdaVect = [0.95 0.98 0.999];
Q = 100;
for qqq=1:Q
    %% signal
    s = ss{qqq};
    ifi = [ff1{qqq} ff2{qqq}].';
    ami = [aa1{qqq} aa2{qqq}].';
    
    N = length(s);
    t = linspace(0,(N-1)/Fs,N);
    M = ceil(N/2);
    
    
    %% TYCOON
    hypn = 0;
    for lambda = lambdaVect
        hypn = hypn+1;
        F = zeros(M,N); % TF matrix
        alpha = zeros(1,N); % Alpha vector

        nbtmu = 8;
        tMuVect = logspace(0,-7,nbtmu);

        BigAlpha = zeros(N,nbtmu);
        BigF = zeros(M,N,nbtmu);

        gamma = 5e-4;

        nbxp = 0;
        stop_eps = 5e-4;

        for tmu = tMuVect
            nbxp = nbxp+1;
            Fold = F;
            alphaold = alpha;
            mu = tmu*lambda;
            nu = tmu*(1-lambda);
            alpha = zeros(1,N); % Alpha vector
            [F,alpha] = tycoon(s,Fs,mu,nu,gamma,F,alpha,stop_eps,0);
        end
        
        omega = linspace(0,Fs/2,M);
        D_tycoon(qqq,hypn) = performIND(ifi,ami,F,omega);
        
    end
    
end

%% Results

fprintf('    TF Rep.          | Mean(D) | STD(D)\n')
for n = 1:hypn
    fprintf('Tycoon lambda=%.3f  |  %.2f   | %.2f\n', lambdaVect(n), mean(D_tycoon(:,n)), std(D_tycoon(:,n)))
end