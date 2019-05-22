clear variables; close all; clc;

addpath('../SynthSig');
addpath('../algorithmTYCOON/');
addpath(genthpath('../OtherMethods'));

load('DoubleComp100');

Q = 100;
for qqq=1:Q
    %% signal
    s = ss{qqq};
    ifi = [ff1{qqq} ff2{qqq}].';
    ami = [aa1{qqq} aa2{qqq}].';
    
    
    %% TYCOON
    N = length(s);
    t = linspace(0,(N-1)/Fs,N);

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

        BigF(:,:,nbxp) = F;
        BigAlpha(:,nbxp) = alpha(:);
    end
    
    omega = linspace(0,Fs/2,M);
    D_tycoon(qqq) = performIND(ifi,ami,F,omega);
    %% Synchrosqueezing STFT
    [tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 151, 1, 6, 1, 0, 0, 0) ;
    omega = Fs*tfrtic ;
    omegasq = Fs*tfrsqtic ;
    
    D_stft(qqq) = performIND(ifi,ami,tfr,omega);
    D_sststft(qqq) = performIND(ifi,ami,tfrsq+eps,omegasq);
    
    %% Synchrosqueezing CWT
    opts = struct();
    opts.motherwavelet = 'Cinfc' ;
    opts.CENTER = 1 ;
    opts.FWHM = 0.2 ;
    tt = linspace(t(1), t(end), length(t)*3) ;
    sss = interp1(t, s, tt, 'spline', 'extrap') ;
    [tfrsqC, ~, tfrsqticC] = ConceFT_CWT(tt', sss', 0, 5, 5e-3, 1, opts, 0, 0) ;
    
    F2 = tfrsqC(:, 1:3:end) ; 
    omega = tfrsqticC;
    D_sstcwt(qqq) = performIND(ifi,ami,F2,omega);

    %% EMD
    EMDn = 6 ;
    allmode = eemd(s,0,1) ;
    allmode = allmode(:, 2:EMDn+1) ;

    % smooth the result if really need it.
    EMDif = zeros(size(allmode)) ; 
    EMDam = zeros(size(allmode)) ;
    for ii = 1:EMDn
        xhat = hilbert(allmode(:,ii)) ;
        EMDam(:, ii) = abs(xhat) ;
        tmp = phase(xhat) ; 
        tmpif = (tmp(2:end)-tmp(1:end-1))*Fs/2/pi ;
        EMDif(1:end-1,ii) = tmpif ; 
        EMDif(end,ii) = EMDif(end-1,ii) ;
    end
    
    omega = linspace(0,Fs/2,M);
    c0 = omega(2)- omega(1);
    tvPSemd = zeros(size(F)) ;
    for kk = 1: size(F, 2)
        for iii = 1: EMDn
            if EMDif(kk,iii)/c0 > 0 
                tvPSemd(min(size(F,1),max(1,round(EMDif(kk,iii)./c0))), kk) = EMDam(kk,iii) + tvPSemd(min(size(F,1),max(1,round(EMDif(kk,iii)./c0))), kk) ; 
            end
        end
    end
    
    D_emd(qqq) = performIND(ifi,ami,tvPSemd,omega);
end

%% Results

fprintf('TF Rep. | Mean(D) | STD(D)\n')
fprintf('Tycoon  |  %.2f   | %.2f\n',mean(D_tycoon), std(D_tycoon))
fprintf('STFT    |  %.2f   | %.2f\n',mean(D_stft), std(D_stft))
fprintf('SST-STFT|  %.2f   | %.2f\n',mean(D_sststft), std(D_sststft))
fprintf('SST-CWT |  %.2f   | %.2f\n',mean(D_sstcwt), std(D_sstcwt))
fprintf('EMD-HT  |  %.2f   | %.2f\n',mean(D_emd), std(D_emd))