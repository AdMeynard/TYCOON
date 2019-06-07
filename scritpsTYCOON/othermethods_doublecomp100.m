clear variables; close all; clc;

addpath('../SynthSig');
addpath('../AlgorithmTYCOON/');
addpath(genpath('../OtherMethods'));
addpath('../PerfEvaluation/');

load('DoubleComp100');

Q = 100;
for qqq=1:Q
    %% signal
    s = ss{qqq};
    ifi = [ff1{qqq} ff2{qqq}].';
    ami = [aa1{qqq} aa2{qqq}].';
    
    N = length(s);
    M = ceil(N/2);
    t = linspace(0,(N-1)/Fs,N);

    %% Synchrosqueezing STFT
    
    supp = 6; % control the variance of the Gaussian window
    [tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 151, 1, supp, 1, 0, 0, 0) ;
    omega = Fs*tfrtic ;
    omegasq = Fs*tfrsqtic ;
    D_stft2(qqq) = performIND(ifi,ami,tfr,omega);
    D_sststft2(qqq) = performIND(ifi,ami,tfrsq+eps,omegasq);
    
    supp = 12; % control the variance of the Gaussian window
    [tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 151, 1, supp, 1, 0, 0, 0) ;
    omega = Fs*tfrtic ;
    omegasq = Fs*tfrsqtic ;
    D_stft3(qqq) = performIND(ifi,ami,tfr,omega);
    D_sststft3(qqq) = performIND(ifi,ami,tfrsq+eps,omegasq);
    
    supp = 24; % control the variance of the Gaussian window
    [tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 151, 1, supp, 1, 0, 0, 0) ;
    omega = Fs*tfrtic ;
    omegasq = Fs*tfrsqtic ;
    D_stft4(qqq) = performIND(ifi,ami,tfr,omega);
    D_sststft4(qqq) = performIND(ifi,ami,tfrsq+eps,omegasq);
    
    supp = 48; % control the variance of the Gaussian window
    [tfr, tfrtic, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 151, 1, supp, 1, 0, 0, 0) ;
    omega = Fs*tfrtic ;
    omegasq = Fs*tfrsqtic ;
    D_stft5(qqq) = performIND(ifi,ami,tfr,omega);
    D_sststft5(qqq) = performIND(ifi,ami,tfrsq+eps,omegasq);
    
    %% Synchrosqueezing CWT
    tt = linspace(t(1), t(end), length(t)*3) ;
    sss = interp1(t, s, tt, 'spline', 'extrap') ;
    opts = struct();
    opts.motherwavelet = 'Cinfc' ;
    opts.CENTER = 1 ;
    
    % PREMIERE LARGEUR
    opts.FWHM = 0.1 ;
    [tfrsqC, ~, tfrsqticC] = ConceFT_CWT(tt', sss', 0, 5, 5e-3, 1, opts, 0, 0) ;
    F2 = tfrsqC(:, 1:3:end) ; 
    omega = tfrsqticC;
    D_sstcwt1(qqq) = performIND(ifi,ami,F2,omega);
    
    % DEUXIEME LARGEUR
    opts.FWHM = 0.2 ;
    [tfrsqC, ~, tfrsqticC] = ConceFT_CWT(tt', sss', 0, 5, 5e-3, 1, opts, 0, 0) ;
    F2 = tfrsqC(:, 1:3:end) ; 
    omega = tfrsqticC;
    D_sstcwt2(qqq) = performIND(ifi,ami,F2,omega);
    
    % TROISIEME LARGEUR
    opts.FWHM = 0.4 ;
    [tfrsqC, ~, tfrsqticC] = ConceFT_CWT(tt', sss', 0, 5, 5e-3, 1, opts, 0, 0) ;
    F2 = tfrsqC(:, 1:3:end) ; 
    omega = tfrsqticC;
    D_sstcwt3(qqq) = performIND(ifi,ami,F2,omega);
    
    % 4eme LARGEUR
    opts.FWHM = 0.6 ;
    [tfrsqC, ~, tfrsqticC] = ConceFT_CWT(tt', sss', 0, 5, 5e-3, 1, opts, 0, 0) ;
    F2 = tfrsqC(:, 1:3:end) ; 
    omega = tfrsqticC;
    D_sstcwt4(qqq) = performIND(ifi,ami,F2,omega);

    %% EMD
    EMDn = 1 ;
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
    tvPSemd = zeros(M,N) ;
    for kk = 1: N
        for iii = 1: EMDn
            if EMDif(kk,iii)/c0 > 0 
                tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) = EMDam(kk,iii) + tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) ; 
            end
        end
    end
    D_emd1(qqq) = performIND(ifi,ami,tvPSemd,omega);
    
    EMDn = 2 ;
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
    tvPSemd = zeros(M,N) ;
    for kk = 1: N
        for iii = 1: EMDn
            if EMDif(kk,iii)/c0 > 0 
                tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) = EMDam(kk,iii) + tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) ; 
            end
        end
    end
    D_emd2(qqq) = performIND(ifi,ami,tvPSemd,omega);
    
    EMDn = 4 ;
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
    tvPSemd = zeros(M,N) ;
    for kk = 1: N
        for iii = 1: EMDn
            if EMDif(kk,iii)/c0 > 0 
                tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) = EMDam(kk,iii) + tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) ; 
            end
        end
    end
    D_emd3(qqq) = performIND(ifi,ami,tvPSemd,omega);
    
    EMDn = 6;
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
    tvPSemd = zeros(M,N) ;
    for kk = 1: N
        for iii = 1: EMDn
            if EMDif(kk,iii)/c0 > 0 
                tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) = EMDam(kk,iii) + tvPSemd(min(M,max(1,round(EMDif(kk,iii)./c0))), kk) ; 
            end
        end
    end
    D_emd4(qqq) = performIND(ifi,ami,tvPSemd,omega);
end

%% Results

fprintf('TF Rep.  | Mean(D) | STD(D)\n')
fprintf('STFT1    |  %.2f   | %.2f\n',mean(D_stft2), std(D_stft2))
fprintf('STFT2    |  %.2f   | %.2f\n',mean(D_stft3), std(D_stft3))
fprintf('STFT3    |  %.2f   | %.2f\n',mean(D_stft4), std(D_stft4))
fprintf('STFT4    |  %.2f   | %.2f\n',mean(D_stft5), std(D_stft5))
fprintf('SST-STFT1|  %.2f   | %.2f\n',mean(D_sststft2), std(D_sststft2))
fprintf('SST-STFT2|  %.2f   | %.2f\n',mean(D_sststft3), std(D_sststft3))
fprintf('SST-STFT3|  %.2f   | %.2f\n',mean(D_sststft4), std(D_sststft4))
fprintf('SST-STFT4|  %.2f   | %.2f\n',mean(D_sststft5), std(D_sststft5))
fprintf('SST-CWT1 |  %.2f   | %.2f\n',mean(D_sstcwt1), std(D_sstcwt1))
fprintf('SST-CWT2 |  %.2f   | %.2f\n',mean(D_sstcwt2), std(D_sstcwt2))
fprintf('SST-CWT3 |  %.2f   | %.2f\n',mean(D_sstcwt3), std(D_sstcwt3))
fprintf('SST-CWT4 |  %.2f   | %.2f\n',mean(D_sstcwt4), std(D_sstcwt4))
fprintf('EMD-HT1  |  %.2f   | %.2f\n',mean(D_emd1), std(D_emd1))
fprintf('EMD-HT2  |  %.2f   | %.2f\n',mean(D_emd2), std(D_emd2))
fprintf('EMD-HT3  |  %.2f   | %.2f\n',mean(D_emd3), std(D_emd3))
fprintf('EMD-HT4  |  %.2f   | %.2f\n',mean(D_emd4), std(D_emd4))