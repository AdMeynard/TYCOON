clear all; close all; clc;

addpath('../SynthSig');
addpath('../AlgorithmTYCOON/');
addpath(genpath('../OtherMethods'));
addpath('../PerfEvaluation/');

load('DoubleComp');
N = length(s);
t = linspace(0,(N-1)/Fs,N);

s0 = s; % Clear signal
ifi = [if1 if2].';
ami = [am1 am2].';

M = ceil(N/2);
F = zeros(M,N); % TF matrix

alpha = zeros(1,N); % Alpha vector

nbtmu = 8;
tMuVect = logspace(0,-7,nbtmu);

BigAlpha = zeros(N,nbtmu);
BigF = zeros(M,N,nbtmu);

lambda = 0.98;
gamma = 5e-4;

nbxp = 0;
stop_eps = 5e-4;

DEBUG = 0;

snrlevel = linspace(30,0,5);
for indXP = 1:length(snrlevel)
    sigma_noise = std(s0)*10^(-snrlevel(indXP)/20);
    s = s0 + sigma_noise*randn(N,1);
    
    %% TYCOON
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
    
    omega = linspace(0,Fs/2,M);
    D_tycoon(indXP) = performIND(ifi,ami,F,omega);
    
    ssynth = real(Fs/2/(M-1)*real(sum(F)));
    snrTYC(indXP) = 20*log10(std(s0)/(std(ssynth(:)-s0(:))));
    
    %% Synchrosqueezing STFT
    supp = 12;
    [tfr, tfrtic, ~, ~,~] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 101, 1, supp, 1, 0, 0, 0) ;
    supp = 24;
    [~, ~, tfrsq, ~, tfrsqtic] = ConceFT_STFT(s, 0, 0.5, 0.0005, 1, 101, 1, supp, 1, 0, 0, 0) ;
    
    omega = Fs*tfrtic ;
    omegasq = Fs*tfrsqtic ;
    
    D_stft(indXP) = performIND(ifi,ami,tfr,omega);
    D_sststft(indXP) = performIND(ifi,ami,tfrsq,omegasq);
    
     %% Synchrosqueezing CWT
    opts = struct();
    opts.motherwavelet = 'Cinfc' ;
    opts.CENTER = 1 ;
    opts.FWHM = 0.6 ;
    tt = linspace(t(1), t(end), length(t)*3) ;
    sss = interp1(t, s, tt, 'spline', 'extrap') ;
    [tfrsqC, ~, tfrsqticC] = ConceFT_CWT(tt', sss', 0, 5, 5e-3, 1, opts, 0, 0) ;
    
    F2 = tfrsqC(:, 1:3:end) ; 
    omega = tfrsqticC;
    D_sstcwt(indXP) = performIND(ifi,ami,F2,omega);

    %% EMD
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
    tvPSemd = zeros(size(F)) ;
    for kk = 1: size(F, 2)
        for iii = 1: EMDn
            if EMDif(kk,iii)/c0 > 0 
                tvPSemd(min(size(F,1),max(1,round(EMDif(kk,iii)./c0))), kk) = EMDam(kk,iii) + tvPSemd(min(size(F,1),max(1,round(EMDif(kk,iii)./c0))), kk) ; 
            end
        end
    end
    
    D_emd(indXP) = performIND(ifi,ami,tvPSemd,omega);
end

%% Figure

figure;
plot(snrlevel, D_tycoon,'-ro', snrlevel, D_stft,'--k+', snrlevel, D_sststft,'b-.', snrlevel, D_sstcwt,'g:', snrlevel, D_emd,'-cx', 'linewidth',2);
xlabel('SNR (dB)'); ylabel('Distance $D$','interpreter','latex'); grid on; set(gca,'FontSize',18);
legend('TYCOON','TFCT','SST-STFT','SST-CWT','EMD-HT')

figure;
plot(snrlevel,snrTYC,snrlevel,snrlevel,'linewidth',2)