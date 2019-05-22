function [tfr, tfrtic, tfrsq, Cm, ConceFTmean, ConceFTmedian, tfrsqtic] = Cmedian(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT) ;

%
% Usage: 
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
%
% MT = 1: ordinary SST; MT > 1: ConceFT
% alpha: resolution in the frequency axis
% WinLen, dim, supp: for hermf.m
%
% Example:
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);


N = length(x) ;

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       	%% generate the window for short time Fourier transform (STFT)
[h, Dh, ~] = hermf(WinLen, dim, supp) ;


%=======================================

[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', 0, 0) ;


%=======================================

if MT > 1

	%% Conceft

	[~, ~, Cm, tfrsqtic] = Cmedianbase(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT);

end




%========================

ConceFTmean = zeros(size(tfrsq)) ;
ConceFTmedian0 = zeros([size(tfrsq) MT]) ;

if MT > 1

    %% Conceft

    for ii = 1: MT
        fprintf('\b\b\b\b') ;   tmp = sprintf('%4d',ii) ; fprintf([tmp]) ;
        rv = randn(1, dim) ; rv = rv ./ norm(rv) ;
        rh = rv * h ;
        rDh = rv * Dh ;

        [~, ~, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, rh', rDh', 0, 0);

        ConceFTmean = ConceFTmean + tfrsq ;
        ConceFTmedian0(:,:,ii) = abs(tfrsq) ;
    end

    ConceFTmean = ConceFTmean ./ MT ;

	ConceFTmedian = [] ;
	for ii = 1: size(ConceFTmedian0,1)
		for jj = 1:size(ConceFTmedian0,2)
			ConceFTmedian(ii, jj) = median(ConceFTmedian0(ii,jj,:)) ;
		end
	end

    fprintf('\n') ;

end




end
