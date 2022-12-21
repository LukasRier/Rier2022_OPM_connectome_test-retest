function [po, fs] = get_PSD(data,f)
%% [po, fs] = get_PSD(data,f)
% Calculate power spectrum from data (Number of channels x number of
% samples) with sampling frequency f
% po....power spectrum for each channel
% fs....frequency axis values
%

S = []; 
eD = [];
S.trialength = 5*f;
for n = 1:floor(length(data)/S.trialength)
    eD(:,:,n) = data(:,((n-1)*S.trialength) + 1:n*S.trialength);
end
% set flattop window 
nepochs = size(eD,3);
F = [0:.1:600]';
pow = zeros(size(F,1),size(eD,1),nepochs);
wind = window(@flattopwin,size(eD,2));
% loop over windows dc correcting along the way
for j = 1:nepochs
    Btemp = eD(:,:,j)';
    mu = mean(Btemp);
    zf = bsxfun(@minus,Btemp,mu);
    nchan = size(zf,2);
    fzf = zf;
    [pxx,fs] = periodogram(fzf,wind,F,f);
    pow(:,:,j) = sqrt(pxx);
end
Nchans = size(data,1);
chans = [1:Nchans];
po = median(pow(:,chans,:),3);