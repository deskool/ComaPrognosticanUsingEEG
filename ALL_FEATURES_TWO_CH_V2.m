% AUTHOR:  Mohammad Mahdi Ghassemi, PhD Candidate, 
%          Department of Electrical Engineering and Computer Science, 
%          Massachusetts Institute of Technology
%          ghassemi@mit.edu

% PURPOSE: This file takes EEG data and generates features that
%          measure connectivity.

% INPUTS:  x  - a matrix of inputs which is [samples x channels]
%          Fs - the sampling rate of the signal
%          i  - channel number 
%          i2 - another channel number
%          featurenum - the feature you want to generate

% OUTPUTS:
%          out    - a structure containing the extracted features.
%          header - the header for the output file.

function [out, header] = ALL_FEATURES_TWO_CH( x, i,i2, featurenum, Fs )
%% FEATURES OF COMPLEXITY
header = {'COH_DEL',
          'COH_ALL',
          'PLI',
          'ACOR_MAG',
          'ACOR_LAG',
          'MI',
          'GRANGER'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THESE ARE ENTIRE EEG FEATURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE IN THE DELTA BAND - 0.36 second
if featurenum == 1
nearestpow2(length(x));
nfft = 500;
[Cxy F] = mscohere(x(:,i),x(:,i2),hanning(nfft),nfft/2,nfft,Fs);
out = mean(Cxy(find(F >= 0.5 & F<=4)));


%% COHERENCE IN THE DELTA BAND - 0.36 second
elseif featurenum == 2
nearestpow2(length(x));
nfft = 500;
[Cxy F] = mscohere(x(:,i),x(:,i2),hanning(nfft),nfft/2,nfft,Fs);
out = mean(Cxy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Phase Lag Index - 0.15 seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%http://onlinelibrary.wiley.com/doi/10.1002/hbm.20346/epdf
elseif featurenum == 3
hxi = hilbert(x(:,i));
hxj = hilbert(x(:,i2));
% calculating the INSTANTANEOUS PHASE
inst_phasei = atan(angle(hxi));
inst_phasej = atan(angle(hxj));

out = abs(mean(sign(inst_phasej - inst_phasei)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cross-Correlation - 0.15 seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif featurenum == 4
[acor,lag] = (xcorr(x(:,i),x(:,i2)));
[~, uu] = max(abs(acor));
out(4) = (max(abs(acor)) - mean(abs(acor)))/std(abs(acor));

elseif featurenum == 5
[acor,lag] = (xcorr(x(:,i),x(:,i2)));
[~, uu] = max(abs(acor));
out = lag(uu)/Fs; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mutual Information - 0.15 seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif featurenum == 6
Min = -200;
Max = 200;
step = 2;
edges = Min:step:Max;
values = mean([edges(2:end);edges(1:end-1)]);
xi = discretize(x(:,1),edges,values);
yi = discretize(x(:,2),edges,values);

out = mutual_information(xi',yi');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Granger Causality - 0.15 seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif featurenum == 7
out = granger_cause(x(:,1),x(:,2),0.05,Fs);
end

end

