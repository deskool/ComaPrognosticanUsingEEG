% AUTHOR:  Mohammad Mahdi Ghassemi, PhD Candidate, 
%          Department of Electrical Engineering and Computer Science, 
%          Massachusetts Institute of Technology
%          ghassemi@mit.edu

% PURPOSE: This file takes EEG data and generates features that
%          measure complexity and category.

% INPUTS:  x  - a matrix of inputs which is [samples x channels]
%          Fs - the sampling rate of the signal
%          i - channel number
%          featurenum - the feature you want to generate

% OUTPUTS:
%          out - a structure containing the extracted features.

% EXAMPLE: Fs = 100; x = randi([-200 200],500,22);

function [out, header] = ALL_FEATURES_ONE_CH( x, i, featurenum, Fs )
%% IN THIS CASE X is a matrix of the EEGs

warning('off','all')
%% FEATURES
header = {'Shannon'
    'Tsalis_q1'
    'Tsalis_q2'
    'Tsalis_q3'
    'Tsalis_q4'
    'Tsalis_q5'
    'Tsalis_q6'
    'Tsalis_q7'
    'Tsalis_q8'
    'Tsalis_q9'
    'Tsalis_q10'
    'subband_IQ'
    'Cepstrum_1'
    'Cepstrum_2'
    'Lyapunov'
    'Fractal_Dim'
    'Hjorth_mobility'
    'Hjorth_complexity'
    'false_nearest_neighbor'
    'ARMA_1'
    'ARMA_2'
    'MedianFrequency'
    'StandardDev'
    'AlphaDeltaRatio'
    'Regularity'
    'DeltaPow'
    'ThetaPow'
    'AlphaPow'
    'BetaPow'
    'GammaPow'
    'MuPow'
    'LOW_VOLTAGE_5'
    'LOW_VOLTAGE_10'
    'LOW_VOLTAGE_20'
    'NORMAL'
    'DIFF_SLOW'
    'EPI_NUM_SPIKE'
    'EPI_PEAK_AFTER_DELTA'
    'EPI_NUM_SHARP_WAVE'
    'burst_length_mean'
    'burst_length_std'
    'supression_length_mean'
    'supression_length_std'
    'burst_supression_total'
    'num_bursts'
    'num_sup'
    'burst_delta'
    'burst_theta'
    'burst_alpha'
    'burst_beta'
    'burst_gamma'
    'burst_mu'
    };

%% COMPUTE THE SHANNON ENTROPY
if featurenum == 1
    bin_min = -200; bin_max = 200; binWidth = 2;
    [out, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    
    %% COMPUTE THE TSALIS ENTROPY
elseif featurenum == 2
    bin_min = -200; bin_max = 200; binWidth = 2; q = 1
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
    
elseif featurenum == 3
    bin_min = -200; bin_max = 200; binWidth = 2; q = 2
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 4
    bin_min = -200; bin_max = 200; binWidth = 2; q = 3
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 5
    bin_min = -200; bin_max = 200; binWidth = 2; q = 4
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 6
    bin_min = -200; bin_max = 200; binWidth = 2; q = 5
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 7
    bin_min = -200; bin_max = 200; binWidth = 2; q = 6
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 8
    bin_min = -200; bin_max = 200; binWidth = 2; q = 7
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 9
    bin_min = -200; bin_max = 200; binWidth = 2; q = 8
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 10
    bin_min = -200; bin_max = 200; binWidth = 2; q = 9
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
elseif featurenum == 11
    bin_min = -200; bin_max = 200; binWidth = 2; q = 10
    [~, prob] = CRI_ShannonEntropy(x(:,i), bin_min, bin_max, binWidth);
    out = Tsallis_entro(prob',(q + 0.01));
    
    
    %% Sub-Band Information Quantity
elseif featurenum == 12
    bin_min = -200;
    bin_max = 200;
    binWidth = 2;
    
    waveletFunction = 'db8';
    [C,L] = wavedec(x(:,1),5,waveletFunction);
    % Calculation The Coificients Vectors
    cGamma = detcoef(C,L,1); %GAMMA 50
    cBeta = detcoef(C,L,2);  %BETA 25
    cAlpha = detcoef(C,L,3); %ALPHA 12.5
    cTheta = detcoef(C,L,4); %THETA 6.25
    cDelta = detcoef(C,L,5); %DELTA  3.125
    %Step 3:Compute the entropy of each of the sub-bands.
    [sbe(1), prob] = CRI_ShannonEntropy(cGamma, bin_min, bin_max, binWidth);
    [sbe(2), prob] = CRI_ShannonEntropy(cBeta, bin_min, bin_max, binWidth);
    [sbe(3), prob] = CRI_ShannonEntropy(cAlpha, bin_min, bin_max, binWidth);
    [sbe(4), prob] = CRI_ShannonEntropy(cTheta, bin_min, bin_max, binWidth);
    [sbe(5), prob] = CRI_ShannonEntropy(cDelta, bin_min, bin_max, binWidth);
    
    %Sub-band information quantity
    out = mean(sbe)
    
    %%  Cepstrum coefficients - 0.28 seconds
    %The cepstrum can be seen as information about rate of change in
    %the different spectrum bands. The cepstrum coefficients have been
    %shown to be
    
elseif featurenum == 13
    % The number of components
    num_lcp = 2;
    
    %Take the autocorrelation
    hac = dsp.Autocorrelator;
    hac.MaximumLagSource = 'Property';
    hac.MaximumLag = num_lcp;           % Compute autocorrelation lags between [0:9]
    a = step(hac, x(:,i));
    
    % Run levinson solver to find LPC coefficeints.
    hlevinson = dsp.LevinsonSolver;
    hlevinson.AOutputPort = true;   % Output polynomial coefficients
    LPC = step(hlevinson, a);       % Compute LPC coefficients
    
    %Now convert to cepstral coefficents.
    hlpc2cc = dsp.LPCToCepstral('CepstrumLength',round(1.5*num_lcp));
    
    CC{i} = step(hlpc2cc, LPC); % Convert LPC to CC.
    out = CC{i}(2);
    
elseif featurenum == 14
    num_lcp = 2;
    
    %Take the autocorrelation
    hac = dsp.Autocorrelator;
    hac.MaximumLagSource = 'Property';
    hac.MaximumLag = num_lcp;           % Compute autocorrelation lags between [0:9]
    a = step(hac, x(:,i));
    
    % Run levinson solver to find LPC coefficeints.
    hlevinson = dsp.LevinsonSolver;
    hlevinson.AOutputPort = true;   % Output polynomial coefficients
    LPC = step(hlevinson, a);       % Compute LPC coefficients
    
    %Now convert to cepstral coefficents.
    hlpc2cc = dsp.LPCToCepstral('CepstrumLength',round(1.5*num_lcp));
    
    CC{i} = step(hlpc2cc, LPC); % Convert LPC to CC.
    out = CC{i}(3);
    
    
    
    
    %% Lyapunov Exponents - 43 secs
elseif featurenum == 15
    out = lyapunov(x(:,1),Fs);
    
    %% Higuchi Fractal Dimention
    %FRactal dimention is a non-local measure that describes the complexity of
    %the fundamental patterns in the signal. It is also understood as a measure
    %of self-similarity of the signal. That is, if we zoom in on a very small
    %portion of the signal, are the properties the same as the larger portion
    %of the signal?
elseif featurenum == 16
    kmax = 10;
    out = Higuchi_FD(x(:,i), kmax);
    
    %% Hjorth Parameters of mobility and complexity
elseif featurenum == 17
    [out] = HjorthParameters(x(:,i));
    
elseif featurenum == 18
    [~, out] = HjorthParameters(x(:,i));
    
    %% False Nearest Neighbors
    % Average Mutual Information
    %There exist good arguments that if the time delayed mutual information exhibits a marked minimum at a certain value of tex2html_wrap_inline6553, then this is a good candidate for a reasonable time delay.
elseif featurenum == 19
    npts = 1000;
    maxdims = 50;
    max_delay = 200;
    distance_thresh = 0.5;
    
    tt = ami(x(:,i),npts,max_delay);
    [~,idx] = min(tt);
    nn = false_neighbors_kd(x(:,i), idx, maxdims, npts, 1, 0, 0, 0);
    mindim = find((mean(nn.z./nn.d > distance_thresh) < 0.1) == 1);
    out = mindim(1);
    
    %% Autoregressive Models on each segment - 8.87
    %http://www.sciencedirect.com/science/article/pii/S1388245700003795
    
elseif featurenum == 20
    mod{i} = arima(2,0,0);
    arma_mod{i} = estimate(mod{i},x(:,i));
    out = arma_mod{i}.AR{1};
    
elseif featurenum == 21
    mod{i} = arima(2,0,0);
    arma_mod{i} = estimate(mod{i},x(:,i));
    out = arma_mod{i}.AR{2};
    

    %% Median Frequency of the Signal
elseif featurenum == 22
    [out] = MedianFrequency(x(:,i),1,50)
    
    %% COMPUTE THE STANDARD DEVIAITON OF THE SIGNAL
elseif featurenum == 23
    out = std(x(:,i));
    
    %% ALPHA TO DELTA RATIO.  - 0.35
     
elseif featurenum == 24
    window=Fs*2;
    noverlap=32;
    h = spectrum.welch('Hamming',window);
    
    hpsd = psd(h,x(:,i),'Fs',Fs);
    Pw = hpsd.Data;
    Fw = hpsd.Frequencies;
    out = bandpower(x(:,i),Fs,[8,13])/bandpower(x(:,i),Fs,[0.5,4]);
    
    %% REGULARITY - 0.02
    %In this technique, we first squared the signal and applied a moving-average
    %filter with a window of 0.5 seconds to create a nonnega- tive smooth signal.
    %The window length of the moving average was set at 0.5 seconds.
    
    %x = [zeros(990,1);ones(10,1);zeros(990,1);ones(10,1);]
    %x = [ones(10000,1)]
elseif featurenum == 25
    %square the signal
    in_x = x(:,i).^2;
    %find the filter length in samples - we want 0.5 seconds.
    num_wts = Fs/2;
    
    a = 1;
    wts = ones(1,num_wts)/num_wts;
    q = filter(wts,a,in_x);
    %Subsequently, we sorted the values of the smoothed signal in ?descending? order
    q = sort(q,'descend');
    N = length(q);
    u = 1:N;
    %COMPUTE THE REG
    %var(q)
    out = sqrt(sum(u.^2.* q')/(sum(q)*N^2*1/3));
    
    
    
    %% BAND POWERS -  0.03
elseif featurenum == 26
    out = bandpower(x(:,i),Fs,[0.5,4]);
elseif featurenum == 27
    out = bandpower(x(:,i),Fs,[4,7]);
elseif featurenum == 28
    out = bandpower(x(:,i),Fs,[8,15]);
elseif featurenum == 29
    out = bandpower(x(:,i),Fs,[16,31]);
elseif featurenum == 30
    out = bandpower(x(:,i),Fs,[32,49]);
elseif featurenum == 31
    out = bandpower(x(:,i),Fs,[8,12]);
    
    %% EEG STATES
    
    %% FIND LOW VOLTAGE - 0.08 seconds
elseif featurenum == 32
    LOW_VOLTAGE_5 = mean(abs(x(:,i)) < 5);
    out = LOW_VOLTAGE_5;
elseif featurenum == 33
    LOW_VOLTAGE_10 = mean(abs(x(:,i)) < 10);
    out = LOW_VOLTAGE_10;
elseif featurenum == 34
    LOW_VOLTAGE_20 = mean(abs(x(:,i)) < 20);
    out = LOW_VOLTAGE_20;
    
    %% FIND DIFFUSED SLOWLY < 8Hz peak OR NORMAL >= 8Hz peak  - 2.52 seconds
elseif featurenum == 35
    for j = 1:(Fs/2 - 1)
        pow(j) = bandpower(x(:,i),Fs,[j-1 j]);
    end
    [a ind] = max(pow);
    
    NORMAL = ind >=8;
    out = NORMAL;
    
elseif featurenum == 36
    for j = 1:(Fs/2 - 1)
        pow(j) = bandpower(x(:,i),Fs,[j-1 j]);
    end
    [a ind] = max(pow);
    
    DIFF_SLOW = ind < 8;
    out = DIFF_SLOW;
    
    %% Epileptiform - 0.39 Seconds
    %You’re basically looking for high amplitude (usually >100 uv) sharp transients.
    %So we will count the number of peaks that are 3 std from mean, and 70 ms
    %in width or less. We are also looking for slow waves following the spike.
elseif featurenum == 37
    stds_away = 3;
    dealta_jumps = 2;
    %remove the mean from the signal and look for the spikes
    %Find a spikes lasting lt 70 ms -14 samples
    [pks, locs]=findpeaks(x(:,i)-mean(x(:,i)),'MaxPeakWidth',14,'SortStr', 'descend');
    
    %make sure that it is more than 4 standard deviations from the median
    pk_index = locs(find(stds_away * std(x(:,i)) < pks));
    NUM_SPIKE = length(pk_index);
    
    %check for a stronger delta in the one second following the shark peak.
    try
        for j = 1:length(locs)
            %is the power in the preceeding second less than in the seconds before the spike.
            delta_jump = bandpower(x(locs(j)-Fs:locs(j),i),Fs,[0.5 4])*dealta_jumps < bandpower(x(locs(j):locs(j)+Fs,i),Fs,[0.5 4]);
        end
        PEAK_AFTER_DELTA = length(delta_jump);
    catch
        PEAK_AFTER_DELTA = 0;
    end
    
    [pks, locs]=findpeaks(x(:,i)-mean(x(:,i)),'MinPeakWidth',14, 'MaxPeakWidth', 40, 'SortStr', 'descend');
    pk_index = locs(find(stds_away * std(x(:,i)) < pks));
    NUM_SHARP_WAVE = length(pk_index);
    out = NUM_SPIKE;
    
elseif featurenum == 38
    stds_away = 3;
    dealta_jumps = 2;
    %remove the mean from the signal and look for the spikes
    %Find a spikes lasting lt 70 ms -14 samples
    [pks, locs]=findpeaks(x(:,i)-mean(x(:,i)),'MaxPeakWidth',14,'SortStr', 'descend');
    
    %make sure that it is more than 4 standard deviations from the median
    pk_index = locs(find(stds_away * std(x(:,i)) < pks));
    NUM_SPIKE = length(pk_index);
    
    %check for a stronger delta in the one second following the shark peak.
    try
        for j = 1:length(locs)
            %is the power in the preceeding second less than in the seconds before the spike.
            delta_jump = bandpower(x(locs(j)-Fs:locs(j),i),Fs,[0.5 4])*dealta_jumps < bandpower(x(locs(j):locs(j)+Fs,i),Fs,[0.5 4]);
        end
        PEAK_AFTER_DELTA = length(delta_jump);
    catch
        PEAK_AFTER_DELTA = 0;
    end
    
    [pks, locs]=findpeaks(x(:,i)-mean(x(:,i)),'MinPeakWidth',14, 'MaxPeakWidth', 40, 'SortStr', 'descend');
    pk_index = locs(find(stds_away * std(x(:,i)) < pks));
    NUM_SHARP_WAVE = length(pk_index);
    out = PEAK_AFTER_DELTA;
    
    
    
elseif featurenum == 39
    stds_away = 3;
    dealta_jumps = 2;
    %remove the mean from the signal and look for the spikes
    %Find a spikes lasting lt 70 ms -14 samples
    [pks, locs]=findpeaks(x(:,i)-mean(x(:,i)),'MaxPeakWidth',14,'SortStr', 'descend');
    
    %make sure that it is more than 4 standard deviations from the median
    pk_index = locs(find(stds_away * std(x(:,i)) < pks));
    NUM_SPIKE = length(pk_index);
    
    %check for a stronger delta in the one second following the shark peak.
    try
        for j = 1:length(locs)
            %is the power in the preceeding second less than in the seconds before the spike.
            delta_jump = bandpower(x(locs(j)-Fs:locs(j),i),Fs,[0.5 4])*dealta_jumps < bandpower(x(locs(j):locs(j)+Fs,i),Fs,[0.5 4]);
        end
        PEAK_AFTER_DELTA = length(delta_jump);
    catch
        PEAK_AFTER_DELTA = 0;
    end
    
    [pks, locs]=findpeaks(x(:,i)-mean(x(:,i)),'MinPeakWidth',14, 'MaxPeakWidth', 40, 'SortStr', 'descend');
    pk_index = locs(find(stds_away * std(x(:,i)) < pks));
    NUM_SHARP_WAVE = length(pk_index);
    out = NUM_SHARP_WAVE;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Burst Supression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THESE ARE THE TIME DOMAIN INFORMATION ON THE BS.
    %burst_length_mean
elseif featurenum == 40
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= mean(bur(:,2) - bur(:,1))/Fs;
    
    
    %burst_length_std
elseif featurenum == 41
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= std(bur(:,2) - bur(:,1))/Fs;
    
    
    %supression_length_mean
elseif featurenum == 42
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= mean(sup(:,2) - sup(:,1))/Fs;
    
    
    %supression_length_std
elseif featurenum == 43
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= std(sup(:,2) - sup(:,1))/Fs;
    
    
    %burst_supression_total
elseif featurenum == 44
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= sum(bur(:,2) - bur(:,1))/length(x);
    
    %num_bursts
elseif featurenum == 45
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= length(bur);
    
    %num_sup
elseif featurenum == 46
    [ bur, sup ] = Burst_supression( x(:,i), Fs )
    out= length(sup);
    
elseif featurenum == 47
    for i = 1:length(bur)
        this_burst = x(bur(i,1):bur(i,2));
        del(i) = bandpower(this_burst,Fs,[0.5,4]);
        the(i) = bandpower(this_burst,Fs,[4,7]);
        alp(i) = bandpower(this_burst,Fs,[8,15]);
        bet(i) = bandpower(this_burst,Fs,[16,31]);
        gam(i) = bandpower(this_burst,Fs,[32,49]);
        mu(i) = bandpower(this_burst,Fs,[8,12]);
    end
    
    %THIS IS SOME INFORMATION ON THE BANDPOWERS
    out=mean(del)
    
elseif featurenum == 48
    for i = 1:length(bur)
        this_burst = x(bur(i,1):bur(i,2));
        del(i) = bandpower(this_burst,Fs,[0.5,4]);
        the(i) = bandpower(this_burst,Fs,[4,7]);
        alp(i) = bandpower(this_burst,Fs,[8,15]);
        bet(i) = bandpower(this_burst,Fs,[16,31]);
        gam(i) = bandpower(this_burst,Fs,[32,49]);
        mu(i) = bandpower(this_burst,Fs,[8,12]);
    end
    out=mean(the)
    
elseif featurenum == 49
    for i = 1:length(bur)
        this_burst = x(bur(i,1):bur(i,2));
        del(i) = bandpower(this_burst,Fs,[0.5,4]);
        the(i) = bandpower(this_burst,Fs,[4,7]);
        alp(i) = bandpower(this_burst,Fs,[8,15]);
        bet(i) = bandpower(this_burst,Fs,[16,31]);
        gam(i) = bandpower(this_burst,Fs,[32,49]);
        mu(i) = bandpower(this_burst,Fs,[8,12]);
    end
    out=mean(alp)
    
elseif featurenum == 50
    for i = 1:length(bur)
        this_burst = x(bur(i,1):bur(i,2));
        del(i) = bandpower(this_burst,Fs,[0.5,4]);
        the(i) = bandpower(this_burst,Fs,[4,7]);
        alp(i) = bandpower(this_burst,Fs,[8,15]);
        bet(i) = bandpower(this_burst,Fs,[16,31]);
        gam(i) = bandpower(this_burst,Fs,[32,49]);
        mu(i) = bandpower(this_burst,Fs,[8,12]);
    end
    out=mean(bet)
    
elseif featurenum == 51
    for i = 1:length(bur)
        this_burst = x(bur(i,1):bur(i,2));
        del(i) = bandpower(this_burst,Fs,[0.5,4]);
        the(i) = bandpower(this_burst,Fs,[4,7]);
        alp(i) = bandpower(this_burst,Fs,[8,15]);
        bet(i) = bandpower(this_burst,Fs,[16,31]);
        gam(i) = bandpower(this_burst,Fs,[32,49]);
        mu(i) = bandpower(this_burst,Fs,[8,12]);
    end
    out=mean(gam)
    
elseif featurenum == 52
    for i = 1:length(bur)
        this_burst = x(bur(i,1):bur(i,2));
        del(i) = bandpower(this_burst,Fs,[0.5,4]);
        the(i) = bandpower(this_burst,Fs,[4,7]);
        alp(i) = bandpower(this_burst,Fs,[8,15]);
        bet(i) = bandpower(this_burst,Fs,[16,31]);
        gam(i) = bandpower(this_burst,Fs,[32,49]);
        mu(i) = bandpower(this_burst,Fs,[8,12]);
    end
    out=mean(mu)
    
end

end
