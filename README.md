# Features For EEG Analysis
This repository contains code to generate features that quantify the complexity, category and connectivity of EEG waveforms.

# Code
All code was written in MATLAB 2016a

"ALL_FEATURES_ONE_CH_V2.m" contains the code used to generate features of complexity and category.
"ALL_FEATURES_TWO_CH_V2.m" contains the code used to generate features of connectivity.

# Features
Our features of EEG complexity include: Shannon Entropy, Tsalis Entropy (with q ranging from 1:10), Cepstrum coefficients, Sub-band Information quantity, Lyaponov exponent, fractal dimension, Hjorth mobility/complexity, false nearest neighbor embedding dimension and the coefficients of a second order auto-regressive moving model.

Our selected features of connectivity include: coherence in delta band, coherence in all bands, phase lag index, cross correlation magnitude, cross correlation lag, mutual information and granger causality.

Our selected features of category include: standard deviation, signal regularity, EEG frequency band power (delta, theta, alpha, beta, gamma, mu), the ratio of alpha-to-delta band power, signal amplitude less than 5 microV,  signal amplitude less than 10 microV, signal amplitude less than 20 microV, "normal" EEG, diffuse slowing, number of epileptiform spikes, epileptiform peaks followed by increase in delta band power. We quantified burst suppression using the following features: burst length (mean and std), suppression length (mean and std), number of bursts, number of suppressions.

Please see the paper for a full description of the features.
