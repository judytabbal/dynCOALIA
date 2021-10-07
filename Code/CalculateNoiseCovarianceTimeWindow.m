function [ noiseCov ] = CalculateNoiseCovarianceTimeWindow( timeWindow )

% Compute noise covariance matrix
% inputs: timewindow, nb_channels*nb_samples
% outputs: noiseCov, noise covariance matrix, nb_channels*nb_channels

% This code was originally developped by Sahar Yassine and Ahmad Mheich based on the brainstorm toolbox codes 
% (Tadel, F., Baillet, S., Mosher, J.C., Pantazis, D., Leahy, R.M., 2011.
% Brainstorm: A user-friendly application for MEG/EEG analysis. Computational Intelligence and Neuroscience 2011. https://doi.org/10.1155/2011/879716)

% contact: saharyassine94@gmail.com
%          mheich.ahmad@gmail.com

%% Noise Covariance
[~,nbSamples]=size(timeWindow);
%remove DC Offset
baseTimeWindow=RemoveDCOffset(timeWindow);
noiseCov=baseTimeWindow*baseTimeWindow';
noiseCov=noiseCov./nbSamples;

end

