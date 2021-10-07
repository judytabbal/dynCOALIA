function [ kernel ] = ComputeWMNE( noiseCov,gain,gridLoc,gridOrient,weightExp,weightLimit,SNR )
% Compute the kernel matrix by solving the inverse problem
% inputs: - noiseCov: the noise covariance matrix, nb_channels*nb_channels
%         - gain: the gain Matrix (leadfield matrix ) obtained after
%         computing the head model, nb_channels*(3*nb_regions), nb_regions denotes the number
%         of cortical sources, nb_channels denotes the number of EEG channels.
%         - gridLoc: the position (x,y,x) of the sources, nb_regions*3
%         - gridOrient: the orientation of the sources, nb_regions*3
%         - weighExp: parameter from brainstorm to compute the wMNE
%         - weightLimit: parameter from brainstorm to compute the wMNE
%         - SNR: signal to noise ratio
%
% Output: - kernel: nb_regions*nb_channels, nb_regions denotes the number
%           of cortical sources, nb_channels denotes the number of EEG channels.

% This code was originally developped by Sahar Yassine and Ahmad Mheich based on the brainstorm toolbox codes 
% (Tadel, F., Baillet, S., Mosher, J.C., Pantazis, D., Leahy, R.M., 2011.
% Brainstorm: A user-friendly application for MEG/EEG analysis. Computational Intelligence and Neuroscience 2011. https://doi.org/10.1155/2011/879716)

% contact: saharyassine94@gmail.com
%          mheich.ahmad@gmail.com

%% Compute the inverse problem

if(size(noiseCov,2)~= size(gain,1))
    if(size(noiseCov,2)<size(gain,1))
        gain=gain(1:size(noiseCov,2),:);
    else
        error('unmatched dimensions between the gain matrix and the noise covariance');
    end
end
numberOfChannels=size(noiseCov,1);
montageMatrix=eye(numberOfChannels)-(ones(numberOfChannels,numberOfChannels)./numberOfChannels);
CNoise=(montageMatrix*noiseCov)*montageMatrix';
CNoise=(CNoise+CNoise')./2;
gain=montageMatrix*gain;
% CNoise=diag(diag(CNoise));

[alpha,gainWQ]=sourceModelAssumption(gain,gridLoc,gridOrient,weightExp,weightLimit);
[~,iW]=TruncateAndRegularizeCovariance(CNoise);
L=iW*gainWQ;
[UL,SL2] = svd((L*L'));
SL2 = diag(SL2);
SL = sqrt(SL2); 
SNR=SNR*SNR;
lambda=SNR/mean(SL.^2);

kernel = lambda * L' * (UL * diag(1./(lambda * SL2 + 1)) * UL');
kernel = kernel * iW;
for line=1:size(kernel,1)
    kernel(line,:)=alpha(line)*kernel(line,:);
end


end

