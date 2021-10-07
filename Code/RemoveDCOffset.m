function [ noOffsetEpoch ] = RemoveDCOffset( originalEpoch )

% removing DC offset
% inputs: originalEpoch: data, nb_channels*nb_samples
% output: noOffsetEpoch, nb_channels*nb_samples

% This code was originally developped by Sahar Yassine and Ahmad Mheich based on the brainstorm toolbox codes 
% (Tadel, F., Baillet, S., Mosher, J.C., Pantazis, D., Leahy, R.M., 2011.
% Brainstorm: A user-friendly application for MEG/EEG analysis. Computational Intelligence and Neuroscience 2011. https://doi.org/10.1155/2011/879716)

% contact: saharyassine94@gmail.com
%          mheich.ahmad@gmail.com

originalEpoch=double(originalEpoch);
[nbChannels,nbSamples]=size(originalEpoch);
noOffsetEpoch=zeros(nbChannels,nbSamples);

for i=1:nbChannels
    noOffsetEpoch(i,:)=originalEpoch(i,:)-mean(originalEpoch(i,:));
end
end

