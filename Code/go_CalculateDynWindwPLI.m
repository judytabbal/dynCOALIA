function wPLIs = go_CalculateDynWindwPLI(dataEpochs,srate,totalWindow,step,freqBand,trial_length,ROIs_labels)

% Compute Connectivity matrix using the weighted phase lag index wPLI implemented in fieldtrip
% Computed for a subject, repetitions needed for wPLI are trials of the
% subject 

%wPLIs: output 3D matrix of dim [nb_win x nROIs x nROI] where nb_win=nb windows, nROI=number of ROIs in atlas

%dataEpochs: cell of length ntrials, dim: nROIs x time_samples
%srate: sampling rate 
%totalWindow: window width in sec (ex 0.6s)
%step: window step size in sec  (ex 0.05s if 91.67% overlapping)
%freqBand: freqency band used for (ex [13 30])
%trial_length: length of trials in sec (calculated as abs(prestim)+poststim)
%ROIs_labels: 1*nROIs cell containing the labels of atlas used


numberOfTrials=length(dataEpochs);
numberOfChannels=size(dataEpochs{1},1);
numberOfSamples=size(dataEpochs{1},2);

nb_win=length((0+totalWindow/2):step:(trial_length-totalWindow/2));

windowSamples=floor(totalWindow*srate);
stepSamples=floor(step*srate);

%vector times for extracting windowed signals
ti = 1+(0:(nb_win-1))*stepSamples;
tf = windowSamples+(0:(nb_win-1))*stepSamples;

wPLIs=zeros(nb_win,numberOfChannels,numberOfChannels);

%For every window, compute the wPLI averaged over repetitions (trials)
for i=1:nb_win    
    for ep=1:numberOfTrials
        ftData.time{ep} =  totalWindow*(i-1):1/srate:totalWindow*i-1/srate;
        ftData.trial{ep} = dataEpochs{ep}(:,ti(i):tf(i));
    end
    ftData.fsample = srate;   
    ftData.label = ROIs_labels;

    % cross-spectrum using FieldTrip
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'powandcsd';
    cfg.taper = 'hanning';
    cfg.foilim = freqBand; 
    cfg.tapsmofrq = 2;
    cfg.pad = 'nextpow2';
    cfg.keeptrials = 'yes';
    freq = ft_freqanalysis(cfg, ftData);

    % wpli using FieldTrip
    cfg = [];
    cfg.method = 'wpli';
    conn_wPLI = ft_connectivityanalysis(cfg,freq);
    conn_wPLI = mean(conn_wPLI.wplispctrm,2); %avg over all freq bins
    conn_wPLI = abs(conn_wPLI); 
    
    %convert channels combination format 1D (of the first dim of conn_wPLI) to chan*chan format 2D
    wPLI=zeros(numberOfChannels,numberOfChannels);
    matrix=ones(numberOfChannels,numberOfChannels);
    for c=1:numberOfChannels
        matrix(c,c)=0;
    end
    index = find(tril(matrix)~=0);
    tmp = zeros(numberOfChannels,numberOfChannels);
    tmp(index) = conn_wPLI;
    tmp = tmp + transpose(tmp);
    wPLI = tmp;
    
    wPLIs(i,:,:)=wPLI;
end

end