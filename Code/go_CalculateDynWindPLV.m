function [ plvs,N ] = go_CalculateDynWindPLV(dataEpochs,samplingRate,fmin,fmax,totalWindow,overLapWindow,step,n_wind )

%plvs: computed for each subject, output format cell length Nt, each 3D [N x nROI x nROI] where
%N=nb windows, Nt=nb of trials for this subject
%N: number of windows used
%dataEpochs: cell of length ntrials, dim: nROI x samples
%totalWindow: window length in sec (ex computed for 6 gamma cycles)
%overLapWindow: length of overlapping in sec(ex 90%*totalWindow)
%step: total-overlap (in sec)

numberOfTrials=length(dataEpochs);
for ep=1:numberOfTrials
    numberOfChannels=size(dataEpochs{ep},1);
    numberOfSamples(ep)=size(dataEpochs{ep},2);
end

filterOrder=floor(0.1*samplingRate);

% filtering with FIR
b1 = fir1(filterOrder,[fmin fmax]/(samplingRate/2));
filteredData=cell(numberOfTrials,1);

for k = 1:numberOfChannels
    for trial=1:numberOfTrials
        filteredData{trial}(k,:) = filtfilt(b1,1,double(dataEpochs{trial}(k,:)));
    end
end

for nb_tr=1:numberOfTrials
    for channelCount = 1:numberOfChannels
        filteredData{nb_tr}(channelCount,:) = angle(hilbert(squeeze(filteredData{nb_tr}(channelCount,:))));
    end
end

fInterest=fmin+(fmax-fmin)/2;
windowSamples=floor(totalWindow*samplingRate);
overLapSamples= floor(overLapWindow*samplingRate);
stepSamples=floor(step*samplingRate);

tmp = (0+totalWindow/2):step:(numberOfSamples/1000-totalWindow/2);
N=n_wind;

plvs=[];

for trial=1:numberOfTrials
    plv=zeros(N,numberOfChannels,numberOfChannels);
    No=floor(windowSamples/2)+1;
    for count=1:N
        for channelCount = 1:numberOfChannels-1
            channelData = squeeze(filteredData{trial}(channelCount, No-floor(windowSamples/2):No+floor(windowSamples/2)-1));
            for compareChannelCount = channelCount+1:numberOfChannels
                compareChannelData = squeeze(filteredData{trial}(compareChannelCount, No-floor(windowSamples/2):No+floor(windowSamples/2)-1));
                diff=channelData(:, :) - compareChannelData(:, :);
                diff=diff';
                plv(count,channelCount,compareChannelCount) =abs(sum(exp(1i*diff)))/length(diff);  
            end
        end
        No=No+stepSamples;
    end
    plv = squeeze(plv);
    plvs{trial}=plv;
end
end