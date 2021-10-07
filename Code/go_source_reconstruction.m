function filters=go_source_reconstruction(opt,data,subvol,subgrid,sources_orients)

%Calculate Lead fields (Forward Model)
cfg             = [];
cfg.headmodel   = subvol;
cfg.elec        = data.elec;
cfg.grid.pos    = subgrid.pos;
cfg.normalise   = 'yes';
cfg.rankreduce  = 3; 
lf              = ft_prepare_leadfield(cfg); 

% load('E:\DynamicEEG\from BS\DIRECT_Methods\lf_15002to66LFP.mat');
% lf=lf_15002to66LFP;

if(opt.meth==1) %wMNE
    
    %Timelock for noise covariance estimation
    cfg                     = [];
    cfg.covariance          = 'yes';
    cfg.window              = [-opt.prestim -0.0]; %baseline for noise cov (computed from prestim in seconds to 0)
    tlk_noise                = ft_timelockanalysis(cfg,data);
    noise_cov_ft=tlk_noise.cov;

    for nn=1:length(data.trial)
        data_tr_all(:,:,nn)=data.trial{nn};
    end
    eeg_for_noise=mean(data_tr_all,3);
%     noise_cov=CalculateNoiseCovarianceTimeWindow(eeg_for_noise(:,1:opt.prestim*data.fsample));

    filters = ComputeWMNE(noise_cov_ft,cell2mat(lf.leadfield),lf.pos(lf.inside,:),sources_orients(lf.inside,:),opt.weightExp,opt.weightLimit,opt.SNR);

elseif(opt.meth==2) %eLORETA
    cfg                     = [];
    cfg.method              = 'eloreta';
    cfg.grid                = lf;
    cfg.grid.mom            = transpose(sources_orients);
    cfg.headmodel           = subvol;
    cfg.eloreta.keepfilter  = 'yes';
    cfg.eloreta.keepmom     = 'no';
    cfg.eloreta.lambda      = 0.05;
    src                     = ft_sourceanalysis(cfg,data);
    filters=cell2mat(transpose(src.avg.filter));

elseif(opt.meth==3) %LCMV Beamforming
    
    % Generate Covariance for Beamforming
    cfg                     = [];
    cfg.covariance          = 'yes';
    cfg.window              = [-opt.prestim -0.0]; %baseline for noise cov (computed from prestim in seconds to 0)
    cfg.keeptrials          = 'yes';
    cfg.vartrllength        = 2;
    timelock                = ft_timelockanalysis(cfg,data);

    tmp = squeeze(mean(timelock.cov));
    lambda = 0.05*max(svd(tmp)) - min(svd(tmp));

    % Generate Beamformer weights
    cfg                 = [];
    cfg.method          = 'lcmv';
    cfg.grid            = lf;
    cfg.vol             = subvol;
    cfg.keepfilter      = 'yes';
    cfg.lcmv.fixedori   = 'yes';
    cfg.lcmv.lambda     = lambda;  % Heavy regularisation to blurr data a bit
    src                 = ft_sourceanalysis(cfg,timelock);

    % extract the beaforming filters (weights) for next part of the script;
    filters = cell2mat(src.avg.filter);
end
    
end

% % uncomment below only for source plot visualization (indicate interval of time to average in cfg.latency)
% cfg                 = [];
% cfg.method          = 'mne';
% cfg.grid            = lf;
% cfg.elec            = data.elec;
% cfg.headmodel       = subvol;
% cfg.mne.snr         = 3; %regularisation param for noise cov matrix 
% src                 = ft_sourceanalysis(cfg,data);
% 
% nROIs=66;
% nSamples=1168;
% nTrials=200;
% scout_source=zeros(nROIs,nSamples);
% for tr=1:nTrials
%     scout_source=scout_source+abs(filters*data.trial{tr}); 
% end
% scout_source=scout_source./nTrials;
% src.aa=scout_source;
% 
% cfg = [];
% cfg.funparameter = 'aa';
% cfg.maskparameter = 'aa';
% cfg.method = 'surface';
% cfg.surffile= 'D:\Brainstorm\brainstorm_db\Protocol05\anat\Subject01\tess_cortex_pial_low.mat';
% cfg.latency = [0.190 0.340];
% cfg.avgovertime = 'yes';
% max_avg=max(mean(src.aa(:,floor((cfg.latency(1)+0.3)*data.fsample):floor((cfg.latency(2)+0.3)*data.fsample)),2));
% cfg.opacitylim = [0.2*max_avg max_avg];
% ft_sourceplot(cfg, src);