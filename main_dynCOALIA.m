% 1- Define your path base and some params of your data

path_base='E:\DynamicEEG\dynCOALIA';
n_sub=20; %nb of subjects
n_trials=100; %nb of trials per subject
n_parcels=66; %nb of parcels (ROIs)
n_chans=257; %nb of channels (EEG)
bpfreq=[30 40]; %frequency band of interest
fs=1024; %sampling frequency
pre_samples=1024; %nb samples before trial onset 
post_samples=1024; %nb samples after trial onset
nb_samples=2049; %total number of samples in each trial (again, if variable between subjects, define this param inside the for loop subjects, preferable to be consistent between trials for a single subject for later concatenation)
onset=1024+1;


%% 2-Calculate Headmodel (BEM)using OpenMEEG

load([path_base 'Inputs\mri_colin_realign.mat']); 

cfg           = [];
cfg.output    = {'brain','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, mri_colin_realign); %ctf/mm

cfg=[];
cfg.method='headshape';
brain=ft_read_headshape([path_base 'Inputs\tess_innerskull_colin.mat']);
brain_mm=ft_convert_units(brain,'mm');
cfg.headshape=brain_mm;
cfg.numvertices = [3000];
bnd(1)=ft_prepare_mesh(cfg,segmentedmri);

cfg=[];
cfg.method='headshape'; 
skull=ft_read_headshape([path_base 'Inputs\tess_outerskull_colin.mat']);
skull_mm=ft_convert_units(skull,'mm');
cfg.headshape=skull_mm;
cfg.numvertices = [3000];
bnd(2)=ft_prepare_mesh(cfg,segmentedmri);

cfg=[];
cfg.method='headshape';
head=ft_read_headshape([path_base 'Inputs\tess_head_colin.mat']);
head_mm=ft_convert_units(head,'mm');
cfg.headshape=head_mm;
cfg.numvertices = [3000];
bnd(3)=ft_prepare_mesh(cfg,segmentedmri);

figure();
ft_plot_mesh(bnd(1), 'edgecolor', 'none', 'facecolor', 'r')
ft_plot_mesh(bnd(2), 'edgecolor', 'none', 'facecolor', 'g')
ft_plot_mesh(bnd(3), 'edgecolor', 'none', 'facecolor', 'b')
alpha 0.3

cfg        = [];
cfg.method ='openmeeg'; % You can also specify 'openmeeg', 'bemcp', or another method.
subvol        = ft_prepare_headmodel(cfg, bnd); %ctf/mm


%% 3- Calculate Lead fields (Forward Model)

%Add openmeeg path and default fieldtrip path
setenv('PATH', [path_base '\Toolboxes\OpenMEEG\bin'])
addpath([path_base '\Toolboxes\fieldtrip-20190224\external\openmeeg'])
tmp = which('ft_defaults');
if isempty(tmp)
    addpath([path_base '\Toolboxes\fieldtrip-20190224']); %add defaults fieldtrip
    ft_defaults
end

%Load necessary input variables
load([path_base 'Inputs\scout_scs_66LFP.mat']); 
scout_scs=scout_scs_66LFP; source_Orient = transpose(scout_scs.orients); %scout position and orientation of the 66 desikan ROIs
load([path_base 'Inputs\subgrid_66LFP.mat']);  %subject grid source pos
load([path_base 'Inputs\lf_66LFP.mat']); ftLeadfield=lf_66LFP; %leadfield computed on 66 desikan ROIs

leadfields_const = zeros(257,66);

% constrain the orientation of the sources to the normal to the surface
for i=1:66
    leadfields_const(:,i) = ftLeadfield.leadfield{i}*source_Orient(:,i);
end

% load data (make your own data following this structure and load it
% instead)
% data_all_LFP: LFP simulated data for all subj and trials
load([path_base 'Inputs\data_all_LFP.mat']); 

for sub=1:n_sub
    for j=1:n_trials   
        %filter LFP data in the band frequency of interest
        filterorder=0.1;
        b1 = fir1(floor(filterorder*fs),[bpfreq(1) bpfreq(2)]/(fs/2));
        for kk = 1:n_parcels
            B_demean_deref_filt(kk,:) = filtfilt(b1,1,double(data_all_LFP{sub}{j}(kk,:))); %filter LFP
        end
        %generate EEG simulated data using the calculated forward model
        eeg_comp_lf66_filt{sub}(:,:,j)=leadfields_const*B_demean_deref_filt; 
    end
end


%% 4- Define Sources and dynamic Functional Connectivity (dFC) structures

%Automatically add folders and subfolders (Code and Inputs) necessary for
%loading variables and executing codes
folder1=[path_base '\Inputs'];
folder2=[path_base '\Code'];
addpath(genpath(folder1));
addpath(genpath(folder2));

% Construct Source Structure
source.weightExp   = 0.5; %param for depth weighting
source.weightLimit = 10; %param for depth weighting limit
source.SNR         = 3; %Signal to Noise Ratio
source.meth        = 1; %1:wMNE,2:eLORETA

% Construct dynamic Functional Connectivity (dFC) Structue
dFC.conn_method = 'plv_dyn'; %'plv_dyn' or 'wPLI' or 'aec'
dFC.bpfreq      = bpfreq; %frequency band of interest
dFC.window.size = 0.17; % sliding window length in seconds (for example calculated for 6cycles,CentralFreq=35 ==> 6/35=0.17s)
dFC.window.step = 0.017; % step between windows in seconds (for example 90% overlapping=10/100*window_size)

load([path_base '\Inputs\elec_BS_colin.mat']); %fieldtrip format electrode already computed from EGI257 system: extracted from BrainStorm
load([path_base '\Inputs\desikan_labels_66_LFP.mat']); %cell labels for the scout regions: extracted from BrainStorm


%% 5- Compute Sources FILTERS 

label={};
for i_chan=1:nb_chan-1
    label{i_chan}=['E' int2str(i_chan)];
end
label{nb_chan}='Cz';   

data=[];
filters=[];
   
for j=1:nb_trials
    data.label=label;
    data.fsample=fs;   
    data.trial{j}= eeg_comp_lf66_filt{1}(:,:,j);
    data.time{j}=[-pre_samples:1:post_samples]/fs; 
end

data.elec=elec_BS_mm;
    
%Source Reconstruction
cfg             = [];
cfg.prestim     = pre_samples/fs; %prestimulus in seconds for noise cov computation
cfg.weightExp   = source.weightExp; %param for depth weighting
cfg.weightLimit = source.weightLimit; %param for depth weighting limit
cfg.SNR         = source.SNR; %Signal to Noise Ratio
cfg.meth        = source.meth; %wMNE or eLORETA method
filters   = go_source_reconstruction(cfg,data,subvol,subgrid_66LFP,scout_scs.orients);


%% Compute connectivities for reconstructed sources (after obtaining filters)

%Loop over all subjects
for sub_ind = 1:nb_sub
    data=[];
    cmat=[];
    
    for j=1:nb_trials
        data.label=label;
        data.fsample=fs;
        data.trial{j}= eeg_comp_lf66_filt{sub_ind}(:,:,j);
        data.time{j}=[-pre_samples:1:post_samples]/fs; 
    end
    
    data.elec=elec_BS_mm;   
    data_all_eeg{sub_ind}=data.trial;
    
    for tr=1:nb_trials
        data_all_rec_sources{sub_ind}{tr}=filters*data_all_eeg{sub_ind}{tr};
    end
    
    %Dynamic FC Computation
    cfg             = [];
    cfg.window.size = dFC.window.size; %sliding window length in seconds 
    cfg.window.step = dFC.window.step; %step between windows in seconds 
    cfg.bpfreq      = dFC.bpfreq; %frequency band of interest
    cfg.prestim     = pre_samples/fs; %prestimulus time in seconds
    cfg.poststim    = post_samples/fs; %poststimulus time in seconds
    cfg.conn_method = dFC.conn_method; %3 choices: 'plv_inst_pn' for instantaneousPLV or 'plv_dyn' for windowedPLV or 'wPLI' for windowedwPLI
    cfg.labels      = scout_labels;
    cmat            = go_dynamicFC(cfg,filters,data);
    
    cmat_allsub{sub_ind}=cmat; %structure containing connectivity matrices of all

end


%% Do the decomposition to reconstruct dynamic Brain Network States (JADE/FastICA/PSAUD/PCA/NMF/Kmeans)

% to specify by the user
val=6; %1:JADE,2:FastICA,3:PSAUD,4:PCA,5:NMF,6:Kmeans
num_comp=6; %number of states

cmat_all=[];

for ii = 1:n_sub    
    cmat=cmat_allsub{ii};
    
    cmat_sub = [];
    cmat_2d  = [];
    cmat_red = [];
    
    % loop across trials
    for jj = 1:cmat.n_trials
        tmp = cmat.connectivity{jj};
        cmat_sub = cat(3,cmat_sub,tmp);
    end   
    
    % Do a DC correction on each connection
    cmat_sub = cmat_sub-repmat(mean(cmat_sub,3),[1 1 size(cmat_sub,3)]);   
    % flatten tensor into 2.5d (as the machine learning kids call it)
    if(size(cmat_sub,3)==length(cmat.time)-1)
        inst=1;
        cmat_2d = reshape(cmat_sub,n_parcels*n_parcels,(length(cmat.time)-1)*cmat.n_trials);
    else
        inst=0;
        cmat_2d = reshape(cmat_sub,n_parcels*n_parcels,length(cmat.time)*cmat.n_trials);
    end
    % as each slice of the tensor is symmetric, remove redundant data to save some memory
    index = find(tril(squeeze(cmat_sub(:,:,1)))~=0);
    cmat_red = cmat_2d(index,:);     
    % concatenate onto the group level connectome data
    cmat_all = cat(2,cmat_all,cmat_red);
    
    sub_trials(ii) = cmat.n_trials;   
end

disp('DONE')

if(val==1) %ICA-JADE
    [B] = jader(cmat_all,NCs);
    A=B'; 
    sig=B*cmat_all; 
elseif(val==2) %ICA-InfoMax
    % add the fastica toolbox to path if its not already present
    tmp = which('fastica');
    if isempty(tmp)
        [status] = ft_hastoolbox('fastica', 1, 0);
        if ~status
            error('fastica toolbox failed to load correctly')
        end
    end
    [iq, A, W, sig, sR]=icasso(cmat_all,10,'g','tanh','lastEig',NCs,'numOfIC',NCs,'approach','defl','vis','off');
elseif(val==3) %ICA-PSAUD
    tau=1;
    P=NCs;
    Nit=20;
    alpha_min=0;
    alpha_max=50;
    [F,S]=P_SAUD(cmat_all,tau,P,Nit,alpha_min,alpha_max);
    A=F; %estimated mixing matrix
    sig=S; %estimated signal matrix
elseif(val==4) %PCA
    [coeff, A, latent, tsquared, explained] = pca(cmat_all,'NumComponents', NCs);
    sig=coeff';
elseif(val==5) %NMF
    [A,sig] = nnmf(cmat_all,NCs,'rep',10);
elseif(val==6) %Kmeans
    [idxx,C]=kmeans(cmat_all',NCs,'Distance','cityblock','MaxIter',200,'Start','Sample','Replicates',1);
    A=C';
end

% post processing of the results 
if(val==6) %kmeans
    signal=[];
    idx_2d = reshape(idxx,size(idxx,1)/sum(sub_trials),sum(sub_trials));
    for nn=1:NCs
        iddd=find(idx_2d==nn);
        idx_2d_new=idx_2d; idx_2d_new(:)=0; idx_2d_new(iddd)=1;
        idx_2d_new_all(nn,:,:)=idx_2d_new;
    end
    for ii=1:NCs
        idx_2d_new_all_normal(ii,:,:)=zscore(idx_2d_new_all(ii,:,:));
    end
    results.signals = idx_2d_new_all_normal;
else
    tmp = reshape(sig,NCs,size(sig,2)/sum(sub_trials),sum(sub_trials));
    results.signals = tmp;
end

for ii = 1:NCs
    tmp = zeros(n_parcels,n_parcels);
    tmp(index) = A(:,ii);
    tmp = tmp + transpose(tmp);
    results.maps(:,:,ii) = tmp;
end

results.NCs         = NCs;
results.n_trials    = sum(sub_trials);
results.sub_trials  = sub_trials;
if(inst==1)
    results.time    = round(cmat.time(1:end-1)*100)/100;
else
    results.time    = round(cmat.time*100)/100;
end
if(val~=6)
    results.sig=sig;
end
results.mixing=A;
results.cmat_all=cmat_all;

results.maps=abs(results.maps);

