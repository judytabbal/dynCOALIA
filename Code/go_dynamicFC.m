function cmat=go_dynamicFC(opt,filters,data)

n_parcels = size(filters,1);
opt.trial_length=abs(opt.prestim)+abs(opt.poststim);

% determine the number of individual windows per trial
tmp = (0+opt.window.size/2):opt.window.step:(opt.trial_length-opt.window.size/2);
n_windows = length(tmp);
n_trials  = length(data.trial);
n_samples = floor(opt.window.size * data.fsample);
n_shifts  = floor(opt.window.step * data.fsample);

% work out which sample numbers are the start and end of each window
ti = 1+(0:(n_windows-1))*n_shifts;
tf = n_samples+(0:(n_windows-1))*n_shifts;

% start building the cmat structure
cmat                    = [];
% cmat.connectivity       = cell(1,n_trials);

disp('Generating connectivity matrices:')
ft_progress('init', 'text', 'Please wait...')

% loop through all trials to generate connectomes
for ii = 1:n_trials
    ft_progress(ii/n_trials, 'Processing trial %d from %d', ii, n_trials);
    % generate virtual electrodes
    VE = [];
    VE.raw = filters*data.trial{ii};
    data_input{ii}=VE.raw;
end

eegData=zeros(n_parcels,size(data_input{1},2),n_trials);
for tr=1:n_trials
    eegData(:,:,tr)=data_input{tr};
end

fmin=opt.bpfreq(1);
fmax=opt.bpfreq(2);

tmp2 = data.time{1};
for ii = 1:n_windows
    time_wind(ii) = mean(tmp2(ti(ii):tf(ii)));
end

ROIs_labels=opt.labels;

switch opt.conn_method
     case 'plv_inst_pn'
         filtSpec.range=[fmin fmax]; filtSpec.order=100;
         [plvs_inst] = pn_eegPLV(eegData, data.fsample, filtSpec);
         plvs_inst_perm=permute(plvs_inst,[2 3 1]);
         for w=1:opt.trial_length*data.fsample
             plvs_inst_sym(:,:,w)=squeeze(plvs_inst_perm(:,:,w))+squeeze(plvs_inst_perm(:,:,w))';
         end
         cmat.connectivity{1}=plvs_inst_sym;
         cmat.time=tmp2;
         cmat.n_trials=1;
    case 'plv_dyn'
         [plvs,N]=go_CalculateDynWindPLV(data_input,data.fsample,fmin,fmax,opt.window.size,opt.window.size-opt.window.step,opt.window.step,n_windows);
         for t=1:n_trials
             plvs_perm{t}=permute(plvs{t},[2 3 1]);
             for w=1:n_windows
                 plvs_sym{t}(:,:,w)=squeeze(plvs_perm{t}(:,:,w))+squeeze(plvs_perm{t}(:,:,w))';
             end
             cmat.connectivity{t}=plvs_sym{t};
         end
         cmat.time=time_wind;
         cmat.n_trials = n_trials;
    case 'wPLI'
         wPLIs_ft=go_CalculateDynWindwPLI(data_input,data.fsample,opt.window.size,opt.window.step,[fmin fmax],opt.trial_length,ROIs_labels);
         wPLIs_ft_perm=permute(wPLIs_ft,[2 3 1]);
         cmat.connectivity{1}=wPLIs_ft_perm;
         cmat.time=time_wind;
         cmat.n_trials=1;
    case 'pdd'
         for trr=1:n_trials
             [pdd, pddnode, meanpdd] = dyn_pdd(data_input{trr}(:,900:end),data.fsample,floor(1/floor((fmin+fmax)/2)*data.fsample),floor((fmin+fmax)/2));
             cmat.connectivity{trr}=pdd;
         end
 tmp2=[-126/data.fsample:1:opt.poststim]; %-126:1:post for pdd and iac
         cmat.time=tmp2(26:end-26);
         cmat.n_trials = n_trials;
    case 'iac'
         for trr=1:n_trials
             [dyn_IAC, mean_IAC] = dyn_env(data_input{trr}(:,900:end));
             cmat.connectivity{trr}=dyn_IAC;
         end
  tmp2=[-126/data.fsample:1:opt.poststim]; %-126:1:post for pdd and iac
        cmat.time=tmp2(26:end-26);
        cmat.n_trials = n_trials;
    case 'aec'
         for trr=1:n_trials
             for jj = 1:n_windows
                 VE.windowed = transpose(data_input{trr}(:,ti(jj):tf(jj))); 
                 VE.H        = abs(hilbert(VE.windowed));
                 cmat_trial(:,:,jj) = corrcoef(VE.H) - eye(n_parcels);
             end
             cmat.connectivity{trr}=cmat_trial;
         end
         cmat.time=time_wind;
         cmat.n_trials = n_trials;
    end

ft_progress('close')
disp('DONE')

% add final information
cmat.window             = opt.window;
cmat.bpfreq             = opt.bpfreq;
cmat.n_parcels          = n_parcels;
cmat.conn_type          = opt.conn_method;