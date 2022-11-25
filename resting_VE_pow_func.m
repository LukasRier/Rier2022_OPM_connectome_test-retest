function resting_VE_pow_func(sub,ses,project_dir)
restoredefaultpath
cleaning_only = 0;
close all
clc
% 
% sub = '001';
% ses = '001';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([project_dir,'scripts',filesep,'fieldtrip-20190212'])
addpath([project_dir,'scripts'])
addpath([project_dir,'scripts',filesep,'Beamformer',filesep,''])
ft_defaults;
datadir = [project_dir,'data',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_type = 'task-movie';
filename = ['sub-',sub,'_',exp_type,'_','run-',ses];
path.main = [datadir,'sub-',sub,filesep];
path.ICA = [datadir,'derivatives',filesep,'ICA',filesep,'sub-',sub,filesep];
path.cleaning = [datadir,'derivatives',filesep,'cleaning',filesep,'sub-',sub,filesep];
path.meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
path.VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
path.pow = [datadir,'derivatives',filesep,'POW',filesep,'sub-',sub,filesep];
path.noise = [path.main,'noise',filesep];

if ~exist(path.ICA,'dir'); mkdir(path.ICA);end
if ~exist(path.cleaning,'dir'); mkdir(path.cleaning);end
if ~exist(path.VEs,'dir'); mkdir(path.VEs);end
if ~exist(path.pow,'dir'); mkdir(path.pow);end

path.data = [path.main,'meg',filesep];
path.mri = [path.main,'anat',filesep];
files.mri = ['sub-',sub,'.nii'];
files.meshes = ['sub-',sub,'_meshes.mat'];
files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
files.pow = ['sub-',sub,'_','run-',ses,'_POW'];
files.voxlox = ['sub-',sub,'_voxlox.mat'];
files.channels = '_channels.tsv';
files.sens_order = '_sensor_order.mat';
files.ICA = [path.ICA,filename];
files.cleaning = [path.cleaning,filename];
files.noise_channels = ['_channels.tsv'];
files.noise = ['sub-',sub,'_task-noise'];

S.mri_file = [path.meshes,'sub-',sub,'_',filesep,'x.nii']; % only need path for beamformer and single shell function. bit hacky, FIX
% global files
cd(path.data)
load([path.data,filename,'_meg.mat'],'fs','data','triggers')
load([datadir,'derivatives',filesep,'helmet',filesep,'M1p5_adult.mat'])
% basefilename = [path.main,'ses-',ses,'/Matt/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load MRI and meshes
mri = ft_read_mri([path.mri,files.mri])
load([path.meshes,files.meshes],'segmentedmri','meshes');
load([path.meshes,files.voxlox],'voxlox');


%% Preproc
% Mean correct
data = data - mean(data,1);

% Notch filter
Wo = 50/(fs/2);  BW = Wo/35;
[b,a] = iirnotch(Wo,BW);

disp('Applying Notch filter')
data = filter(b,a,data,[],1);
Nchans = size(data,2);

% bandpass filter for viewing
disp('Applying 1-150 Hz bandpass filter')

hp = 1;
lp = 150;
[b,a] = butter(4,2*[hp lp]/fs);
data_f = [filtfilt(b,a,data)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of bad channels
%% Get rid of bad channels

if ~exist([path.data,filename,files.channels,'_new'],'file')
    ch_table = readtable([path.data,filename,files.channels],'FileType','text','Delimiter','tab');
    
    ch_table.isx = endsWith(ch_table.name,'[X]');
    ch_table.isy = endsWith(ch_table.name,'[Y]');
    ch_table.isz = endsWith(ch_table.name,'[Z]');
    ch_table.slot_no = zeros(height(ch_table),1);
    % sanity check
    if sum(sum(ch_table{:,11:13},2)) ~= height(ch_table)
        error('Channel orientation [x,y,z] labels might be wrong!')
    end
    
    % bad channels from tsv
    % view all data w tsv bad ones marked
    [ch_table] = Bad_Channels(data_f',ch_table,fs,hp,lp);
    writetable(ch_table,[path.data,filename,files.channels,'_new'],...
        'WriteRowNames',true,'Delimiter','tab','FileType','text')
else
    ch_table = readtable([path.data,filename,files.channels,'_new'],...
        'Delimiter','tab','FileType','text');
end
%% Bad channels in noise recording
if ~exist(sprintf('%s%s%s_new',path.noise,files.noise,files.noise_channels),'file')
    noise_ch_table = readtable(sprintf('%s%s%s',path.noise,files.noise,files.noise_channels),'FileType','text','Delimiter','tab');
    
    noise_ch_table.isx = endsWith(noise_ch_table.name,'[X]');
    noise_ch_table.isy = endsWith(noise_ch_table.name,'[Y]');
    noise_ch_table.isz = endsWith(noise_ch_table.name,'[Z]');
    noise_ch_table.slot_no = zeros(height(noise_ch_table),1);
    % sanity check
    if sum(sum(noise_ch_table{:,11:13},2)) ~= height(noise_ch_table)
        error('[ Noise recording ] : Channel orientation [x,y,z] labels might be wrong!')
    end
    
    % bad channels from tsv
    % view all data w tsv bad ones marked
    [noise_ch_table] = Bad_Channels(noise_data_f',noise_ch_table,fs,hp,lp);
    writetable(noise_ch_table,sprintf('%s%s%s_new',path.noise,files.noise,files.noise_channels),...
        'WriteRowNames',true,'Delimiter','tab','FileType','text')
else
    noise_ch_table = readtable(sprintf('%s%s%s_new',path.noise,files.noise,files.noise_channels),...
        'Delimiter','tab','FileType','text');
end
%%

load([path.data,filename,files.sens_order],'T')
for sl_i = 1:height(T)
    if ~isempty(T{sl_i,1}{1})
        ch_table.slot_no(startsWith(ch_table.name,T{sl_i,1}{1}(1:2))) = sl_i;
    end
end
%% make sure noise and data include same channels
for ch_i = 1:height(ch_table)
    ind_in_noise_data = find(startsWith(noise_ch_table.name,ch_table.name(ch_i)));
    if isempty(ind_in_noise_data)
       ch_table.status(ch_i) = {'bad'};
    elseif strcmp(noise_ch_table.status{ind_in_noise_data},'bad')
       ch_table.status(ch_i) = {'bad'};
    end
end
for ch_i = 1:height(noise_ch_table)
    ind_in_data = find(startsWith(ch_table.name,noise_ch_table.name(ch_i)));
    if isempty(ind_in_data)
       noise_ch_table.status(ch_i) = {'bad'};
    elseif strcmp(ch_table.status{ind_in_data},'bad')
       noise_ch_table.status(ch_i) = {'bad'};
    end
end
clear noise_ch_table

% remove from data matrix
disp("Removing bad channels")
bad_chans = find(startsWith(ch_table.status,'bad'))
ch_table(bad_chans,:) = [];
data_f(bad_chans,:) = [];
% data(:,bad_chans) = [];
% make sensor orientation struct here from table with bad chans removed
% Plot all channels together
%% sensor info
S.sensor_info.pos = [ch_table.Px,ch_table.Py,ch_table.Pz];
S.sensor_info.ors = [ch_table.Ox,ch_table.Oy,ch_table.Oz];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%% Mean field correction  stuff (later?)
% N = sens_info.ors_used;
N = S.sensor_info.ors; % orientation matrix (N_sens x 3)
S.M = eye(length(N)) - N*pinv(N);

%% Epoch data
% segment using trigger
disp("Chopping data into 5s epochs")
tr = diff(triggers(:,1)>2);
start_sample1 = find(tr == 1);
start_sample = start_sample1(1);
end_sample = find(tr == -1);
duration = floor((end_sample - start_sample)/fs);

% data = data(start_sample+1:start_sample+duration*fs,:);
data_f = data_f(:,start_sample+1:start_sample+duration*fs);

% View data again
%   Bad_Channels(data_f',ch_table,fs,hp,lp);

% chop into short segments
epoch_length = 5;
data_f_mat = reshape(data_f,[size(data_f,1),epoch_length*fs,duration/epoch_length]);

%% Put data in FT format
disp("Converting to FieldTrip format")
[data_strct] = makeFTstruct(data_f_mat,fs,ch_table,S.sensor_info);

%% resample for viewing and mark artefacts
disp("Check for bad trials")
if ~exist([files.cleaning,'_vis_artfcts.mat'],'file')
    resamplefs = 150; %Hz
    cfg            = [];
    cfg.resamplefs = resamplefs;
    cfg.detrend    = 'no';
    data_preproc_150   = ft_resampledata(cfg, data_strct);
    
    %%% Get rid of bad trials
    cfg_art          = [];
    cfg_art.viewmode = 'vertical';
    cfg_art = ft_databrowser(cfg_art,data_preproc_150);
    clear data_preproc_150
    
    vis_artfcts = cfg_art.artfctdef.visual.artifact * ...
        (data_strct.fsample/resamplefs);
    save([files.cleaning,'_vis_artfcts.mat'],'vis_artfcts')
else
    load([files.cleaning,'_vis_artfcts.mat'],'vis_artfcts')
end

% automatic artifact rejection
thresh_val = 3;
auto_artfcts = get_bad_segments(data_f_mat,thresh_val);
fprintf('Found %d artifacts using a threshold of %d std. deviations.\n',...
    size(auto_artfcts,1),thresh_val);
clear data_f_mat
% combine artifacts and reject bad trials
cfg = [];
cfg.artfctdef.visual.artifact = [vis_artfcts;auto_artfcts];
cfg.artfctdef.reject  = 'complete';
data_vis_clean = ft_rejectartifact(cfg,data_strct);

fprintf('\nRejected %d of %d epochs of length %1.2f s.\n',...
    size(data_strct.trial,2)-size(data_vis_clean.trial,2),size(data_strct.trial,2),epoch_length);

clear data_strct
%% ICA
lay = Helmet_info.lay;
disp("ICA artifact rejection")
if ~exist([files.ICA,'_bad_ICA_comps.mat'],'file') || ~exist([files.ICA,'_ICA_data.mat'],'file')
    
    % Resample for faster ICA
    cfg            = [];
    cfg.resamplefs = 150;
    cfg.detrend    = 'no';
    data_ica_150   = ft_resampledata(cfg, data_vis_clean);
    
    % Run ICA on 150 Hz data or load previous unmixing matrix
    cfg            = [];
    if ~exist([files.ICA,'_ICA_data.mat'],'file')
        cfg.method = 'runica';
    else
        load([files.ICA,'_ICA_data.mat'],'comp150')
        cfg.unmixing   = comp150.unmixing;
        cfg.topolabel  = comp150.topolabel;
    end
    comp150    = ft_componentanalysis(cfg, data_ica_150);
    
    % Inspect components for rejection or load file with bad component list
    if ~exist([files.ICA,'_bad_ICA_comps.mat'],'file')
        disp("Choose bad components")
        [bad_comps] = plot_ICA_comps(comp150,ch_table,lay,[]);
        
        save([files.ICA,'_bad_ICA_comps.mat'],'bad_comps')
        close(gcf)
    else
        disp("Loading saved bad components")
        load([files.ICA,'_bad_ICA_comps.mat'],'bad_comps')
    end
    
    % only keep unmixing matrix and topolabel for component removal
    tokeep = {'unmixing','topolabel'};
    fns=fieldnames(comp150);
    toRemove = fns(~ismember(fns,tokeep));
    comp150 = rmfield(comp150,toRemove);
    save([files.ICA,'_ICA_data.mat'],'comp150')
    
    clear data_ica_150
else
    disp("Loading bad coponents, topographies and old unmixing matrix")
    load([files.ICA,'_bad_ICA_comps.mat'],'bad_comps')
    load([files.ICA,'_ICA_data.mat'],'comp150')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Perform ICA on original 1200 Hz data using unmixing matrix and topolabel
cfg            = [];
cfg.unmixing   = comp150.unmixing;
cfg.topolabel  = comp150.topolabel;
comp1200        = ft_componentanalysis(cfg, data_vis_clean);

% Plot comps again to confirm they are correct
disp("Confirm bad ICA components")
% [bad_comps] = plot_ICA_comps(comp1200,ch_table,lay,bad_comps)
save([files.ICA,'_bad_ICA_comps.mat'],'bad_comps');

% Remove components from data
cfg           = [];
cfg.component = bad_comps;
data_ica_clean    = ft_rejectcomponent(cfg, comp1200,data_vis_clean);

% Plot  example sensor time courses pre and post ica
figure(111);clf
plot(data_vis_clean.trial{1,1}(contains(data_vis_clean.label,'LR [Z]'),:),'Color',[0.5,0.5,0.5])
hold on
plot(data_ica_clean.trial{1,1}(contains(data_vis_clean.label,'LR [Z]'),:),'Color',[0.9,0.5,0.5])
xlabel('t/s')
legend('Pre ICA','Post ICA')

%% Reconstitute data from FT structure

data_f_clean = [data_ica_clean.trial{1,:}];

% apply mean field correction
disp("Applying mean field correction")
data_f_clean = S.M*data_f_clean;

%% Further steps:
%
% Filter to band of interest
% Get source model info and beamformer location etc
% Beamform... spit out VEs

%% Source positions
for n = 1:78
    sourcepos1 = ft_warp_apply(mri.transform,[voxlox(1,n) voxlox(2,n) voxlox(3,n)]);
    sourcepos(:,n) = sourcepos1'/1000; % convert to metres
    sourcepos(:,n)';
end

[bf_outs_shell] = run_beamformer('shell',sourcepos',S,0,[],1);
lead_fields_shell_xyz = bf_outs_shell.LF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert orientation of sources to polar
meshes = ft_convert_units(meshes,'m');

X = meshes.pnt; Origin = mean(X,1);
Ndips = length(sourcepos);
for n = 1:Ndips;
    thispos = sourcepos(:,n);
    [phi,theta1,r] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,thispos(3) - Origin(3));
    theta = pi/2 - theta1;
    Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
    Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
    Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
    Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make a plot of the geometry...
figure(1);
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
scatter3(sourcepos(1,:),sourcepos(2,:),sourcepos(3,:),'ro','linewidth',3)
view([130,30])
fig = gcf;
fig.Color = [1,1,1];
plot3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),'o')
quiver3(S.sensor_info.pos(ch_table.isx==1,1),S.sensor_info.pos(ch_table.isx==1,2),S.sensor_info.pos(ch_table.isx==1,3),...
    S.sensor_info.ors(ch_table.isx==1,1),S.sensor_info.ors(ch_table.isx==1,2),S.sensor_info.ors(ch_table.isx==1,3),'r','linewidth',2)
quiver3(S.sensor_info.pos(ch_table.isy==1,1),S.sensor_info.pos(ch_table.isy==1,2),S.sensor_info.pos(ch_table.isy==1,3),...
    S.sensor_info.ors(ch_table.isy==1,1),S.sensor_info.ors(ch_table.isy==1,2),S.sensor_info.ors(ch_table.isy==1,3),'g','linewidth',2)
quiver3(S.sensor_info.pos(ch_table.isz==1,1),S.sensor_info.pos(ch_table.isz==1,2),S.sensor_info.pos(ch_table.isz==1,3),...
    S.sensor_info.ors(ch_table.isz==1,1),S.sensor_info.ors(ch_table.isz==1,2),S.sensor_info.ors(ch_table.isz==1,3),'b','linewidth',2)
plot3(Origin(1),Origin(2),Origin(3),'bo','linewidth',4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% take a random lead field and plot it...
figure(2);
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
ft_plot_topo3d(double(S.sensor_info.pos(ch_table.isz==1,:)),Lead_fields(ch_table.isz==1,1,16))
alpha(gca,0.5)
plot3(S.sensor_info.pos(ch_table.isx==1,1),S.sensor_info.pos(ch_table.isx==1,2),S.sensor_info.pos(ch_table.isx==1,3),'go','linewidth',3)
scatter3(sourcepos(1,16),sourcepos(2,16),sourcepos(3,16),'r','linewidth',4)
% quiver3(S.sensor_info.pos(ch_table.isx==1,1),S.sensor_info.pos(ch_table.isx==1,2),S.sensor_info.pos(ch_table.isx==1,3),...
%     S.sensor_info.ors(ch_table.isx==1,1).*Lead_fields(ch_table.isx==1,2,16) + S.sensor_info.ors(ch_table.isy==1,1).*Lead_fields(ch_table.isy==1,2,16) + S.sensor_info.ors(ch_table.isz==1,1).*Lead_fields(ch_table.isz==1,2,16),...
%     S.sensor_info.ors(ch_table.isx==1,2).*Lead_fields(ch_table.isx==1,2,16) + S.sensor_info.ors(ch_table.isy==1,2).*Lead_fields(ch_table.isy==1,2,16) + S.sensor_info.ors(ch_table.isz==1,2).*Lead_fields(ch_table.isz==1,2,16),...
%     S.sensor_info.ors(ch_table.isx==1,3).*Lead_fields(ch_table.isx==1,2,16) + S.sensor_info.ors(ch_table.isy==1,3).*Lead_fields(ch_table.isy==1,2,16) + S.sensor_info.ors(ch_table.isz==1,3).*Lead_fields(ch_table.isz==1,2,16),'r','linewidth',2)
quiver3(S.sensor_info.pos(ch_table.isz==1,1),S.sensor_info.pos(ch_table.isz==1,2),S.sensor_info.pos(ch_table.isz==1,3),...
    S.sensor_info.ors(ch_table.isz==1,1).*Lead_fields(ch_table.isz==1,2,16),...
    S.sensor_info.ors(ch_table.isz==1,2).*Lead_fields(ch_table.isz==1,2,16),...
    S.sensor_info.ors(ch_table.isz==1,3).*Lead_fields(ch_table.isz==1,2,16),'r','linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except path files sub ses epoch_length data_ica_clean data_f_clean cleaning_only fs sourcepos Lead_fields Ndips hp lp
%%
%% filter the OPM data to band of interest
if ~cleaning_only
    
    hp;
    lp;
    if ~exist(sprintf('%s%s_%d_%d_Hz_Z.mat',path.VEs,files.VEs,hp,lp),'file') 
    data_f = data_f_clean';
    %% Beamform
    C = cov(data_f);
    mu = 0.05;
    Cr = C + mu*max(svd(C))*eye(size(C));
    Cr_inv = inv(Cr);
    for n = 1:Ndips
        this_L = Lead_fields(:,:,n);
        W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
        iPower_v = this_L'*Cr_inv*this_L;
        [v,d] = svd(iPower_v);
        [~,id] = min(diag(d));
        lopt = this_L*v(:,id); % turn to nAm amplitude
        w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
        VE(:,n) = (w*data_f'./sqrt(w*w'));
    end
    save(sprintf('%s%s_%d_%d_Hz_Z.mat',path.VEs,files.VEs,hp,lp),'VE')

    VE = (VE - mean(VE,1))./std(VE,[],1);
    save(sprintf('%s%s_%d_%d_Hz_Z_standard.mat',path.VEs,files.VEs,hp,lp),'VE')
    end
    load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.VEs,files.VEs,hp,lp),'VE')
    
    hpfs = [2,4, 8,13,30,35,40];
    lpfs = [4,8,12,30,40,45,48];
    
    N_fs = length(hpfs);
    for f_i = 1:N_fs
        hp = hpfs(f_i);
        lp = lpfs(f_i);
        [b,a] = butter(4,2*[hp lp]/fs);
        VE_f = [filtfilt(b,a,VE)];
        
        hilb_envs = abs(hilbert(VE_f));
        ENV_std_over_mean = std(hilb_envs,[],1)./mean(hilb_envs,1);
        filt_vars= var(VE_f,[],1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        fprintf('Sub: %s | Session: %s | Freq %3d/%3d (%3d-%3d Hz) | pow and env dynamics\n',sub,ses,f_i,length(hpfs),hp,lp);
        save(sprintf('%s%s_%d_%d_Hz_Z.mat',path.pow,files.pow,hp,lp),'filt_vars','ENV_std_over_mean')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    load(sprintf('%s%s_%d_%d_Hz_Z_standard.mat',path.VEs,files.VEs,hp,lp),'VE')
 
    N_fs = length(hpfs);
    for f_i = 1:N_fs
        hp = hpfs(f_i);
        lp = lpfs(f_i);
        [b,a] = butter(4,2*[hp lp]/fs);
        VE_f = [filtfilt(b,a,VE)];
        
        hilb_envs = abs(hilbert(VE_f));
        ENV_std_over_mean = std(hilb_envs,[],1)./mean(hilb_envs,1);
        filt_vars= var(VE_f,[],1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        fprintf('Sub: %s | Session: %s | Freq %3d/%3d (%3d-%3d Hz) | pow and env dynamics\n',sub,ses,f_i,length(hpfs),hp,lp);
        save(sprintf('%s%s_%d_%d_Hz_Z.mat',path.pow,files.pow,hp,lp),'filt_vars','ENV_std_over_mean')
    end
end
end
