%% Data Cleaning Test

function data_cleaning_test_noICA(sub,ses)
restoredefaultpath
cleaning_only = 0;
close all
clc

% sub = '004';
% ses = '001';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_dir = 'R:\OPMMEG\Projects\movie\';
addpath([project_dir,'scripts\fieldtrip-20190212'])
addpath([project_dir,'scripts'])
addpath([project_dir,'scripts\Beamformer\'])
ft_defaults;
datadir = [project_dir,'data\'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_type = 'task-movie';
filename = ['sub-',sub,'_',exp_type,'_','run-',ses];
path.main = [datadir,'sub-',sub,'\'];
path.ICA = [datadir,'derivatives\ICA\sub-',sub,'\'];
path.cleaning = [datadir,'derivatives\cleaning\sub-',sub,'\'];
path.meshes = [datadir,'derivatives\sourcespace\sub-',sub,'\'];
path.VEs = [datadir,'derivatives\VEs\sub-',sub,'\'];
path.AEC = [datadir,'derivatives\AEC\sub-',sub,'\'];

if ~exist(path.ICA,'dir'); mkdir(path.ICA);end
if ~exist(path.cleaning,'dir'); mkdir(path.cleaning);end
if ~exist(path.VEs,'dir'); mkdir(path.VEs);end
if ~exist(path.AEC,'dir'); mkdir(path.AEC);end

path.data = [path.main,'meg\'];
path.mri = [path.main,'anat\'];
files.mri = ['sub-',sub,'.nii'];
files.meshes = ['sub-',sub,'_meshes.mat'];
files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
files.AEC = ['sub-',sub,'_','run-',ses,'_AEC'];
files.voxlox = ['sub-',sub,'_voxlox.mat'];
files.channels = '_channels.tsv';
files.sens_order = '_sensor_order.mat';
files.ICA = [path.ICA,filename];
files.cleaning = [path.cleaning,filename];
S.mri_file = [path.meshes,'sub-',sub,'_\x.nii']; % only need path for beamformer and single shell function. bit hacky, FIX
% global files
cd(path.data)
hpfs = [4,8,13,30,50,1];
lpfs = [8,12,30,50,100,100];

for f_i = 1:length(hpfs)
load([path.data,filename,'_meg.mat'],'fs','data','triggers')
load([datadir,'derivatives\helmet\M1P5_adult.mat'])
% basefilename = [path.main,'ses-',ses,'\Matt\'];
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


    
    % hp = 13;
    % lp = 30;
    hp = hpfs(f_i);
    lp = lpfs(f_i);
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
    
    load([path.data,filename,files.sens_order],'T')
    for sl_i = 1:height(T)
        if ~isempty(T{sl_i,1}{1})
            ch_table.slot_no(startsWith(ch_table.name,T{sl_i,1}{1}(1:2))) = sl_i;
        end
    end
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
    
    % combine artifacts and reject bad trials
    cfg = [];
    cfg.artfctdef.visual.artifact = [vis_artfcts;auto_artfcts];
    cfg.artfctdef.reject  = 'complete';
    data_vis_clean = ft_rejectartifact(cfg,data_strct);
    
    fprintf('\nRejected %d of %d epochs of length %1.2f s.\n',...
        size(data_strct.trial,2)-size(data_vis_clean.trial,2),size(data_strct.trial,2),epoch_length);
    
    
    %% Reconstitute data from FT structure
    
    data_f_clean = [data_vis_clean.trial{1,:}];
    
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
        sourcepos(:,n)'
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
    ft_plot_topo3d(double(S.sensor_info.pos(ch_table.isz==1,:)),Lead_fields(ch_table.isz==1,2,16))
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
    %%
    %% filter the OPM data to band of interest
    
    if ~exist(sprintf('%s%s_%d_%d_Hz_noICA.mat',path.AEC,files.AEC,hp,lp),'file')
        
        data_f = data_f_clean';
        duration = epoch_length*size(data_vis_clean.trial,2);
        
        %% Beamform
        C = cov(data_f);
        mu = 0.05;
        Cr = C + mu*max(svd(C))*eye(size(C));
        Cr_inv = inv(Cr);
        for n = 1:length(sourcepos)
            this_L = Lead_fields(:,:,n);
            W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
            iPower_v = this_L'*Cr_inv*this_L;
            [v,d] = svd(iPower_v);
            [~,id] = min(diag(d));
            lopt = this_L*v(:,id); % turn to nAm amplitude
            w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
            VE(:,n) = w*data_f';
        end
        save(sprintf('%s%s_%d_%d_Hz_noICA.mat',path.VEs,files.VEs,hp,lp),'VE')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        down_f = 10;
        Nlocs = 78;
        for seed = 1:78;
            Xsig = squeeze(VE(:,seed));
            for test = 1:78;
                if seed == test
                    AEC(seed,test) = NaN;
                else
                    Ysig = squeeze(VE(:,test));
                    X_win = (Xsig - mean(Xsig));
                    Y_win = (Ysig - mean(Ysig));
                    %%regress leakage
                    beta_leak = (pinv(X_win)*Y_win);
                    Y_win_cor = Y_win - X_win*beta_leak;
                    %%calculate envelopes
                    H_X = abs(hilbert(X_win));
                    H_X_d = mean(reshape(H_X,fs/down_f,duration*down_f,1));
                    %%calculate envelopes
                    H_Y = abs(hilbert(Y_win_cor));
                    H_Y_d = mean(reshape(H_Y,fs/down_f,duration*down_f,1));
                    AEC(seed,test) = corr(H_X_d',H_Y_d');
                end
            end
            clc;disp(sprintf('running beta connectivity region %d',seed));
        end
        AEC = 0.5*(AEC + AEC');
        figure()
        subplot(121)
        imagesc(AEC);colorbar;
        subplot(122)
        go_netviewer_Matt(AEC,0.7)
        save(sprintf('%s%s_%d_%d_Hz_noICA.mat',path.AEC,files.AEC,hp,lp),'AEC')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
end