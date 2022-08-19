%% Run beamformer after exporting data from Analysis GUI
clearvars -except OPM_vars
close all
clc
restoredefaultpath
warning off
p = matlab.desktop.editor.getActiveFilename;
[p1,p2] = fileparts(p);
p3 = strsplit(p1,filesep);
p4 = strjoin(p3(1:5),filesep);
addpath([p4 '\MatLab files'])
addpath([p4 '\MatLab files\tools'])
addpath([p4 '\MatLab files\analysis_gui\app_functions']);
set(0,'DefaultFigureWindowStyle','docked')
if isempty(which('ft_defaults'))
    addpath([p4 '\MatLab files\fieldtrip-20190212'])
    ft_defaults
end
addpath([p4 '\MatLab files\fieldtrip-20190212\external\spm12'])

%% Inputs
% Load participant information
inputs = inputdlg({'Enter volunteer number:','Date','Filename QZFM_#'},...
    'Inputs',[1 35],{'02524','20210121','3'});
vol_nr      = inputs{1}; % subject ID
filename    = ['QZFM_',inputs{3}];
date        = inputs{2};

% Skanect
path.skanect = [p4 '\data\skanect\' date '\' vol_nr '\'];
files.skanect = ['sens_info_' filename '.mat'];
% MRI
path.mri = [p4 '\data\MRI\' vol_nr '\'];
files.mri = [vol_nr '.mri'];
files.meshes = 'meshes.mat';
% QZFM data
path.data = [p4 '\data\QZFM\' date '\'];
files.data = [filename '.tdms'];

% Results folder
path.results = [path.data,files.data(1:end-5),'_outputs\'];
if ~exist(path.results)
    mkdir(path.results)
end
cd(path.data)
if ~exist(['OPM_vars_',inputs{3},'.mat'])
    fprintf(' 1) Run Analysis GUI \n 2) Remove bad trials, bad channels, etc ... \n 3) Export data \n');
    cd([p4 '\MatLab files\analysis_gui'])
    run('analysis_GUI.mlapp')
else
    load(['OPM_vars_',inputs{3},'.mat'])
end

%% Data in structure
OPM_data_struct = OPM_vars.OPM_data_struct;
QZFM_data = OPM_vars.QZFM_data;

%% Generate meshes from fieldtrip, also do MRI segmentation
cd(path.mri)
if ~exist(files.meshes)
    cfg = [];
    cfg.method = 'fieldtrip';
    [meshes, segmentedmri] = go_MRI2Meshes(cfg,files.mri); % meshes output in metres
    save([path.mri,files.meshes],'meshes','segmentedmri')
else
    load(files.meshes)
end

%% Read MRI, extract brain using the segmented mri output, downsample to 4mm
path.anat = [path.mri '\ctf_space\'];
files.anat = 'anat.nii';
files.brain = 'brain.nii';
files.brainDS = 'brain_4mm.nii';
if ~exist([path.anat files.anat]) || ...
        ~exist([path.anat files.brain]) || ...
        ~exist([path.anat files.brainDS])
    disp('Creating new .nii files...')
    mri = ft_read_mri([path.mri files.mri]);
    
    brain = mri;
    brain.anatomy = brain.anatomy.*segmentedmri.brain;
    
    cfg = [];
    cfg.parameter = 'anatomy';
    cfg.downsample = 4;
    brainDS = ft_volumedownsample(cfg,brain);
    
    disp('Writing .nii files...')
    if ~exist([path.mri 'ctf_space\'])
        mkdir([path.mri 'ctf_space\'])
    end
    ft_write_mri([path.anat files.anat],mri.anatomy,'dataformat','nifti','transform',mri.transform);
    ft_write_mri([path.anat files.brain],brain.anatomy,'dataformat','nifti','transform',brain.transform);
    ft_write_mri([path.anat files.brainDS],brainDS.anatomy,'dataformat','nifti','transform',brainDS.transform);
    disp('Done!')
else
%     addpath([p4 '\MatLab files\fieldtrip-20190212\external\spm12\external\freesurfer'])
    disp('Loading .nii files...')
    mri = ft_read_mri([path.mri files.mri]);
    anat = ft_read_mri([path.anat files.anat]);
    brain = ft_read_mri([path.anat files.brain]);
    brainDS = ft_read_mri([path.anat files.brainDS]);
    disp('Done!')
end

%% Generate the source locations
% Make brain mask
brainDS.inside = brainDS.anatomy>0;
% Find voxels
voxID = find(brainDS.inside);
[x,y,z] = ind2sub(brainDS.dim,voxID);
% Generate CTF coorinates for them.
sourcepos = ft_warp_apply(brainDS.transform,[x y z]);
sourcepos = sourcepos/1000; % convert to metres
f1 = figure(1);
f1.Name = sprintf('Coverage of %s',num2str(vol_nr));
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3))
view([180,0])
fig = gcf;
fig.Color = [1,1,1];

%% Load in OPM information
Nlocs = OPM_vars.Nlocs;
loc_names = OPM_vars.loc_names;
pos_no = OPM_vars.pos_no;

% Skanect files
load([path.skanect files.skanect])
[~,filename_skanect] = fileparts(files.skanect);
% Sensor channels
if ~isfield(sens_info,'pos_used')
    if ~exist([path.data files.data(1:end-5) '_sensor_order.mat'])
        sensorgui_all
        save([path.data files.data(1:end-5) '_sensor_order.mat'],'T')
    else
        load([path.data files.data(1:end-5) '_sensor_order.mat'])
    end
    
    for n = 1:Nlocs
        Loc_ind = OPM_data_struct.(['Sensor_',loc_names{n}]);
        this_axes = QZFM_data.axis(Loc_ind);
        for nn = 1:length(Loc_ind)
            sens_info.pos_used(Loc_ind(nn),:) = sens_info.pos(pos_no(n),:);
            if strcmpi(this_axes{nn},'y')
                sens_info.ors_used(Loc_ind(nn),:) = sens_info.ors_tan(pos_no(n),:);
%                 sens_info.ors_tan_used(n,:) = sens_info.ors_tan(pos_no(n),:);
            elseif strcmpi(this_axes{nn},'z')
                sens_info.ors_used(Loc_ind(nn),:) = sens_info.ors(pos_no(n),:);
%                 sens_info.ors_rad_used(n,:) = sens_info.ors(pos_no(n),:);
%                 sens_info.pos_rad_used(n,:) = sens_info.pos(pos_no(n),:);
            end
        end
        sens_info.sensornamesinuse{n} = loc_names{n};
    end
    save([path.skanect files.skanect],'sens_info','-append')
end

%% 
% *************************************************************************
% **** Rotate Gen-3 sensors' tangential orientation 180 degrees ***********
% Should only do it if Gen 3s are present now, and will check it is in Y
% mode before flipping anything
sens_flip = {'AA','L1','L2'};
for n = 1:length(sens_flip)
    if ismember(sens_flip{n},loc_names)
        Loc_ind = OPM_data_struct.(['Sensor_',sens_flip{n}]);
        if strcmpi(QZFM_data.axis{Loc_ind(1)},'Y')
            tan = sens_info.ors_used(Loc_ind(1),:);
            sens_info.ors_used(Loc_ind(1),:) = tan*-1;
        end
    end
end
% *************************************************************************
% **** Flip ALL sensors orientations **************************************
% sens_flip = loc_names;
% for n = 1:length(sens_flip)
%     Loc_ind = OPM_data_struct.(['Sensor_',sens_flip{n}]);
%     tan = sens_info.ors_used(Loc_ind(1),:);
%     rad = sens_info.ors_used(Loc_ind(2),:);
%     sens_info.ors_used(Loc_ind(1),:) = tan*-1;
%     sens_info.ors_used(Loc_ind(2),:) = rad*-1;
% end
% *************************************************************************
%%
figure(1)
scatter3(sens_info.pos_used(:,1),sens_info.pos_used(:,2),sens_info.pos_used(:,3))
text(unique(sens_info.pos_used(:,1),'stable'),unique(sens_info.pos_used(:,2),'stable'),...
    unique(sens_info.pos_used(:,3),'stable'),...
    sens_info.sensornamesinuse,'color','k')
hold on; axis equal;
quiver3(sens_info.pos_used(:,1),sens_info.pos_used(:,2),sens_info.pos_used(:,3),...
    sens_info.ors_used(:,1),sens_info.ors_used(:,2),sens_info.ors_used(:,3))
fig = gcf;
fig.Color = [1,1,1];

% Layout
layout_name = QZFM_data.Layout;
load(strcat('R:\Matlab_Files\fieldtrip-20190212\template\layout\',layout_name))
lay.label = QZFM_data.Sensor_info.Name;
f2 = figure;    
f2.Name = 'Sensor layout';
ft_plot_lay(lay,'box','no')
fig = gcf;
fig.Color = [1,1,1];

%% Extract OPM data
QZFM_data = OPM_vars.QZFM_data;
% Unchopped data
OPM_data1 = QZFM_data.data;
% Chopped data and concatenated
OPM_data = OPM_vars.OPM_data;
OPM_data_f = OPM_vars.OPM_dataf;
% Time
time = OPM_vars.time;
% Trial time vector
trial_time = OPM_vars.trial_time;
% Sampling frequency (Hz)
f = floor(OPM_vars.f);
% Number of channels
Nchans = OPM_vars.Nchans;
sensornamesinuse = OPM_vars.sensornamesinuse;
% Trigger channel
trigger = OPM_vars.trigger;
% Duration of trial
duration = OPM_vars.duration;
% Trigger offset (s)
trig_offset = OPM_vars.trig_offset;
% Find spikes
ind = OPM_vars.ind;
% Number of trials
Ntrials = OPM_vars.Ntrials;
% Plot trigger channel
figure(123)
title('Trigger channel(s)')
plot(QZFM_data.time,trigger)
hold on
plot(QZFM_data.time(ind),ones(length(ind)),'o')
xlabel('Time (s)');ylabel('Voltage (V)')
title([num2str(Ntrials) ' trials found'])

% Filter raw data
hp = OPM_vars.hp;
lp = OPM_vars.lp;
OPM_data1f = nut_filter3(OPM_data1,'butter','bp',4,hp,lp,f,1);

% Plot
f2 = figure;
f2.Name = 'All channels';
plot(time,OPM_data_f+(5e3*[0:1:Nchans-1]))
xlim([time(1), time(end)])
ylim([-5e3,5e3*Nchans])
f2.CurrentAxes.YTick = 0:5e3:5e3*Nchans-1;
f2.CurrentAxes.YTickLabel = sensornamesinuse;
fig = gcf;
fig.Color = [1,1,1];
xlabel('Time (s)');ylabel('Channels')
title(sprintf('All channels - Filtered data [%u-%u] Hz',hp,lp))
clc

% Plot on layout
ylimit = 500; % y axis limit (pT)
f3 = figure;
f3.Name = 'Raw data';
p = uipanel('Title','Axes - Raw data','Units','normalized','Position',...
    [0.85 0.85 0.1 0.12],'BackgroundColor','w','BorderType','none');
ax = axes(p);
grid on
xlabel('Time (s)');ylabel('Field (pT)')
ylim([-ylimit ylimit])
xlim([min(time) max(time)])
for n = 1:Nlocs
    Loc_ind = OPM_data_struct.(['Sensor_',loc_names{n}]);
    p = uipanel('Title',loc_names(n),'Units','normalized','Position',...
        [0.5+lay.pos(pos_no(n),1) (0.48+lay.pos(pos_no(n),2))*0.95 0.07 0.09],...
        'BorderType','none','BackgroundColor','w');
    ax = axes(p);
    plot(time,OPM_data(:,Loc_ind)/1e3) % in pT
    hold on
    xlim([time(1) time(end)])
    ylim([-ylimit ylimit])
    grid on
    ax.XTickLabel = {};ax.YTickLabel = {};
end
fig = gcf;
fig.Color = [1,1,1];

%% TFS
do_TFS = 0;

highpass = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110];
lowpass = [4 6 8 10 13 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
fre = highpass + ((lowpass - highpass)./2);

% Control window
conwin = [4 5].*f-trig_offset*f;

if do_TFS
    figure(200);
    p = uipanel('Title','Axes','Units','normalized','Position',...
        [0.85 0.85 0.1 0.12],'BackgroundColor','w','BorderType','none');
    ax = axes(p);
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    ylim([min(fre) max(fre)])
    xlim([min(trial_time) max(trial_time)])
    axis on
    colorbar;caxis([-.4 .4])
    set(gcf,'color',[1 1 1])
    
    figure(201);
    p = uipanel('Title','Axes','Units','normalized','Position',...
        [0.85 0.85 0.1 0.12],'BackgroundColor','w','BorderType','none');
    ax = axes(p);
    xlabel('Time (s)');ylabel('Frequency (Hz)')
    ylim([min(fre) max(fre)])
    xlim([min(trial_time) max(trial_time)])
    axis on
    colorbar;caxis([-.4 .4])
    set(gcf,'color',[1 1 1])
    
    for n = 1:Nlocs
        Loc_ind = OPM_data_struct.(['Sensor_',loc_names{n}]);
        this_axes = QZFM_data.axis(Loc_ind);
        for nn = 1:length(Loc_ind)
            if strcmpi(this_axes{nn},'y')
                fg1 = figure(200);
                fg1.Name = 'TFS for Y axis';
            elseif strcmpi(this_axes{nn},'z')
                fg2 = figure(201);
                fg2.Name = 'TFS for Z axis';
            end
            % This channel data
            OPM_ch = OPM_data1(:,Loc_ind(nn))';
            % Filter data within bands and calculate envelope
            OPM_ch_fb = zeros(length(OPM_ch),length(fre));
            for fb = 1:length(highpass)
                fprintf('\n Band %u/%u ',fb,length(highpass))
                filt_VE = nut_filter3(OPM_ch','butter','bp',3,highpass(fb),lowpass(fb),f,1)';
                OPM_ch_fb(:,fb) = abs(hilbert(filt_VE));
            end
            OPM_mean = zeros(duration*f,length(fre));
            for fb = 1:length(highpass)
                OPM_ch_fb_trials = [];
                % Chop data
                for i = 1:Ntrials
                    OPM_ch_fb_trials = cat(1,OPM_ch_fb_trials,OPM_ch_fb(ind(i):ind(i)+(duration*f-1),fb));
                end
                OPM_ch_filt = reshape(OPM_ch_fb_trials,duration*f,Ntrials);
                
                % Average across trials
                OPM_mean(:,fb) = mean(OPM_ch_filt,2);
            end
            meanrest = mean(OPM_mean(conwin(1):conwin(2),:),1);
            meanrestmat = repmat(meanrest,size(OPM_mean,1),1);
            TFS = (OPM_mean'-meanrestmat')./meanrestmat';
            % Plot
            p = uipanel('Title',sensornamesinuse(Loc_ind(nn)),...
                'Units','normalized','Position',...
                [0.5+lay.pos(pos_no(n),1) (0.48+lay.pos(pos_no(n),2))*0.95 0.07 0.09],...
                'BackgroundColor','w','BorderType','none');
            ax = axes(p);
            
            pcolor(trial_time,fre,TFS);shading interp
            axis fill
            caxis([-.4 .4])
            axis off
            title(sprintf('Ch %s',num2str(Loc_ind(nn))))
            pause(0.1)
        end
    end
end

%% Plot single TFS over active region to determine on/off windows
% Select manually
ch1 = find(ismember(sensornamesinuse,'L2 (Z)'));
OPM_ch = OPM_data1(:,ch1)';
OPM_ch_fb = zeros(length(OPM_ch),length(fre));
for fb = 1:length(highpass)
    fprintf('\n Band %u/%u ',fb,length(highpass))
    filt_VE = nut_filter3(OPM_ch','butter','bp',4,highpass(fb),lowpass(fb),f,1)';
    OPM_ch_fb(:,fb) = abs(hilbert(filt_VE));
end
OPM_mean = zeros(duration*f,length(fre));
for fb = 1:length(highpass)
    OPM_ch_fb_trials = [];
    % Chop data
    for i = 1:Ntrials
        OPM_ch_fb_trials = cat(1,OPM_ch_fb_trials,OPM_ch_fb(ind(i)+1:ind(i)+(duration*f),fb));
    end
    OPM_ch_filt = reshape(OPM_ch_fb_trials,duration*f,Ntrials);
    OPM_mean(:,fb) = mean(OPM_ch_filt,2);
end
meanrest = mean(OPM_mean(conwin(1):conwin(2),:),1);
meanrestmat = repmat(meanrest,size(OPM_mean,1),1);
TFS = (OPM_mean'-meanrestmat')./meanrestmat';

% plot
f4 = figure;
f4.Name = 'TFS for selected sensor';
pcolor(trial_time,fre,TFS);shading interp
axis square;
title(strcat(sensornamesinuse{ch1}))
set(gca,'fontsize',14);set(gcf,'color',[1 1 1]);
xlabel('Time (s)');ylabel('Frequency (Hz)')
cb = colorbar; cb.Label.String = 'Relative change from baseline';
caxis([-0.5 0.5])
set(gcf,'color',[1 1 1]);


% Hilbert envelope
OPM_ch_f = nut_filter3(OPM_ch','butter','bp',4,10,20,f,1)';
H_OPM_ch_f = abs(hilbert(OPM_ch_f));
% Chop data
H_OPM_ch_f_trials = [];
for i = 1:Ntrials
    H_OPM_ch_f_trials = cat(1,H_OPM_ch_f_trials,H_OPM_ch_f(ind(i)+1:ind(i)+(duration*f)));
end
mean_H_OPM_ch1 = mean(H_OPM_ch_f_trials,1);
mean_H_OPM_ch = (mean_H_OPM_ch1 - mean(mean_H_OPM_ch1(conwin(1):conwin(2))))./mean(mean_H_OPM_ch1(conwin(1):conwin(2)))*100;
% Plot
fg = figure;
fg.Name = 'Oscillatory amplitude';
plot(trial_time,mean_H_OPM_ch)
grid on
xlabel('Time (s)');ylabel('Relative change (%)')
set(gcf,'color',[1 1 1]);
hold on

ON_winds = inputdlg('What is the On window (in secconds)? ');
ON_winds = str2num(cell2mat(ON_winds));
OFF_winds = inputdlg('What is the Off window (in secconds)? ');
OFF_winds = str2num(cell2mat(OFF_winds));
plot(trial_time(ON_winds(1)*f:ON_winds(2)*f),mean_H_OPM_ch(ON_winds(1)*f:ON_winds(2)*f),'r')
plot(trial_time(OFF_winds(1)*f:OFF_winds(2)*f),mean_H_OPM_ch(OFF_winds(1)*f:OFF_winds(2)*f),'r')

%% Filter between different band
% OPM_data_f1 = nut_filter3(OPM_data1,'butter','bp',4,13,30,f,1);
% % Chop data into trials
% OPM_data_f = [];
% for n = 1:Ntrials
%     OPM_data_f = cat(1,OPM_data_f,OPM_data_f1(ind(n):ind(n)+(duration*f-1),:));
% end

%% Make covariance matrix and regularise
% Regularisation parameter (%)
mu = 0.05;
% Data in teslas
OPM_data_T = OPM_data'*1e-15;
OPM_data_f_T = OPM_data_f'*1e-15;
OPM_data_T_mat = reshape(OPM_data_T, [Nchans,duration*f,Ntrials]);
OPM_data_f_T_mat = reshape(OPM_data_f_T, [Nchans,duration*f,Ntrials]);
% Covariance matrix
C = cov(OPM_data_f_T');
maxEV = max(svd(C));
Noise_Cr = min(svd(C)).*eye(size(C));
C = C + mu.*maxEV.*eye(size(C));
condnr = cond(C);
Cinv = inv(C);
f5 = figure;
f5.Name = 'Covariance matrix';
imagesc(C)
axis square
title(sprintf('\n mu = %.3f, cond(C) = %.2f',...
    mu,condnr))
% Active and Control covariance
ONwin   = [ON_winds]-trig_offset;
OFFwin  = [OFF_winds]-trig_offset;
OPM_data_f_T_mat_A = OPM_data_f_T_mat(:,ONwin(1)*f+1:ONwin(2)*f,:);
OPM_data_f_T_mat_C = OPM_data_f_T_mat(:,OFFwin(1)*f+1:OFFwin(2)*f,:);
OPM_data_f_T_A = reshape(OPM_data_f_T_mat_A,Nchans,[]);
OPM_data_f_T_C = reshape(OPM_data_f_T_mat_C,Nchans,[]);
Ca = cov(OPM_data_f_T_A');
Cc = cov(OPM_data_f_T_C');

%% Beamforming
fwd_method = input('Single sphere (quickest), local spheres (medium), Single Shell, or BEM (slowest but more accurate?)? (Single, Local, BEM, Shell):   ','s') ;
% Single sphere
if strcmpi(fwd_method,'single')
    origin = mean(sourcepos);
    global mean_sphere_origin
    mean_sphere_origin = origin; % (m)  
    BF_pos = sourcepos;
    T1 = cell(length(BF_pos),1);
    h = waitbar(0,'Running Beamformer (Working...)');
    for ii = 1:length(BF_pos)
        % compute the beamformer projected timecourse using the weights
        ctfmeg_coord = BF_pos(ii,:).*100; % (cm)
        W1 = Quick_beamformer_target_OPM(C,Cinv,Noise_Cr,...
            ctfmeg_coord(1),ctfmeg_coord(2),ctfmeg_coord(3),...
            sens_info.pos_used*100,sens_info.ors_used);
        % Get the stat
        Qa = W1'*Ca*W1;
        Qc = W1'*Cc*W1;
        T1{ii} = (Qa-Qc)./(2*Qc);
        waitbar(ii/length(BF_pos),h)
    end
    delete(h)
    % Save T stat
    files.results = ['T_',num2str(hp),'-',num2str(lp),'_dual_SS_VS'];
    if ~exist([path.results,files.results,'.mat'])
        save([path.results,files.results,'.mat'],'mu','condnr','T1')
    else
        save([path.results,files.results,'.mat'],'mu','condnr','T1','-append')
    end
% Local spheres
elseif strcmpi(fwd_method,'local')
    % Create local spheres
    addpath([p4 '\MatLab files\Local_Spheres\']);
    if exist([path.skanect,'local_spheres_info_',files.data(1:end-5),'.mat'])
        load([path.skanect,'local_spheres_info_',files.data(1:end-5),'.mat']);
    else
        brain_mesh.faces = meshes(1).tri;
        brain_mesh.vertices = meshes(1).pnt;
        for n = 1:Nlocs
            Loc_ind = OPM_data_struct.(['Sensor_',loc_names{n}]);
            this_axes = QZFM_data.axis(Loc_ind);
            for nn = 1:length(Loc_ind)
                if strcmpi(this_axes{nn},'y')
                % Always use radial orientations for local sphere calculation
                elseif strcmpi(this_axes{nn},'z')
                    sensors.pos(n,:) = sens_info.pos(pos_no(n),:);
                    sensors.ors(n,:) = sens_info.ors(pos_no(n),:);
                end
            end
        end
        [so, sr] = create_local_spheres(brain_mesh,sensors);
        save([path.skanect,'local_spheres_info_',files.data(1:end-5),'.mat'],...
            'so','sr','brain_mesh','sensors');
    end
    % Plot spheres
    figure
    ft_plot_mesh(meshes(1),'facecolor',[.5 .5 .5],'facealpha',0.3,'edgecolor','none')
    hold on
    [sx,sy,sz] = sphere;
    [origin,rad_sing] = sphereFit(brain_mesh.vertices);
    sx = sx*rad_sing + origin(1);
    sy = sy*rad_sing + origin(2);
    sz = sz*rad_sing + origin(3);
    surf(sx,sy,sz,'Edgecolor','none','Facecolor',[0.2 0.2 0.2],'FaceAlpha',0.1);
    axis equal
    s1 = [];
    s = [];
    ccc = jet(size(sensors.pos,1));
    for n = 1:size(sensors.pos,1)
        s = scatter3(sensors.pos(n,1),sensors.pos(n,2),sensors.pos(n,3),'MarkerFaceColor',[ccc(n,:)]);
        [sx,sy,sz] = sphere;
        sx = sx*sr(n,:) + so(n,1);
        sy = sy*sr(n,:) + so(n,2);
        sz = sz*sr(n,:) + so(n,3);
        s1 = surf(sx,sy,sz,'Edgecolor','none','Facecolor',[ccc(n,:)],'FaceAlpha',0.3);
    end
    % duplicate the sphere origins for Dual mode
    for n = 1:Nlocs
        Loc_ind = OPM_data_struct.(['Sensor_',loc_names{n}]);
        this_axes = QZFM_data.axis(Loc_ind);
        for nn = 1:length(Loc_ind)
            origins(Loc_ind(nn),:) = so(n,:);
        end
    end
    origin_mean = mean(origins); % Mean of local sphere origins
    
    % Beamforming
    addpath([p4 '\Matlab_files\tools\']);
    vn = patchnormals(brain_mesh);
    BF_pos = sourcepos;
    T1 = cell(length(BF_pos),1);
    S.pos = sens_info.pos_used;
    S.ors = sens_info.ors_used;
    S.C = C;
    S.Noise_Cr = Noise_Cr;
    S.Cinv = Cinv;
    S.ls_origins = origins;
    S.tan_method = 'norm';
    h = waitbar(0,'Running Beamformer (Working...)');
    for ii = 1:length(BF_pos)
        % Find nearest surface normal if using it for tangential calculation
        source2mesh_dist = pdist2(BF_pos(ii,:),brain_mesh.vertices);
        [~,s2m_idx] = min(source2mesh_dist);
        nrst_norm = vn(s2m_idx,:); % Nearest surface normal to the dipole location 
        S.dip_loc = BF_pos(ii,:);
        S.tan_ref = nrst_norm;
        % Find lead fields and compute BF weights
        [Outs] = BF_local_spheres_target(S);
        % Get the stat
        W1 = Outs.Weights;
        Qa = W1'*Ca*W1;
        Qc = W1'*Cc*W1;
        T1{ii} = (Qa-Qc)./(2*Qc);
%       Z{ii} = [W1'*C*W1]./[W1'*eye(size(W1,1))*W1];
        maxLF_LS(:,ii) = max(abs(Outs.lead_fields'));
        waitbar(ii/length(BF_pos),h)
    end
    delete(h)
    % Frobenius norm of lead fields
    for n = 1:size(maxLF_LS,2)
        fro_norm_scalp_LS(n) = norm(maxLF_LS(:,n),'fro');
    end
    figure
    trisurf(meshes(1).tri,meshes(1).pnt(:,1),meshes(1).pnt(:,2),...
        meshes(1).pnt(:,3),fro_norm_scalp_LS*1e15)
    shading interp
    axis equal
    cb = colorbar;caxis([0 200]);cb.Label.String = 'Frobenius norm (fT)';
    hold on
    ft_plot_mesh(meshes(3),'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none')
    view(180,0)
    fig = gcf;
    fig.Color = [1,1,1];    
    % Save T stat
    files.results = ['T_',num2str(hp),'-',num2str(lp),'_dual_LS_VS'];
    if ~exist([path.results,files.results,'.mat'])
        save([path.results,files.results,'.mat'],'mu','condnr','T1')
    else
        save([path.results,files.results,'.mat'],'mu','condnr','T1','-append')
    end
    
elseif strcmpi(fwd_method,'BEM')
    % Create BEM model
    addpath(genpath([p4 '\MatLab files\BEM']))
    % Boundary meshes are described as a cell array Nx1, where each cell is a
    % struct that contains fields "p" for points (vertices), and "e" for
    % elements (faces, triangles). The triangles need to have CCW orientation
    % (when looking at the triangle from outside, the vertex indexing grows in
    % counterclockwise direction).
    
    % The convention of this BEM framework is that the innermost mesh has number
    % 1 and outermost mesh M; this is mandatory.
    
    % I recommend at least 2500 vertices for the innermost mesh; the outer
    % meshes can be a bit more coarse but not much (say, stay over 1500).
%     mri = ft_read_mri([path.mri files.mri]);
    if ~exist([path.mri 'BEM_meshes.mat'])
        cfg           = [];
        cfg.output    = {'brain','skull','scalp'};
        segmentedmri  = ft_volumesegment(cfg, mri);
        
        cfg             = [];
        cfg.tissue      = {'brain', 'skull', 'scalp'};
        cfg.numvertices = [2500, 2500, 2500];
        BEM_mesh            = ft_prepare_mesh(cfg, segmentedmri);
        save([path.mri '\BEM_meshes.mat'],'segmentedmri','BEM_mesh')
    else
        load([path.mri 'BEM_meshes.mat'])
    end
    
    % Define boundary meshes indexing from inside to outside --- scalp surface
    % must be the last one!
    innerskull.p = BEM_mesh(1).pos;
    innerskull.e = BEM_mesh(1).tri;
    outerskull.p = BEM_mesh(2).pos;
    outerskull.e = BEM_mesh(2).tri;
    scalp.p = BEM_mesh(3).pos;
    scalp.e = BEM_mesh(3).tri;
    bmeshes = {innerskull,outerskull,scalp};
    
    % These meshes are nested; you could check this to make sure that meshes are
    % in correct order
    bmeshes = hbf_SortNestedMeshes(bmeshes);
    
    % Set conductivities
    ci = [1 1/50 1]*.33; %conductivity inside each surface --- remember the order!
    
    % Setup sensor locs
    coils.p = sens_info.pos_used;
    coils.n = sens_info.ors_used;
    
    % Create BEM model
    if ~ismember('TBvol',who('-file',[path.skanect files.skanect]))
        TBvol = create_BEM(bmeshes,coils,ci);
        save([path.skanect files.skanect],'TBvol','-append')
    else
        load([path.skanect files.skanect],'TBvol')
        disp('Loaded Volume')
    end
    
    % Create lead fields (takes a fair while...)
    if ~ismember('LF_BEM',who('-file',[path.skanect files.skanect]))
        LF_BEM = zeros(size(coils.p,1),90,size(sourcepos,1));
        BEM_idx = zeros(size(sourcepos,1),1);
        LV_BEM = zeros(size(coils.p,1),size(sourcepos,1));
        for n = 1:size(sourcepos,1)
            disp((n./size(sourcepos,1))*100)
            dip_loc = sourcepos(n,:);
            [LF_BEM(:,:,n),UnitMDIP(:,:,n)] = Forward_OPM_BEM(dip_loc,bmeshes,coils,TBvol);
            clc
        end
        disp('Saving...')
        save([path.skanect files.skanect],'LF_BEM','UnitMDIP','-append')
        disp('Done!')
    else
        disp('Loading lead fields...')
        load([path.skanect files.skanect],'LF_BEM','UnitMDIP')
        disp('Done!')
    end
    
    % Beamforming
    % Compute beamformer weights
    for n = 1:size(sourcepos,1)
        clc
        lead_fields = LF_BEM(:,:,n);
        for angle = 1:size(lead_fields,2)
            Aa(:,angle) = (Cinv*lead_fields(:,angle))/(lead_fields(:,angle)'*...
                Cinv*lead_fields(:,angle));
            Zang_a(angle) = (Aa(:,angle)'*C*Aa(:,angle))/(Aa(:,angle)'*Noise_Cr*Aa(:,angle));
            X(angle) = (Aa(:,angle)'*C*Aa(:,angle))/(Aa(:,angle)'*(1e-30*Noise_Cr)*Aa(:,angle));
        end
        I = find(Zang_a==max(Zang_a));
        W1(:,n) = Aa(:,I);
        leads_vec(:,n) = lead_fields(:,I);
        Best_OR = UnitMDIP(I,:,n);
        best_ang(n) = (180/size(lead_fields,2))*I;
        Qa = W1(:,n)'*Ca*W1(:,n);
        Qc = W1(:,n)'*Cc*W1(:,n);
        T1{n} = (Qa-Qc)./(2*Qc);
        Z{n} = [W1(:,n)'*C*W1(:,n)]./[W1(:,n)'*eye(size(W1(:,n),1))*W1(:,n)];
        fprintf('Done %u/%u \n',n,length(sourcepos))
    end
    
    % Save stuff    
    % Save T stat
    files.results = ['T_',num2str(hp),'-',num2str(lp),'_dual_BEM_VS'];
    if ~exist([path.results,files.results,'.mat'])
        save([path.results,files.results,'.mat'],'mu','condnr','T1')
    else
        save([path.results,files.results,'.mat'],'mu','condnr','T1','-append')
    end
elseif strcmpi(fwd_method,'shell')
    addpath('C:\Users\ppxrh1\The University of Nottingham\OPM-MEG - Documents\MatLab files\BEM\')
    sensor_info.pos = sens_info.pos_used;
    sensor_info.ors = sens_info.ors_used;
    label = ([reshape([sens_info.sensornamesinuse; sens_info.sensornamesinuse],Nchans,1) ...
        OPM_vars.QZFM_data.axis]);
    for n = 1:size(label,1)
        sensor_info.label{n} = strjoin(label(n,:));
    end
    [LF_single_shell,L_reshaped] = Forward_OPM_shell_FT([path.mri,files.mri],sensor_info,sourcepos);
    inside_idx = find(LF_single_shell.inside);
    origin = mean(sourcepos);
    for n = 1:size(sourcepos,1)
        if LF_single_shell.inside(n)
            nn = find(inside_idx == n);
            % project to two tangential fields
            R = sourcepos(n,:) - origin;
            [etheta,ephi] = calctangent(R);
            % tan_ors = null(nrst_norm)';
            tan_ors = [etheta;ephi];
            ltan = L_reshaped(:,:,nn)*tan_ors';
            % get optimal combination for max SNR
            [v,d] = svd(ltan'*Cinv*ltan);
            [~,id] = min(diag(d));
            lopt = ltan*v(:,id); % turn to nAm amplitude
            W1(:,n) = (Cinv*lopt/(lopt'*Cinv*lopt));
            Qa = W1(:,n)'*Ca*W1(:,n);
            Qc = W1(:,n)'*Cc*W1(:,n);
            T1{n} = (Qa-Qc)./(2*Qc);
            Z{n} = [W1(:,n)'*C*W1(:,n)]./[W1(:,n)'*eye(size(W1(:,n),1))*W1(:,n)];

        else
            W1(:,n) = zeros(size(sensor_info.pos,1),1);
            T1{n} = 0;
            Z{n} = 0;
        end
        fprintf('Done %u/%u \n',n,length(sourcepos))
    end
    
    % Save stuff    
    % Save T stat
    files.results = ['T_',num2str(hp),'-',num2str(lp),'_dual_Shell_VS'];
    if ~exist([path.results,files.results,'.mat'])
        save([path.results,files.results,'.mat'],'mu','condnr','T1')
    else
        save([path.results,files.results,'.mat'],'mu','condnr','T1','-append')
    end
end
%% Remember the upside down anatomical problem, fieldtrip can fix that too.
if ~exist([path.mri 'acpc_space\'])
    mkdir([path.mri 'acpc_space\'])
end
% CTF space is a total nightmare with MRI visualisation software, so we
% need to put it in a normal cordinate space, such as ACPC.
mri1 = ft_convert_coordsys(mri, 'acpc');
% Do the final flipping, and then write out hi-res anatomical.
mri2 = align_ijk2xyz(mri1);
if ~exist([path.mri 'acpc_space\anat.nii'])
    ft_write_mri([path.mri 'acpc_space\anat.nii'],mri2.anatomy,'dataformat',...
        'nifti','transform',mri2.transform);
end
% for the 4mm t-stat image, just doing the same thing as above never seems
% to work, so we need to downsample the ACPC space to 4mm and use all the
% transformation matricies to take the the 4mm CTF image into 4mm ACPC
cfg = [];
cfg.downsample = 4;
mri4 = ft_volumedownsample(cfg,mri1);
% This insane multiplication is required to get the 4mm ctf-space image into
% the 4mm-acpc space
kernel=mri4.transform*inv(mri1.transform)*inv(mri1.head2headOrig)*inv(brainDS.transform*inv(brain.transform));

image = brainDS;
image.anatomy(voxID) = cell2mat(T1);
image.transform=kernel*image.transform;
% Do the final flipping, and then write out lo-res functional image.
image2 = align_ijk2xyz(image);
ft_write_mri([path.results,files.results,'.nii'],image2.anatomy,'dataformat','nifti',...
    'transform',image2.transform);
image = brainDS;
image.anatomy(voxID) = cell2mat(T1);
image.transform=kernel*image.transform;
% lets also output the downsample brain, just incase we need to use FSL
% later.
brainDS.transform=kernel*brainDS.transform;
% Do the final flipping, and then write out lo-res functional image.
image2 = align_ijk2xyz(brainDS);
if ~exist([path.mri 'acpc_space\brain_downsampled.nii'])
    ft_write_mri([path.mri 'acpc_space\brain_downsampled.nii'],image2.anatomy,...
        'dataformat','nifti','transform',image2.transform);
end
disp('Done')

%% Virtual electrode
cd(path.results)
filename = uigetfile('*.mat');
load(filename)
files.results = filename(1:end-4);
Tstat = ft_read_mri([path.results files.results,'.nii']);

fg = figure;
fg.Name = files.results;
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,cell2mat(T1),'filled')
colormap hot;caxis([min(cell2mat(T1)), -.2]);cb = colorbar;cb.Label.String = 'Tstat';
axis equal
view([120 20])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;

% Find voxel with min/max value of Tstat
[v l] = min(cell2mat(T1))
dip_loc = sourcepos(l,:)
fg = figure;
views = [-90 0;0 0;-90 90];
for sp = 1:3
    subplot(2,2,sp)
    ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
    hold on
    plot3(dip_loc(1),dip_loc(2),dip_loc(3),'.','markersize',20)
    view(views(sp,:))
end
fig = gcf;
fig.Color = 'w';
% In voxel to view in FSL
dip_vox = find(Tstat.anatomy == v);
[x,y,z] = ind2sub(size(Tstat.anatomy),dip_vox);
fsl_vox = [x y z]-1

% Or input voxel number from FSL view
% fsl_vox = [22 23 49];
% dip_loc = ft_warp_apply(inv(kernel)*Tstat.transform,fsl_vox+1);
% dip_loc = dip_loc/1000 % convert to metres

% Run beamformer for that location
C = cov(OPM_data_f_T');
maxEV = max(svd(C));
Noise_Cr = min(svd(C)).*eye(size(C));
C = C + mu.*maxEV.*eye(size(C));
Cinv = inv(C);

if strcmpi(fwd_method,'single')
    W1 = Quick_beamformer_target_OPM(C,Cinv,Noise_Cr,...
        dip_loc(1)*100,dip_loc(2)*100,dip_loc(3)*100,...
        sens_info.pos_used*100,sens_info.ors_used);
    
elseif strcmpi(fwd_method,'local')
    load([path.skanect,'local_spheres_info_',files.data(1:end-5),'.mat']);
    addpath('R:\Matlab_files\tools\');
    % duplicate the sphere origins for Dual mode
    for n = 1:Nlocs
        Loc_ind = OPM_data_struct.(['Sensor_',loc_names{n}]);
        this_axes = QZFM_data.axis(Loc_ind);
        for nn = 1:length(Loc_ind)
            origins(Loc_ind(nn),:) = so(n,:);
        end
    end
    vn = patchnormals(brain_mesh);
    source2mesh_dist = pdist2(dip_loc,brain_mesh.vertices);
    [~,s2m_idx] = min(source2mesh_dist);
    nrst_norm = vn(s2m_idx,:); % Nearest surface normal to the dipole location
    S.dip_loc = dip_loc;
    S.pos = sens_info.pos_used;
    S.ors = sens_info.ors_used;
    S.ls_origins = origins;
    S.tan_ref = nrst_norm;
    S.tan_method = 'norm';
    S.C = C;
    S.Noise_Cr = Noise_Cr;
    S.Cinv = Cinv;
    addpath('R:\Matlab_files\Local_Spheres\');
    [Outs_VE] = BF_local_spheres_target(S);
    % Get the stat
    W1 = Outs_VE.Weights;
    
elseif strcmpi(fwd_method,'BEM')    
    [LF_BEM_pk] = Forward_OPM_BEM(dip_loc,bmeshes,coils,TBvol);
    for angle = 1:size(LF_BEM_pk,2)
        Aa(:,angle) = (Cinv*LF_BEM_pk(:,angle))/(LF_BEM_pk(:,angle)'*...
            Cinv*LF_BEM_pk(:,angle));
        Zang_a(angle) = (Aa(:,angle)'*C*Aa(:,angle))/(Aa(:,angle)'*Noise_Cr*Aa(:,angle));
        X(angle) = (Aa(:,angle)'*C*Aa(:,angle))/(Aa(:,angle)'*(1e-30*Noise_Cr)*Aa(:,angle));
    end
    I = find(Zang_a==max(Zang_a));
    W1 = Aa(:,I);  
    
elseif strcmpi(fwd_method,'Shell')
    [~,L_reshaped] = Forward_OPM_shell_FT([path.mri,files.mri],sensor_info,dip_loc);
    
    % project to two tangential fields
    R = dip_loc - origin;
    [etheta,ephi] = calctangent(R);
    tan_ors = [etheta;ephi];
    ltan = L_reshaped(:,:,1)*tan_ors';
    
    % get optimal combination for max SNR
    [v,d] = svd(ltan'*Cinv*ltan);
    [~,id] = min(diag(d));
    lopt = ltan*v(:,id); % turn to nAm amplitude
    W1 = (Cinv*lopt/(lopt'*Cinv*lopt));
    
end

% Virtual electrode
VE = W1'*OPM_data_T;
VE_f = W1'*OPM_data_f_T;
% Hilbert envelope
H_VE_f = mean(reshape(abs(hilbert(VE_f)),duration*f,Ntrials),2);
H_VE_f = (H_VE_f - mean(H_VE_f(conwin(1):conwin(2))))./mean(H_VE_f(conwin(1):conwin(2)))*100;
% Plot
fg = figure;
fg.Name = 'Oscillatory amplitude';
plot(trial_time,H_VE_f)
grid on
xlabel('Time (s)');ylabel('Relative change (%)')
title(sprintf('Virtual electrode at location [%.2f,%.2f,%.2f] cm',dip_loc*100))
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;

% TFS of virtual electrode
VE_tfs = W1'*OPM_data1'*1e-15;
VE_tfs_fb = zeros(length(VE_tfs),length(fre));
for fb = 1:length(highpass)
    fprintf('\n Band %u/%u ',fb,length(highpass))
    filt_VE = nut_filter3(VE_tfs','butter','bp',3,highpass(fb),lowpass(fb),f,1)';
    VE_tfs_fb(:,fb) = abs(hilbert(filt_VE));
end
VE_tfs_mean = zeros(duration*f,length(fre));
for fb = 1:length(highpass)
    VE_tfs_fb_trials = [];
    % Chop data
    for i = 1:Ntrials
        VE_tfs_fb_trials = cat(1,VE_tfs_fb_trials,VE_tfs_fb(ind(i)+1:ind(i)+(duration*f),fb));
    end
    VE_tfs_filt = reshape(VE_tfs_fb_trials,duration*f,Ntrials);
    VE_tfs_mean(:,fb) = mean(VE_tfs_filt,2);
end
meanrest = mean(VE_tfs_mean(conwin(1):conwin(2),:),1);
meanrestmat = repmat(meanrest,size(VE_tfs_mean,1),1);
VE_TFS = (VE_tfs_mean'-meanrestmat')./meanrestmat'*100;
% plot
fg = figure;
fg.Name = 'TFS of virtual electrode';
pcolor(trial_time,fre,VE_TFS);shading interp
axis square;
title(sprintf('Virtual electrode at location [%.2f,%.2f,%.2f] cm',dip_loc*100))
xlabel('Time (s)');ylabel('Frequency (Hz)')
cb = colorbar; cb.Label.String = 'Percentage change (%)';
caxis([-80 80])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;
