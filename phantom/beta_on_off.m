
clear all

fname = 'R:\DRS-KidsOPM\phantom\20230309_155929_cMEG_Data\20230309_155929_meg.cMEG';
% load data
Data = read_cMEG_data_newAcq(fname);
% load helmet config file
helmet_config = readtable([fname(1:end-8),'HelmConfig.tsv'],'FileType','text','Delimiter','\t');
% remove useless line at end
helmet_config(end,:)=[];
%%
% remove tabs
for n_i = 1:length(Data.Chan_names)
    Data.Chan_names{n_i}(find(Data.Chan_names{n_i}==sprintf('\t')))=[];
end

% find time, sampling frequency and gain factor
f = Data.samp_frequency;
time = Data.time;
gain_line = split(Data.Session_info{strmatch('OPM Gain:',Data.Session_info)},': ');
gain_conv = 2.7.*str2num(gain_line{2}(1:end-1)); % V/nT

% get triggers
disp('Sorting Triggers')
idx = [];
for n = 1:size(Data.Chan_names,1)
    if strfind(Data.Chan_names{n},'Trigger')
        idx = [idx;n];
    end
end
triggers = Data.data(idx,:);

%get opm data channels
disp('Sorting OPM data')
idx = [];
for n = 1:size(Data.Chan_names,1)
    if strfind(Data.Chan_names{n},'Trigger')
        idx = [idx;n];
    elseif strfind(Data.Chan_names{n},'BNC')
        idx = [idx;n];
    end
end
OPM_names = Data.Chan_names;
OPM_names(idx,:) = [];
for n = 1:size(OPM_names,1)
    tmp = strsplit(OPM_names{n});
    Names_only(n,1) = tmp(1);
    Axis_only(n,1) = tmp(2);
end
OPM_data1 = Data.data;
OPM_data1(idx,:) = [];

% Convert to fT
OPM_data1 = (1e6.*OPM_data1)./gain_conv;

% Mean correct
OPM_data1 = OPM_data1 - mean(OPM_data1,2);

% OPM data in layout order
load([fname(1:end-5) '_sensor_order.mat'])
loc_info = cell(size(T.Name));
OPM_data = [];
sensornames = [];
count = 0;
for n = 1:size(T,1)
    idx = [];
    if ~isempty(T.Name{n})
        sens_name = strsplit(T.Name{n});
        idx = find(strcmpi(sens_name(1),Names_only));
        for m = 1:length(idx)
            count = count + 1;
            OPM_data_mat.(['Sensor_' sens_name{1}]).([Axis_only{idx(m)}(2)]) = OPM_data1(idx(m),:);
            loc_info{n} = ['Sensor ' sens_name{1}];
            sensornames{count} = [sens_name{1} ' ' Axis_only{idx(m)}(2)];
            OPM_data = [OPM_data;OPM_data1(idx(m),:)];
        end
    end
end
clear OPM_data1
sensornames = sensornames';
% 
% % Helmet positions and orientations
for n = 1:size(sensornames,1)
    tmp = split(sensornames{n});
    loc_names{n} = tmp{1};
    loc_axis{n} = tmp{2};
%     loc_idx(n) = unique(find(strcmpi(['Sensor ' loc_names{n}],loc_info)));
    ch_index = find(strcmp(sprintf('%s %s',loc_names{n},loc_axis{n}),helmet_config.Sensor))
    Sens_pos(n,:) = [helmet_config.Px(ch_index),helmet_config.Py(ch_index),helmet_config.Pz(ch_index)];
    Sens_ors(n,:) =[helmet_config.Ox(ch_index),helmet_config.Oy(ch_index),helmet_config.Oz(ch_index)];;
end
figure
scatter3(Sens_pos(:,1),Sens_pos(:,2),Sens_pos(:,3))
quiver3(Sens_pos(:,1),Sens_pos(:,2),Sens_pos(:,3),Sens_ors(:,1),Sens_ors(:,2),Sens_ors(:,3))


% get rid of pre and post task
not_during_experiment = (time < 5) | (time > 410);
time(not_during_experiment) = [];
triggers(:,not_during_experiment) = [];
OPM_data(:,not_during_experiment) = [];
time = time - time(1);

figure

plot(time,triggers(1,:));
N_trials = 100;
% xlim([0,10])
hold on
diff_trig = [0,diff(triggers(1,:))];
plot(time,diff_trig)
% find beginnning
beginning_samples = find(diff_trig>0.5,N_trials);
end_samples = find(diff_trig<-0.5,N_trials);
plot(time(beginning_samples),0.5.*ones(size(beginning_samples)),'ro')
plot(time(end_samples),-0.5.*ones(size(end_samples)),'bo')
trl_duration_samps = floor(mean(end_samples - beginning_samples));
%%
%% Make FT data structure
disp('Making FT structure')
if ~exist([fname(1:end-8),'ch_table.mat'])
    ch_table = table();
    ch_table.name = sensornames;
    sens_info.pos = Sens_pos;
    sens_info.ors = Sens_ors;
    
    for ch_i = 1:height(ch_table)
        ch_table_(ch_i,:) = cat(2, ch_table(ch_i,1), helmet_config(startsWith(helmet_config.Sensor, ch_table.name{ch_i}),2:end-1));
    end
    ch_table = ch_table_;
    ch_table.status = repmat({'good'},height(ch_table),1);
    
    for sl_i = 1:height(T)
        if ~isempty(T{sl_i,1}{1})
            ch_table.slot_no(startsWith(ch_table.name,T{sl_i,1}{1}(1:2))) = sl_i;
        end
    end
    
    get_good_channels(OPM_data,ch_table,f);
    [ch_table] = Bad_Channels(OPM_data',ch_table,f);
    save([fname(1:end-8),'ch_table.mat'],'ch_table')
else 
    load([fname(1:end-8),'_triax_only.mat'],'ch_table')
end
%%
disp("Removing bad channels")
bad_chans_data = [find(startsWith(ch_table.status,'bad'))];
ch_table(bad_chans_data,:) = [];
OPM_data(bad_chans_data,:) = [];

%% sensor info
ch_table.Py = ch_table.Py-0.01;
ch_table.Pz = ch_table.Pz+0.01;
helmet_config.Py = helmet_config.Py - 0.01;
helmet_config.Pz = helmet_config.Pz + 0.01;   

S.sensor_info.pos = [ch_table.Px,ch_table.Py,ch_table.Pz];
S.sensor_info.ors = [ch_table.Ox,ch_table.Oy,ch_table.Oz];

%% Mean field correction
N = S.sensor_info.ors; % orientation matrix (N_sens x 3)
S.M = eye(length(N)) - N*pinv(N);

% source positions
addpath F:\Rdrive\movie\scripts\fieldtrip-20190212
ft_defaults
mri = ft_read_mri('MNI152_T1_1mm.nii');
mri = ft_convert_units(mri,'m');

%%
 
cfg = [];
cfg.downsample = 4;
[downsample] = ft_volumedownsample(cfg, mri);

 
 %%
[sourcepos_vox(:,1),sourcepos_vox(:,2),sourcepos_vox(:,3)] = ind2sub(downsample.dim, find(ones(size(downsample.anatomy))));
sourcepos = ft_warp_apply(downsample.transform,sourcepos_vox);
S.mri_file = 'R:\DRS-KidsOPM\phantom\MNI152_T1_1mm.nii';
addpath F:\Rdrive\movie\scripts\Beamformer
[bf_outs_shell] = Copy_of_run_beamformer('shell',sourcepos,S,0,[],1);
%%
lead_fields_shell_xyz = bf_outs_shell.LF;

load meshes.mat
X = meshes.pnt; Origin = mean(X,1);
Ndips = length(sourcepos);
outside=[];
for n = 1:Ndips
    thispos = sourcepos(n,:)';
    [phi,theta1,r] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,thispos(3) - Origin(3));
    theta = pi/2 - theta1;
    
    Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
    Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
    Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
    Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);

    if all(Lead_fields(:,:,n) == zeros(size(Lead_fields(:,:,n))),'all')
        outside(n)=true;
    else
        outside(n)=false;
    end
end
Lead_fields(:,:,outside==1) = [];
sourcepos(outside==1,:) = [];
Src_Or_theta(outside==1,:) = [];
Src_Or_phi(outside==1,:) = [];
Ndips = length(sourcepos);
%%  filter data
OPM_data_f = OPM_data;
for harms = [50,100,150,200]
    Wo = harms/(f/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    disp(['Applying Notch filter:',num2str(harms)])
    OPM_data_f = filter(b,a,OPM_data_f',[],1)';
end

hp = 13;
lp = 30;
[b,a] = butter(4,2*[hp lp]/f);
OPM_data_f = [filtfilt(b,a,OPM_data_f')]';

% apply mean field correction
disp("Applying mean field correction")
OPM_data_mfc = S.M*OPM_data_f;

%% chop data
data_f_mat = zeros(size(OPM_data_mfc,1),trl_duration_samps*2,size(beginning_samples,1));
Ca = zeros(size(OPM_data_mfc,1),size(OPM_data_mfc,1),size(beginning_samples,1));
Cc = Ca;
for tr_i = 1:size(beginning_samples,1)
    data_f_mat(:,:,tr_i) = OPM_data_mfc(:,beginning_samples(tr_i):beginning_samples(tr_i)+trl_duration_samps*2-1);

    Ca(:,:,tr_i) = cov(data_f_mat(:,1:trl_duration_samps-1,tr_i)');
    Cc(:,:,tr_i) = cov(data_f_mat(:,trl_duration_samps:trl_duration_samps*2 -1,tr_i)');
end

C =mean((Ca + Cc)./2,3);
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);

 for n = 1:Ndips
        this_L = Lead_fields(:,:,n);

        iPower_v = this_L'*Cr_inv*this_L;
        [v,d] = svd(iPower_v);
        [~,id] = min(diag(d));
        lopt = this_L*v(:,id); % turn to nAm amplitude
        w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
        T_stat(n) = ((w*Ca*w') - (w*Cc*w'))./(2*(w*Cc*w'));
 end
% find max and calculate VE
n = find(T_stat==max(T_stat));
peak_pos = sourcepos(n,:);
peak_or_theta = Src_Or_theta(n,:)./100;
peak_or_phi = Src_Or_phi(n,:)./100;
this_L = Lead_fields(:,:,n);

iPower_v = this_L'*Cr_inv*this_L;
[v,d] = svd(iPower_v);
[~,id] = min(diag(d));
lopt = this_L*v(:,id); % turn to nAm amplitude
w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
control_inds = trl_duration_samps:trl_duration_samps*2 -1;
VE = (w*OPM_data_mfc)./sqrt(w*w');
for tr_i = 1:size(beginning_samples,1)
    VE_mat(:,tr_i) = VE(beginning_samples(tr_i):beginning_samples(tr_i)+trl_duration_samps*2-1);
    Envs(:,tr_i) = abs(hilbert(VE_mat(:,tr_i)));
    Envs(:,tr_i) = (Envs(:,tr_i) - mean(Envs(control_inds,tr_i)))./mean(Envs(control_inds,tr_i));
end
mean_Env = mean(Envs,2);
figure
trl_time = linspace(1,trl_duration_samps*2./f,trl_duration_samps*2);
plot(trl_time,mean_Env)
figure
plot(time,VE./max(VE))
hold on
plot(time,triggers)

%% plot T stat map
fg0 = figure;

subplot(2,3,[1 2 4 5])
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,T_stat,'filled')
colormap hot;
% caxis([cax_val]);
% % cb = colorbar;cb.Label.String = 'Tstat';
axis equal
view([120 20])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;

plot3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),'ko')
phantom_slot_z = startsWith(helmet_config.Sensor,'LE Z');
phantom_slot_x= startsWith(helmet_config.Sensor,'LE X');
phantom_tru_ori = [helmet_config.Ox(phantom_slot_x),...
    helmet_config.Oy(phantom_slot_x),...
    helmet_config.Oz(phantom_slot_x)]./(100/1);

plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r')
phantom_ori = [helmet_config.Ox(phantom_slot_z),...
    helmet_config.Oy(phantom_slot_z),...
    helmet_config.Oz(phantom_slot_z)]./(100/4); % 4cm in radially from sensor

phantom_slot_pos = [helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),helmet_config.Pz(phantom_slot_z)];
arr = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_ori(1),phantom_ori(2),phantom_ori(3));
arr.LineWidth = 2;arr.Color = 'k';
dipole_pos = phantom_slot_pos + phantom_ori;
plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',8,'MarkerFaceColor','g')
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    phantom_tru_ori(1),phantom_tru_ori(2),phantom_tru_ori(3),1);
arr_true.LineWidth = 2;arr_true.Color = 'g';


plot3(peak_pos(1),peak_pos(2),peak_pos(3),'bo','MarkerSize',8,'MarkerFaceColor','b')
sc=2;
arr_th = quiver3(peak_pos(1),peak_pos(2),peak_pos(3),...
    peak_or_theta(1),peak_or_theta(2),peak_or_theta(3),sc);
arr_th.LineWidth = 2;arr_th.Color = 'c';

arr_ph = quiver3(peak_pos(1),peak_pos(2),peak_pos(3),...
    peak_or_phi(1),peak_or_phi(2),peak_or_phi(3),sc);
arr_ph.LineWidth = 2;arr_ph.Color = 'm';

opt_ori = -(peak_or_theta.*v(1,id)+peak_or_phi.*v(2,id));
arr_opt = quiver3(peak_pos(1),peak_pos(2),peak_pos(3),...
    opt_ori(1),opt_ori(2),opt_ori(3),sc);
arr_opt.LineWidth = 3;arr_opt.Color = 'r';

drawnow
pause(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now again in ROI and high res
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mri = ft_read_mri('MNI152_T1_1mm.nii');
mri = ft_convert_units(mri,'m');

%%
cfg = [];
cfg.downsample = 1;
[downsample] = ft_volumedownsample(cfg, mri);

 
 %%
clear sourcepos_vox
[sourcepos_vox(:,1),sourcepos_vox(:,2),sourcepos_vox(:,3)] = ind2sub(downsample.dim, find(ones(size(downsample.anatomy))));
sourcepos = ft_warp_apply(downsample.transform,sourcepos_vox);
within_sphere = (sqrt(sum((sourcepos - dipole_pos).^2,2)) < 35/1000);
sourcepos(~within_sphere,:) = [];
fprintf('Beamforming %d positions\n',size(sourcepos,1));
S.mri_file = 'R:\DRS-KidsOPM\phantom\MNI152_T1_1mm.nii';
addpath F:\Rdrive\movie\scripts\Beamformer
[bf_outs_shell] = Copy_of_run_beamformer('shell',sourcepos,S,0,[],1);
%%
lead_fields_shell_xyz = bf_outs_shell.LF;

load meshes.mat
X = meshes.pnt; Origin = mean(X,1);
Ndips = length(sourcepos);
clear outside Lead_fields Src_Or_phi Src_Or_theta T_stat
Ndips = length(sourcepos);
outside=[];
for n = 1:Ndips
    thispos = sourcepos(n,:)';
    [phi,theta1,r] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,thispos(3) - Origin(3));
    theta = pi/2 - theta1;
    
    Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
    Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
    Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
    Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);

    if all(Lead_fields(:,:,n) == zeros(size(Lead_fields(:,:,n))),'all')
        outside(n)=true;
    else
        outside(n)=false;
    end
end
Lead_fields(:,:,outside==1) = [];
sourcepos(outside==1,:) = [];
Src_Or_theta(outside==1,:) = [];
Src_Or_phi(outside==1,:) = [];
Ndips = length(sourcepos);
for n = 1:Ndips
        this_L = Lead_fields(:,:,n);

        iPower_v = this_L'*Cr_inv*this_L;
        [v,d] = svd(iPower_v);
        [~,id] = min(diag(d));
        lopt = this_L*v(:,id); % turn to nAm amplitude
        w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
        T_stat(n) = ((w*Ca*w') - (w*Cc*w'))./(2*(w*Cc*w'));
 end
% find max and calculate VE
n = find(T_stat==max(T_stat));
peak_pos = sourcepos(n,:);
peak_or_theta = Src_Or_theta(n,:)./100;
peak_or_phi = Src_Or_phi(n,:)./100;
this_L = Lead_fields(:,:,n);

iPower_v = this_L'*Cr_inv*this_L;
[v,d] = svd(iPower_v);
[~,id] = min(diag(d));
lopt = this_L*v(:,id); % turn to nAm amplitude
w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
control_inds = trl_duration_samps:trl_duration_samps*2 -1;
VE = (w*OPM_data_mfc)./sqrt(w*w');
for tr_i = 1:size(beginning_samples,1)
    VE_mat(:,tr_i) = VE(beginning_samples(tr_i):beginning_samples(tr_i)+trl_duration_samps*2-1);
    Envs(:,tr_i) = abs(hilbert(VE_mat(:,tr_i)));
    Envs(:,tr_i) = (Envs(:,tr_i) - mean(Envs(control_inds,tr_i)))./mean(Envs(control_inds,tr_i));
end
mean_Env = mean(Envs,2);
figure
trl_time = linspace(1,trl_duration_samps*2./f,trl_duration_samps*2);
plot(trl_time,mean_Env)
figure
plot(time,VE./max(VE))
hold on
plot(time,triggers)


%% plot z stat map
fg = figure;

subplot(2,4,[1 2 3 5 6 7])
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
s_plt = scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,T_stat,'filled');

colormap hot;
% caxis([cax_val]);
% cb = colorbar;cb.Label.String = 'Tstat';
axis equal
view([120 20])
fig = gcf;
fig.Color = 'w';
ax = gca;
ax.FontSize = 14;

plot3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),'ko')
phantom_slot_z = startsWith(helmet_config.Sensor,'LE Z');
phantom_slot_x= startsWith(helmet_config.Sensor,'LE X');
phantom_tru_ori = [helmet_config.Ox(phantom_slot_x),...
    helmet_config.Oy(phantom_slot_x),...
    helmet_config.Oz(phantom_slot_x)]./(100/1);

plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r',...
    'DisplayName','Phantom Slot')
phantom_ori = [helmet_config.Ox(phantom_slot_z),...
    helmet_config.Oy(phantom_slot_z),...
    helmet_config.Oz(phantom_slot_z)]./(100/4); % 4cm in radially from sensor

phantom_slot_pos = [helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),helmet_config.Pz(phantom_slot_z)];
arr = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_ori(1),phantom_ori(2),phantom_ori(3),'DisplayName','Phantom axis');
arr.LineWidth = 2;arr.Color = 'k';
dipole_pos = phantom_slot_pos + phantom_ori;
plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',8,...
    'MarkerFaceColor','g','DisplayName','Ground truth position')
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    phantom_tru_ori(1),phantom_tru_ori(2),phantom_tru_ori(3),3,...
    'DisplayName','Ground truth direction');
arr_true.LineWidth = 2;arr_true.Color = 'g';

fprintf('Peak to truth distance = %1.3f mm\n',sqrt(sum((dipole_pos-peak_pos).^2))*1000);
plot3(peak_pos(1),peak_pos(2),peak_pos(3),'bo','MarkerSize',8,'MarkerFaceColor','b',...
    'DisplayName','Beamformer Peak')
sc=2;
arr_th = quiver3(peak_pos(1),peak_pos(2),peak_pos(3),...
    peak_or_theta(1),peak_or_theta(2),peak_or_theta(3),sc,...
    'DisplayName','L.F. \theta');
arr_th.LineWidth = 2;arr_th.Color = 'c';

arr_ph = quiver3(peak_pos(1),peak_pos(2),peak_pos(3),...
    peak_or_phi(1),peak_or_phi(2),peak_or_phi(3),sc,'DisplayName','L.F. \phi');
arr_ph.LineWidth = 2;arr_ph.Color = 'm';

opt_ori = -(peak_or_theta.*v(1,id)+peak_or_phi.*v(2,id));
arr_opt = quiver3(peak_pos(1),peak_pos(2),peak_pos(3),...
    opt_ori(1),opt_ori(2),opt_ori(3),sc,'DisplayName','B.F. opt orientation');
arr_opt.LineWidth = 3;arr_opt.Color = 'r';
ll = legend;ll.Position(1) = 0.6;
ax=gca;set(ax.Children(end-4:end),'HandleVisibility','off');
fg_copy = copyobj(fg,0);
s_plt.Visible=0;
beta_phantom_peak.pos = peak_pos;
beta_phantom_peak.ori = opt_ori;
save beta_phantom_peak beta_phantom_peak