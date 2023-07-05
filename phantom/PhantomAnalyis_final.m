clear all
close all
set(0,'DefaultFigureWindowStyle','normal')

% fname = 'F:\Rdrive\movie\scripts\phantom\20230309_164734_cMEG_Data\20230309_164734_meg.cMEG';phantom_sensor = 'LE';
fname = 'F:\Rdrive\movie\scripts\phantom\20230330\20230330_160330_cMEG_Data_brainnoise_2xSNR\20230330_160330_meg.cMEG';
addpath F:\Rdrive\movie\scripts\fieldtrip-20190212
ft_defaults

res_in_mm = 1; % downsample from 1mm to x mm resolution
use_3_dirs = 0; % use 3 direction beamformer instead of 2
phantom_sensor = 'KF';
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
hold on
quiver3(Sens_pos(:,1),Sens_pos(:,2),Sens_pos(:,3),Sens_ors(:,1),Sens_ors(:,2),Sens_ors(:,3))

figure

plot(time,[0,diff(triggers(1,:))]);
xlim([0,10])

% find beginnning
beginning_sample = find(triggers(1,:)>0.1,1);
% chop data
duration = 600 ;%s
OPM_data = OPM_data(:,beginning_sample+1:beginning_sample+duration*f);
triggers = triggers(1,beginning_sample+1:beginning_sample+duration*f);
time = time(beginning_sample+1:beginning_sample+duration*f);time=time-time(1);

figure

plot(time,triggers);
xlim([0,10])

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
    %remove all axes from a sensor even if only one channel is bad and save
    %as ..._triax_only.mat'
    save([fname(1:end-8),'ch_table.mat'],'ch_table')
else
    load([fname(1:end-8),'_triax_only.mat'],'ch_table')
    % %     save([fname(1:end-8),'_triax_only.mat'],'ch_table')
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

mri = ft_read_mri('MNI152_T1_1mm.nii');
mri = ft_convert_units(mri,'m');


%% define ground truth
phantom_slot_z = startsWith(helmet_config.Sensor,[phantom_sensor,' Z']);
phantom_slot_x= startsWith(helmet_config.Sensor,[phantom_sensor,' X']);
dipole_ori = -[helmet_config.Ox(phantom_slot_x),...
    helmet_config.Oy(phantom_slot_x),...
    helmet_config.Oz(phantom_slot_x)]./(100/1);
dipole_ori = (rotz(-20)*dipole_ori')';
phantom_slot_pos = [helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),helmet_config.Pz(phantom_slot_z)];
phantom_axis = [helmet_config.Ox(phantom_slot_z),...
    helmet_config.Oy(phantom_slot_z),...
    helmet_config.Oz(phantom_slot_z)]./(100/4); % 4cm in radially from sensor

dipole_pos = phantom_slot_pos + phantom_axis;

%%

cfg = [];
cfg.downsample = res_in_mm;
[downsample] = ft_volumedownsample(cfg, mri);

%%
[sourcepos_vox(:,1),sourcepos_vox(:,2),sourcepos_vox(:,3)] = ind2sub(downsample.dim, find(ones(size(downsample.anatomy))));
sourcepos = ft_warp_apply(downsample.transform,sourcepos_vox);
if cfg.downsample == 1
    within_sphere = (sqrt(sum((sourcepos - dipole_pos).^2,2)) < 35/1000);
    sourcepos(~within_sphere,:) = [];
end
S.mri_file = 'F:\Rdrive\movie\scripts\phantom\\MNI152_T1_1mm.nii';
addpath F:\Rdrive\movie\scripts\Beamformer
[bf_outs_shell] = Copy_of_run_beamformer('shell',sourcepos,S,0,[],1);
%%
lead_fields_shell_xyz = bf_outs_shell.LF;
%mri file meshes
load meshes.mat

X = meshes.pnt; Origin = mean(X,1);
Ndips = length(sourcepos);
outside=[];
for n = 1:Ndips
    thispos = sourcepos(n,:)';
    [phi,theta1,r] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,...
        thispos(3) - Origin(3));
    theta = pi/2 - theta1;

    Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
    Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
    Src_Or_R(n,:) = [thispos(1) - Origin(1),thispos(2) - Origin(2),...
        thispos(3) - Origin(3)];
    Src_Or_R(n,:) = Src_Or_R(n,:)./norm(Src_Or_R(n,:));
    %     Src_Or_R(n,:) = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)];

    Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
    Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);

    if use_3_dirs
        Lead_fields(:,3,n) = Src_Or_R(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_R(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_R(n,3)*lead_fields_shell_xyz(:,3,n);
    end
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
Src_Or_R(outside==1,:) = [];
lead_fields_shell_xyz(:,:,outside==1) = [];
Ndips = length(sourcepos);


%%
ff = figure;
load meshes.mat
ff.Color = 'w';
subplot(2,3,[1 2 4 5])
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on

view([200,30])

isz = endsWith(ch_table.name,'Z');
isx = endsWith(ch_table.name,'X');
isy = endsWith(ch_table.name,'Y');

plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r')

arr_phantom_ax = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_axis(1),phantom_axis(2),phantom_axis(3));
arr_phantom_ax.LineWidth = 2;arr_phantom_ax.Color = 'k';arr_phantom_ax.LineStyle=':';
h_dip_pos = plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',10,...
    'MarkerFaceColor','k');
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    dipole_ori(1),dipole_ori(2),dipole_ori(3),3);
arr_true.LineWidth = 3;arr_true.Color = 'k';
% take LF from actual loc
actual_location = 1;
if actual_location
    [bf_outs_dip] = Copy_of_run_beamformer('shell',dipole_pos,S,0,[],1);
    lead_fields_shell_dip_xyz = bf_outs_dip.LF;
else
    % Take LF from nearest voxel
    [~,nearest_i] = min(sum((sourcepos - dipole_pos).^2,2));
    lead_fields_shell_dip_xyz = lead_fields_shell_xyz(:,1,nearest_i);
end
true_LF = (dipole_ori(1)*lead_fields_shell_dip_xyz(:,1) + ...
    dipole_ori(2)*lead_fields_shell_dip_xyz(:,2) + ...
    dipole_ori(3)*lead_fields_shell_dip_xyz(:,3))./norm(dipole_ori);
%% dipole leadfields

[phi,theta1,r] = cart2sph(dipole_pos(1) - Origin(1),dipole_pos(2) - Origin(2),...
    dipole_pos(3) - Origin(3));
theta = pi/2 - theta1;

Src_Or_theta_dip = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
Src_Or_phi_dip = [-sin(phi) cos(phi) 0];
Src_Or_R_dip = [dipole_pos(1) - Origin(1),dipole_pos(2) - Origin(2),...
    dipole_pos(3) - Origin(3)];
Src_Or_R_dip = Src_Or_R_dip./norm(Src_Or_R_dip);

Lead_fields_dip(:,1) = Src_Or_theta_dip(1)*lead_fields_shell_dip_xyz(:,1) + ...
    Src_Or_theta_dip(2)*lead_fields_shell_dip_xyz(:,2) + ...
    Src_Or_theta_dip(3)*lead_fields_shell_dip_xyz(:,3);

Lead_fields_dip(:,2) = Src_Or_phi_dip(1)*lead_fields_shell_dip_xyz(:,1) + ...
    Src_Or_phi_dip(2)*lead_fields_shell_dip_xyz(:,2) + ...
    Src_Or_phi_dip(3)*lead_fields_shell_dip_xyz(:,3);

Lead_fields_dip(:,3) = Src_Or_R_dip(1)*lead_fields_shell_dip_xyz(:,1) + ...
    Src_Or_R_dip(2)*lead_fields_shell_dip_xyz(:,2) + ...
    Src_Or_R_dip(3)*lead_fields_shell_dip_xyz(:,3);

figure
set(gcf,'Name','GT Leadfields Z channels')
titles = ["\theta","\phi","R"];
comp_vec = {Src_Or_theta_dip,Src_Or_phi_dip,Src_Or_R_dip};
for s = 1:3
    subplot(1,3,s)
    ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
    hold on
    title(titles(s))
    view([200,30])


    plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
        helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r')

    arr_phantom_ax = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
        phantom_axis(1),phantom_axis(2),phantom_axis(3));
    arr_phantom_ax.LineWidth = 2;arr_phantom_ax.Color = 'k';arr_phantom_ax.LineStyle=':';
    h_dip_pos = plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',10,...
        'MarkerFaceColor','k');
    arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
        dipole_ori(1),dipole_ori(2),dipole_ori(3),3);

    arr_true.LineWidth = 3;arr_true.Color = 'k';

    field_mag_dir = Lead_fields_dip(isz,s);
    DT= boundary(S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),0);
    h=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),field_mag_dir);
    h.FaceColor = 'interp';h.EdgeColor = 'none';h.FaceAlpha = 0.5;colormap parula
    colorbar('Location','south')
    clim([-1,1].*0.5e-14);%colormap hot
    ori = comp_vec{s};comp_col = [0,0,0];comp_col(s) = .5;
    quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
        ori(1),ori(2),ori(3),.05,'Color',comp_col,'LineWidth',5);
end

%%

figure(ff)
field_mag = sqrt(true_LF(isx).^2 + true_LF(isy).^2 + true_LF(isz).^2);
field_mag = true_LF(isz);
DT= boundary(S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),0);
h=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),field_mag);
h.FaceColor = 'interp';h.EdgeColor = 'none';h.FaceAlpha = 0.5;colormap parula
clim([-1,1].*1e-14);%colormap hot

slots_in_use = unique(ch_table.slot_no);
for sl_i = 1:length(slots_in_use)
    x_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isx;
    y_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isy;
    z_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isz;
    xv = [ch_table.Ox(x_ch), ch_table.Oy(x_ch), ch_table.Oz(x_ch)];
    yv = [ch_table.Ox(y_ch), ch_table.Oy(y_ch), ch_table.Oz(y_ch)];
    zv = [ch_table.Ox(z_ch), ch_table.Oy(z_ch), ch_table.Oz(z_ch)];
    qx(sl_i) = quiver3(ch_table.Px(x_ch),ch_table.Py(x_ch),ch_table.Pz(x_ch),...
        xv(1),xv(2),xv(3),0.01,'r');
    qy(sl_i) = quiver3(ch_table.Px(x_ch),ch_table.Py(x_ch),ch_table.Pz(x_ch),...
        yv(1),yv(2),yv(3),0.01,'g');
    qz(sl_i) = quiver3(ch_table.Px(x_ch),ch_table.Py(x_ch),ch_table.Pz(x_ch),...
        zv(1),zv(2),zv(3),0.01,'b');
    LF_mags_true(sl_i,:) = true_LF(x_ch)*xv + true_LF(y_ch)*yv + true_LF(z_ch)*zv;

    arrs(sl_i) = quiver3(ch_table.Px(x_ch),ch_table.Py(x_ch),ch_table.Pz(x_ch),...
        LF_mags_true(sl_i,1),LF_mags_true(sl_i,2),LF_mags_true(sl_i,3),1e12,'r');
end
set(arrs,'LineWidth',2,'Color','b','AutoScaleFactor',2e12);
h.Visible = "on";set([qx,qy,qz],'Visible',0);

% % random error in pos and source angle
% N_sim = 500;
% sd_ang = 10; % degree
% sd_loc = 5/1000; % m
% rand_locs = dipole_pos +  sd_loc.*randn(N_sim,3);
% [dipole_ori_sph(1),dipole_ori_sph(2),dipole_ori_sph(3)] = cart2sph(dipole_ori(1),...
%     dipole_ori(2),dipole_ori(3));
% rand_oris_sph = dipole_ori_sph + [deg2rad(sd_ang).*randn(N_sim,2),zeros(N_sim,1)];
% clear rand_oris
% [rand_oris(:,1), rand_oris(:,2), rand_oris(:,3)] = sph2cart(rand_oris_sph(:,1),...
%     rand_oris_sph(:,2),rand_oris_sph(:,3));
%
% plot3(rand_locs(:,1),rand_locs(:,2),rand_locs(:,3),'kx','MarkerSize',2);
% quiver3(rand_locs(:,1),rand_locs(:,2),rand_locs(:,3),...
%     rand_oris(:,1),rand_oris(:,2),rand_oris(:,3),3,'k','LineWidth',1.5);
% h_dip_pos.MarkerSize =1;
%
% % simulate LFs from random locs and oris
% [bf_outs_rand] = Copy_of_run_beamformer('shell',rand_locs,S,0,[],1);
% lead_fields_shell_rand_xyz = bf_outs_rand.LF;
%
% fCount=71;
% figure(ff)
% f_gif = getframe(ff);
% [im,map] = rgb2ind(f_gif.cdata,256,'nodither');
% im(1,1,1,fCount) = 0;
% k = 1;
% for sim_i = 1:N_sim
%     for sl_i = 1:length(slots_in_use)
%         x_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isx;
%         y_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isy;
%         z_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isz;
%         xv = [ch_table.Ox(x_ch), ch_table.Oy(x_ch), ch_table.Oz(x_ch)];
%         yv = [ch_table.Ox(y_ch), ch_table.Oy(y_ch), ch_table.Oz(y_ch)];
%         zv = [ch_table.Ox(z_ch), ch_table.Oy(z_ch), ch_table.Oz(z_ch)];
%
%         rand_LF = (rand_oris(sim_i,1)*lead_fields_shell_rand_xyz(:,1,sim_i) + ...
%             rand_oris(sim_i,2)*lead_fields_shell_rand_xyz(:,2,sim_i) + ...
%             rand_oris(sim_i,3)*lead_fields_shell_rand_xyz(:,3,sim_i))./norm(rand_oris(sim_i,:));
%         LF_mags(sl_i,:) = rand_LF(x_ch)*xv + rand_LF(y_ch)*yv + rand_LF(z_ch)*zv;
%         px(sl_i) = ch_table.Px(x_ch);
%         py(sl_i) = ch_table.Py(x_ch);
%         pz(sl_i) = ch_table.Pz(x_ch);
%     end
%     arrs_unc(sl_i,sim_i) = quiver3(px',py',pz',...
%         LF_mags(:,1),LF_mags(:,2),LF_mags(:,3),1,'k');
%     sim_similarity(sim_i) = abs(pdist2(LF_mags(:)',LF_mags_true(:)','cosine') -1);
%
%     drawnow
%     f_gif = getframe(ff);
%     im(:,:,1,k) = rgb2ind(f_gif.cdata,map,'nodither');
%     k = k + 1;
% end
% imwrite(im,map,'Animation_rand_sim_2xsnr.gif','DelayTime',0.1,'LoopCount',inf)
%%

% filter data
OPM_data_f = OPM_data;
triggers_f = triggers;
for harms = [50,100,150,200]
    Wo = harms/(f/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    disp(['Applying Notch filter:',num2str(harms)])
    OPM_data_f = filter(b,a,OPM_data_f',[],1)';
    triggers_f = filter(b,a,triggers_f',[],1)';
end

hp = 1;
lp = 100;
[b,a] = butter(4,2*[hp lp]/f);
OPM_data_f = [filtfilt(b,a,OPM_data_f')]';
triggers_f = [filtfilt(b,a,triggers_f')]';

% apply mean field correction
disp("Applying mean field correction")
OPM_data_mfc = S.M*OPM_data_f;

%
%
for ti = 1:720000
    % waitforbuttonpress
    corrs(ti) = corr(true_LF,OPM_data_mfc(:,ti));
end

figure
subplot(2,1,1)
plot(corrs)
subplot(2,1,2)
histogram(abs(corrs))

%% %
for sl_i = 1:length(slots_in_use)
    x_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isx;
    y_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isy;
    z_ch = (ch_table.slot_no == slots_in_use(sl_i)) & isz;

    px(sl_i) = ch_table.Px(x_ch);
    py(sl_i) = ch_table.Py(x_ch);
    pz(sl_i) = ch_table.Pz(x_ch);
    xv(sl_i,:) = [ch_table.Ox(x_ch), ch_table.Oy(x_ch), ch_table.Oz(x_ch)];
    yv(sl_i,:) = [ch_table.Ox(y_ch), ch_table.Oy(y_ch), ch_table.Oz(y_ch)];
    zv(sl_i,:) = [ch_table.Ox(z_ch), ch_table.Oy(z_ch), ch_table.Oz(z_ch)];
end

figure(ff)
% ff=figure
view([200,30])

isz = endsWith(ch_table.name,'Z');
isx = endsWith(ch_table.name,'X');
isy = endsWith(ch_table.name,'Y');

plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r')

arr_phantom_ax = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_axis(1),phantom_axis(2),phantom_axis(3));
arr_phantom_ax.LineWidth = 2;arr_phantom_ax.Color = 'k';arr_phantom_ax.LineStyle=':';
h_dip_pos = plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',10,...
    'MarkerFaceColor','k');
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    dipole_ori(1),dipole_ori(2),dipole_ori(3),3);
arr_true.LineWidth = 3;arr_true.Color = 'k';

sgtitle('"Ground truth" pos/ori and data fieldmap')

% Field vectors for single timepoint

ti=1;
h=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),(OPM_data_mfc(isz,ti)));
h.FaceColor = 'interp';h.EdgeColor = 'none';h.FaceAlpha = 0.5;colormap parula
view([200,30])
% for ti = 1:2399
field = OPM_data_mfc(isx,ti).*xv + OPM_data_mfc(isy,ti).*yv + OPM_data_mfc(isz,ti).*zv;

arrs = quiver3(px',py',pz',...
    field(:,1),field(:,2),field(:,3),1,'r','LineWidth',2);
clim(1e3.*[-4    4])

%% similarity between GT and field at each timepoint
for ti = 1:size(OPM_data_mfc, 2)
    % waitforbuttonpress

    field = OPM_data_mfc(isx,ti).*xv + OPM_data_mfc(isy,ti).*yv + OPM_data_mfc(isz,ti).*zv;
    similarity(ti) = abs(pdist2(field(:)',LF_mags_true(:)','cosine') -1);
end
%%
figure
histogram(abs(similarity),'Normalization','cdf')
hold on
legend('Field data over time')
xlabel("|Cosine similarity|")
ylabel('CDF')
set(gcf,'Color','w')
%%
N_folds = 6;
size(OPM_data_mfc,2)./N_folds;

Z = zeros(Ndips,N_folds);
for f_ind = 1:N_folds
    win(f_ind,:) = (1 + (f_ind - 1) * size(OPM_data_mfc, 2) / N_folds) : (f_ind * size(OPM_data_mfc, 2) / N_folds);
    fold_data = OPM_data_mfc(:,win(f_ind,:));
    triggers_f_fold = triggers_f(win(f_ind,:));
    C = cov(fold_data');
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
        Z(n,f_ind) = (w*Cr*w')./(w*w');
    end
    % find max and calculate VE
    n = find(Z(:,f_ind)==max(Z(:,f_ind)));
    peak_pos(f_ind,:) = sourcepos(n,:);
    peak_or_theta(f_ind,:) = Src_Or_theta(n,:)./100;
    peak_or_phi(f_ind,:) = Src_Or_phi(n,:)./100;
    peak_or_R(f_ind,:) = Src_Or_R(n,:)./100;
    % ROT = [peak_or_theta;peak_or_phi;peak_or_R]';
    % ROT = [peak_or_theta;peak_or_phi;peak_or_R]';
    this_L = Lead_fields(:,:,n);
    % this_L = lead_fields_shell_xyz(:,:,n);

    iPower_v = this_L'*Cr_inv*this_L;
    [v,d] = svd(iPower_v);
    [~,id] = min(diag(d));
    lopt = this_L*v(:,id); % turn to nAm amplitude
    w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
    if size(v,1) == 3
        peak_ori = (peak_or_theta.*v(1,id) + peak_or_phi.*v(2,id) + peak_or_R.*v(3,id));
    else
        peak_ori = (peak_or_theta.*v(1,id) + peak_or_phi.*v(2,id));
    end
    % peak_ori = ROT*v(id,:)';

    VE = (w*fold_data)./sqrt(w*w');
    VE_z = (VE - mean(VE))./std(VE);
    % VE_z = -(VE./(std((VE))));

    triggers_z = (triggers_f_fold - mean(triggers_f_fold))./std(triggers_f_fold);
    % triggers_z = triggers;

    % time([1:300,end-300:end])=[];
    VE_z([1:300,end-300:end])=[];
    triggers_z([1:300,end-300:end])=[];
    %
    % figure
    % plot(time,VE_z)
    % hold on
    % plot(time,triggers_z)

    fftwin_s = 1;
    [pxx_VE(:,f_ind),~] = pwelch(VE_z,fftwin_s*f,[],[],f);
    [pxx_trig(:,f_ind),fxx] = pwelch(triggers_z,fftwin_s*f,[],[],f);
    spec_rms_diff(f_ind) = sqrt(mean((pxx_VE(:,f_ind) - pxx_trig(:,f_ind)).^2));

    figure
    plot(fxx,pxx_VE(:,f_ind))
    hold on
    plot(fxx,pxx_trig(:,f_ind))
    legend('VE', 'Trigger signal')
    xlim([0,60])
    
    similiarity_folds(f_ind) = mean(abs(similarity(win(f_ind,triggers_z>1.5))));
    
    correlation_folds(f_ind) = corr(VE_z',triggers_z');
    figure
    plot(VE_z)
    hold on
    plot(triggers_z)
    xlim([1639.2    2546.4])

    fftwin_s = 10;
    [coh(:,f_ind), fxx]  = mscohere(VE_z,triggers_z,fftwin_s*f,[],[],f);
%     [coh(:,f_ind), fxx]  = mscohere(VE_z,triggers_z,[],[],[],f);
    figure
    plot(fxx,coh(:,f_ind))
    xlim([0,60])
    

end

fprintf('RMS Spectral difference (VE vs TRIG) = %1.5f +/- %1.5f (AU)\n',mean(spec_rms_diff),std(spec_rms_diff));
fprintf('|Cosine Similarity| = %1.5f +/- %1.5f\n',mean(similiarity_folds),std(similiarity_folds));
% fprintf('Angular difference = %1.5f +/- %1.5f\n',acosd(mean(similiarity_folds)),abs(acosd(mean(similiarity_folds)+std(similiarity_folds)) - acosd(mean(similiarity_folds))))
fprintf('Correlation between trig and VE = %1.5f +/- %1.5f \n',mean(correlation_folds),std(correlation_folds));

figure
ci_coh = prctile(coh, [25, 75, 50],2);
plot(fxx,median(coh,2),'Color','k','LineWidth',2,'DisplayName','Coherence')
hold on
ciplot(ci_coh(:,1),ci_coh(:,2),fxx,[0,0,0]);
xlim([0,60])

xlabel('Frequency (Hz)')
ylabel('$|COH|^2$')

%% plot resulting loc
ff2 = figure;
load meshes.mat
ff2.Color = 'w';
ax_tru = subplot(122);
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
view([200,30])

isz = endsWith(ch_table.name,'Z');
isx = endsWith(ch_table.name,'X');
isy = endsWith(ch_table.name,'Y');

plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r')

arr_phantom_ax = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_axis(1),phantom_axis(2),phantom_axis(3));
arr_phantom_ax.LineWidth = 2;arr_phantom_ax.Color = 'k';arr_phantom_ax.LineStyle=':';
h_dip_pos = plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',10,...
    'MarkerFaceColor','k');
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    dipole_ori(1),dipole_ori(2),dipole_ori(3),3);
arr_true.LineWidth = 3;arr_true.Color = 'k';

% field_mag = sqrt(true_LF(isx).^2 + true_LF(isy).^2 + true_LF(isz).^2);
field_mag = true_LF(isz);
DT= boundary(S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),0);
h=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),field_mag);
h.FaceColor = 'interp';h.EdgeColor = 'none';h.FaceAlpha = 0.5;colormap parula
clim([-1,1].*1e-14);%colormap hot

ax_res = subplot(121);
copyobj(ax_tru.Children,ax_res);
view([200,30])

axis equal off
hold on
h.Visible = 0;

a2 = subplot(122);

plot3(peak_pos(f_ind,1),peak_pos(f_ind,2),peak_pos(f_ind,3),'b*','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','b')
arr_opt = quiver3(peak_pos(f_ind,1),peak_pos(f_ind,2),peak_pos(f_ind,3),...
    peak_ori(f_ind,1),peak_ori(f_ind,2),peak_ori(f_ind,3),2);
arr_opt.LineWidth = 3;arr_opt.Color = 'b';
quiver3(peak_pos(f_ind,1),peak_pos(f_ind,2),peak_pos(f_ind,3),peak_or_theta(f_ind,1), peak_or_theta(f_ind,2), peak_or_theta(f_ind,3),2,'r')
quiver3(peak_pos(f_ind,1),peak_pos(f_ind,2),peak_pos(f_ind,3),peak_or_phi(f_ind,1), peak_or_phi(f_ind,2), peak_or_phi(f_ind,3),2,'g')
quiver3(peak_pos(f_ind,1),peak_pos(f_ind,2),peak_pos(f_ind,3),peak_or_R(f_ind,1), peak_or_R(f_ind,2), peak_or_R(f_ind,3),2,'b')

dist_in_mm = sqrt(sum((peak_pos - dipole_pos).^2,2))*1000;

res = sqrt(sum((sourcepos(1,:) -  sourcepos(2,:)).^2))*1000;

angle_diffs = acosd(dot(peak_ori,repmat(dipole_ori,N_folds,1),2)./(repmat(vecnorm(dipole_ori,2,2),N_folds,1).*vecnorm(peak_ori,2,2)));

fprintf('Peak to truth distance = %1.3f +/- %1.3f mm (%1.1f mm res)\n',mean(dist_in_mm),std(dist_in_mm),res);
fprintf('Peak to truth orientation diff = %1.3f +/- %1.3f degs\n',mean(angle_diffs),std(angle_diffs));
%%
h2=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),...
    S.sensor_info.pos(isz,3),lopt(isz));
h2.FaceColor = 'interp';h2.EdgeColor = 'none';h2.FaceAlpha = 0.5;colormap parula
% clim([-1,1].*1e-15)
%% Z stat map
fg0 = figure;
fg0.Units = 'centimeters';
fg0.Position = [[17.9652 14.7373 4 4.5]];
lw = 2;
% subplot(2,3,[1 2 4 5])
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
s_plt=scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,mean(Z,2),'filled');
colormap hot;
% caxis([cax_val]);
% cb = colorbar;cb.Label.String = 'Zstat';
axis equal
view([-90 0])
xlim([-0.0000    0.2306])
ylim([-0.12    0.090])
zlim([-0.0838    0.1272])


fg = gcf;
fg.Color = 'w';
ax = gca;
ax.FontSize = 14;

%plot phantom location
plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',5,'MarkerFaceColor','r')

arr_phantom_ax = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_axis(1),phantom_axis(2),phantom_axis(3));
arr_phantom_ax.LineWidth = 2;arr_phantom_ax.Color = 'k';arr_phantom_ax.LineStyle=':';
h_dip_pos = plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',10,...
    'MarkerFaceColor','k');
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    dipole_ori(1),dipole_ori(2),dipole_ori(3),3);
arr_true.LineWidth = 3;arr_true.Color = 'k';

plot3(peak_pos(:,1),peak_pos(:,2),peak_pos(:,3),'g*','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','g')
arr_opt = quiver3(peak_pos(:,1),peak_pos(:,2),peak_pos(:,3),...
    peak_ori(:,1),peak_ori(:,2),peak_ori(:,3),4);
arr_opt.LineWidth = lw;arr_opt.Color = 'b';

view([-90 0])
xlim([-0.0000    0.2306])
ylim([-0.12    0.090])
zlim([-0.0838    0.1272])

view([0 90])
xlim([-0.2306    0.2306])
ylim([-0.12    0.090])
zlim([-0.0838    0.06575])
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_compare = figure;
lw = 1.5;
f_compare.Color = 'w';f_compare.Units = 'centimeters';
f_compare.Position = [37 13.6260 8 4];
ax_tru = subplot(122);
% head and brain meshes
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none');
ax_tru.Children(1).Visible = 0;
hold on
view([200,30])
% phantom slot red sphere
plot3(helmet_config.Px(phantom_slot_z),helmet_config.Py(phantom_slot_z),...
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',3,'MarkerFaceColor','r')
% dotted line to phantom position
arr_phantom_ax = quiver3(phantom_slot_pos(1),phantom_slot_pos(2),phantom_slot_pos(3),...
    phantom_axis(1),phantom_axis(2),phantom_axis(3));
arr_phantom_ax.LineWidth = 3;arr_phantom_ax.Color = 'k';arr_phantom_ax.LineStyle=':';
% ground truth phantom position
h_dip_pos = plot3(dipole_pos(1),dipole_pos(2),dipole_pos(3),'go','MarkerSize',5,...
    'MarkerFaceColor','k');
% ground truth phantom orientation
arr_true = quiver3(dipole_pos(1),dipole_pos(2),dipole_pos(3),...
    dipole_ori(1),dipole_ori(2),dipole_ori(3),3);
arr_true.LineWidth = 3;arr_true.Color = 'k';

% Field in Z
field_mag = true_LF(isz);
DT= boundary(S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),0);
h=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),field_mag);
h.FaceColor = 'interp';h.EdgeColor = 'none';h.FaceAlpha = 0.5;colormap parula
clim([-1,1].*0.5e-14);%colormap hot
axis equal


ax_res = subplot(121);
ax_res.CameraViewAngleMode = 'manual';
axis off
axis equal
view([200,30])
copyobj(ax_tru.Children,ax_res);
clim([-1,1].*5e-15);%colormap hot
hold on
arrs = quiver3(px',py',pz',...
    field(:,1),field(:,2),field(:,3),1,'r','LineWidth',lw);

hold on
h.Visible = 0;

subplot(ax_tru);

view([200,30])

isz = endsWith(ch_table.name,'Z');
isx = endsWith(ch_table.name,'X');
isy = endsWith(ch_table.name,'Y');


OPM_data_mfc_mean = mean(OPM_data_mfc(:,triggers<-3),2);
% OPM_data_mfc_mean = mean(OPM_data_mfc(:,corrs>0.5),2);
field  = (OPM_data_mfc_mean(isx).*xv + OPM_data_mfc_mean(isy).*yv + OPM_data_mfc_mean(isz).*zv);

h=trisurf(DT,S.sensor_info.pos(isz,1),S.sensor_info.pos(isz,2),S.sensor_info.pos(isz,3),-(OPM_data_mfc_mean(isz)));
h.FaceColor = 'interp';h.EdgeColor = 'none';h.FaceAlpha = 0.5;colormap parula
view([200,30])
field_lf = -(true_LF(isx).*xv + true_LF(isy).*yv + true_LF(isz).*zv);

arrs_lf = quiver3(px',py',pz',...
    field_lf(:,1),field_lf(:,2),field_lf(:,3),1,'b','LineWidth',lw);

clim(5e3.*[-1    1])

for f_ind = 1:N_folds
    fold_data = OPM_data_mfc(:,win(f_ind,:));
    triggers_f_fold = triggers_f(win(f_ind,:));

    OPM_data_mfc_mean = mean(fold_data(:,triggers_f_fold < -0.6*max(abs(triggers_f_fold))),2);
    field_fold  = (OPM_data_mfc_mean(isx).*xv + OPM_data_mfc_mean(isy).*yv + OPM_data_mfc_mean(isz).*zv);
    
    field_corr_fold(f_ind) = corr(field_lf(:),field_fold(:));
end
fprintf('Correlation between FL and fold average field @ 60 %% of max amp = %1.5f +/- %1.5f \n',mean(field_corr_fold),std(field_corr_fold));


views = [117,34;...
    170,15;...
    0,70;...
    200,30];


axpos = ax_tru.Position;
offset=0.02;
cross_pos = [ax_tru.XLim(1)-offset,ax_tru.YLim(2)+offset,ax_tru.ZLim(1)-offset];
% plot3([cross_pos(1),cross_pos(1)],[cross_pos(2),cross_pos(2)-0.05],[cross_pos(3),cross_pos(3)],'g','LineWidth',2)
% plot3([cross_pos(1),cross_pos(1)+0.05],[cross_pos(2),cross_pos(2)],[cross_pos(3),cross_pos(3)],'r','LineWidth',2)
% plot3([cross_pos(1),cross_pos(1)],[cross_pos(2),cross_pos(2)],[cross_pos(3),cross_pos(3)+0.05],'b','LineWidth',2)

quiver3([cross_pos(1)],[cross_pos(2)],[cross_pos(3)],...
    [cross_pos(1)-cross_pos(1)],[cross_pos(2)-cross_pos(2)-0.05],[cross_pos(3)-cross_pos(3)],...
    'g','LineWidth',2,'MaxHeadSize',0.7)
quiver3([cross_pos(1)],[cross_pos(2)],[cross_pos(3)],...
    [cross_pos(1)-cross_pos(1)+0.05],[cross_pos(2)-cross_pos(2)],[cross_pos(3)-cross_pos(3)],...
    'r','LineWidth',2,'MaxHeadSize',0.7)
quiver3([cross_pos(1)],[cross_pos(2)],[cross_pos(3)],...
    [cross_pos(1)-cross_pos(1)],[cross_pos(2)-cross_pos(2)],[cross_pos(3)-cross_pos(3)+0.05]...
    ,'b','LineWidth',2,'MaxHeadSize',0.7)

txt_offset = 0.07;
tr = text(cross_pos(1)+txt_offset,cross_pos(2),cross_pos(3),'R');
tp = text(cross_pos(1),cross_pos(2)-txt_offset,cross_pos(3),'P');
ts = text(cross_pos(1),cross_pos(2),cross_pos(3)+txt_offset,'S');
ax_tru.Position = axpos;
ax_res.XLim = ax_tru.XLim;
ax_res.YLim = ax_tru.YLim;
ax_res.ZLim = ax_tru.ZLim;

%
% for vi = 1
for vi = 1:length(views)
    set([ax_res, ax_tru],'View',views(vi,:))
    if vi ==2;tp.Visible = 0;end

    if vi == 3
        ts.Position(3) =     0.020;
        ts.Position(1) =     -0.1425;
        tr.Position(3) =     -0.0220;
        tr.Position(1) =     -0.0725;
    end
    input('next')
    tp.Visible = 1;
    ts.Position(3) =     -0.0220;
    ts.Position(1) =     -0.1225;
    tr.Position(3) =     -0.0920;
    tr.Position(1) =     -0.0525;
end

% sgtitle('"Ground truth" leadfield (blue) and data fieldmap in Z channels (red; average)')

%%


