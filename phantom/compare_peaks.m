
clear all

fname = 'R:\DRS-KidsOPM\phantom\20230309_164734_cMEG_Data\20230309_164734_meg.cMEG';
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
    load([fname(1:end-8),'ch_table.mat'],'ch_table')
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

load meshes.mat
%%
figure;

subplot(2,4,[1 2 3 5 6 7])
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
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
    helmet_config.Pz(phantom_slot_z),'ro','MarkerSize',15,'MarkerFaceColor','r',...
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
load beta_phantom_peak.mat beta_phantom_peak
load broad_band_phantom_peak.mat broad_band_phantom_peak
peak_pos_beta = beta_phantom_peak.pos
peak_pos_bb = broad_band_phantom_peak.pos;
fprintf('Beta peak to truth distance = %1.3f mm\n',sqrt(sum((dipole_pos-peak_pos_beta).^2))*1000);
plot3(peak_pos_beta(1),peak_pos_beta(2),peak_pos_beta(3),'bo','MarkerSize',8,'MarkerFaceColor','b',...
    'DisplayName','Beamformer Peak beta')
sc=2;

opt_ori_beta = beta_phantom_peak.ori;
arr_opt_beta = quiver3(peak_pos_beta(1),peak_pos_beta(2),peak_pos_beta(3),...
    opt_ori_beta(1),opt_ori_beta(2),opt_ori_beta(3),sc,'DisplayName','B.F. opt orientation | beta peak');
arr_opt_beta.LineWidth = 3;arr_opt_beta.Color = 'b';

fprintf('Noise peak to truth distance = %1.3f mm\n',sqrt(sum((dipole_pos-peak_pos_bb).^2))*1000);
plot3(peak_pos_bb(1),peak_pos_bb(2),peak_pos_bb(3),'co','MarkerSize',8,'MarkerFaceColor','c',...
    'DisplayName','Beamformer Peak noise')

fprintf('Noise peak to beta peak distance = %1.3f mm\n',sqrt(sum((peak_pos_beta-peak_pos_bb).^2))*1000);

opt_ori_broadb = broad_band_phantom_peak.ori;
arr_opt_broadb = quiver3(peak_pos_bb(1),peak_pos_bb(2),peak_pos_bb(3),...
    opt_ori_broadb(1),opt_ori_broadb(2),opt_ori_broadb(3),sc,'DisplayName','B.F. opt orientation | noise peak');
arr_opt_broadb.LineWidth = 3;arr_opt_broadb.Color = 'c';

% arr_th = quiver3(peak_pos_beta(1),peak_pos_beta(2),peak_pos_beta(3),...
%     peak_or_theta(1),peak_or_theta(2),peak_or_theta(3),sc,...
%     'DisplayName','L.F. \theta');
% arr_th.LineWidth = 2;arr_th.Color = 'c';
% 
% arr_ph = quiver3(peak_pos_beta(1),peak_pos_beta(2),peak_pos_beta(3),...
%     peak_or_phi(1),peak_or_phi(2),peak_or_phi(3),sc,'DisplayName','L.F. \phi');
% arr_ph.LineWidth = 2;arr_ph.Color = 'm';



ll = legend;ll.Position(1) = 0.6;
ax=gca;set(ax.Children(end-3:end),'HandleVisibility','off');
