restoredefaultpath
cleaning_only = 0;
close all
clear all
clc
%
sub = '010';
ses = '002';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project_dir = 'R:\OPMMEG\Projects\movie\';
% project_dir = 'F:\Rdrive\movie\';
project_dir = '/net/cador/data_local/Lukas/movie/';
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
path.noiseVEs = [datadir,'derivatives',filesep,'noiseVEs',filesep,'sub-',sub,filesep];
path.pow = [datadir,'derivatives',filesep,'POW',filesep,'sub-',sub,filesep];

if ~exist(path.ICA,'dir'); mkdir(path.ICA);end
if ~exist(path.cleaning,'dir'); mkdir(path.cleaning);end
if ~exist(path.VEs,'dir'); mkdir(path.VEs);end
if ~exist(path.noiseVEs,'dir'); mkdir(path.noiseVEs);end
if ~exist(path.pow,'dir'); mkdir(path.pow);end

path.data = [path.main,'meg',filesep];
path.noise = [path.main,'noise',filesep];
path.mri = [path.main,'anat',filesep];
files.mri = ['sub-',sub,'.nii'];
files.meshes = ['sub-',sub,'_meshes.mat'];
files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
files.noise = ['sub-',sub,'_task-noise'];
files.noiseVEs = ['sub-',sub,'_','run-',ses,'_noiseVE'];
files.pow = ['sub-',sub,'_','run-',ses,'_POW'];
files.voxlox = ['sub-',sub,'_voxlox.mat'];
files.channels = '_channels.tsv';
files.noise_channels = ['_channels.tsv'];
files.sens_order = '_sensor_order.mat';
files.ICA = [path.ICA,filename];
files.cleaning = [path.cleaning,filename];
S.mri_file = [path.meshes,'sub-',sub,'_',filesep,'x.nii']; % only need path for beamformer and single shell function. bit hacky, FIX
% global files
cd(path.data)
load([path.data,filename,'_meg.mat'],'fs','data','triggers')
load([datadir,'derivatives',filesep,'helmet',filesep,'M1p5_adult.mat'])
noisedata = load(sprintf('%s%s.mat',path.noise,files.noise));
noise_fs = noisedata.fs; noise_data = noisedata.data; clear noisedata
% basefilename = [path.main,'ses-',ses,'/Matt/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load MRI and meshes
mri = ft_read_mri([path.mri,files.mri]);
load([path.meshes,files.meshes],'meshes');
load([path.meshes,files.voxlox],'voxlox');


% Preproc
% Mean correct
data = data - mean(data,1);
noise_data = noise_data - mean(noise_data,1);
% Notch filter
for harms = [50,100,150,120]
    Wo = harms/(fs/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    if fs ~= noise_fs; error("Noise and data not sampled at the same frequency!");end
    disp('Applying Notch filter')
    data = filter(b,a,data,[],1);
    noise_data = filter(b,a,noise_data,[],1);
end
% bandpass filter for viewing
disp('Applying 1-150 Hz bandpass filter')

hp = 1;
lp = 150;
[b,a] = butter(4,2*[hp lp]/fs);
data_f = [filtfilt(b,a,data)]';
noise_data_f = [filtfilt(b,a,noise_data)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of bad channels
% Get rid of bad channels in data
ch_table = readtable([path.data,filename,files.channels],'FileType','text','Delimiter','tab');
ch_table.isx = endsWith(ch_table.name,'[X]');
ch_table.isy = endsWith(ch_table.name,'[Y]');
ch_table.isz = endsWith(ch_table.name,'[Z]');
ch_table.slot_no = zeros(height(ch_table),1);
% sanity check
if sum(sum(ch_table{:,11:13},2)) ~= height(ch_table)
    error('Channel orientation [x,y,z] labels might be wrong!')
end
    
noise_ch_table = readtable(sprintf('%s%s%s',path.noise,files.noise,files.noise_channels),'FileType','text','Delimiter','tab');
noise_ch_table.isx = endsWith(noise_ch_table.name,'[X]');
noise_ch_table.isy = endsWith(noise_ch_table.name,'[Y]');
noise_ch_table.isz = endsWith(noise_ch_table.name,'[Z]');
% sanity check
if sum(sum(noise_ch_table{:,11:13},2)) ~= height(noise_ch_table)
    error('Channel orientation [x,y,z] labels might be wrong!')
end
% Get rid of channels not present in either noise or MEG data
for ch_i = 1:height(ch_table)
    ind_in_noise_data = find(startsWith(noise_ch_table.name,ch_table.name(ch_i)));
    if isempty(ind_in_noise_data)
        ch_table.status(ch_i) = {'bad'}
    end
end
for ch_i = 1:height(noise_ch_table)
    ind_in_data = find(startsWith(ch_table.name,noise_ch_table.name(ch_i)));
    if isempty(ind_in_data)
        noise_ch_table.status(ch_i) = {'bad'}
    end
end
% sanity checks
if sum(sum(ch_table{:,11:13},2)) ~= height(ch_table)
    error('[ MEG recording ] : Channel orientation [x,y,z] labels might be wrong!')
end
if sum(sum(noise_ch_table{:,11:13},2)) ~= height(noise_ch_table)
    error('[ Noise recording ] : Channel orientation [x,y,z] labels might be wrong!')
end

% make sure all bad channels are excluded from noise and MEG data
bn = unique([noise_ch_table.name(strcmp(noise_ch_table.status, 'bad'));ch_table.name(strcmp(ch_table.status, 'bad'))]);
for ii = 1:length(bn) 
    ch_table.status(find(strcmp(ch_table.name,bn{ii}))) = {'bad'};
    noise_ch_table.status(find(strcmp(noise_ch_table.name,bn{ii}))) = {'bad'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot spectra to check for bad channels
[psd_data_all,fxx_d] = get_PSD(data_f,fs);

noise_data_f(:,1:1200) = [];noise_data_f(:,end-1200:end) = []; % get rid of edge artifacts
[psd_data_all_n,fxx_n] = get_PSD(noise_data_f,fs);

figure
phs = plot(fxx_n,psd_data_all_n,'r');set(gca,'YScale','Log');xlim([0,150]);hold on
for ii = 1:length(phs);phs(ii).DisplayName = noise_ch_table.name{ii};end
hold on
phs = plot(fxx_d,psd_data_all);set(gca,'YScale','Log');xlim([0,150]);
for ii = 1:length(phs);phs(ii).DisplayName = ch_table.name{ii};end
plotbrowser;



%% remove bad channels (noise)
[noise_ch_table] = Bad_Channels(noise_data_f',noise_ch_table,fs,hp,lp);

%% make sure all bad channels are excluded from noise and MEG data again
bn = unique([noise_ch_table.name(strcmp(noise_ch_table.status, 'bad'));ch_table.name(strcmp(ch_table.status, 'bad'))])
for ii = 1:length(bn) 
    ch_table.status(find(strcmp(ch_table.name,bn{ii}))) = {'bad'};
    noise_ch_table.status(find(strcmp(noise_ch_table.name,bn{ii}))) = {'bad'};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove bad channels (data)
[ch_table] = Bad_Channels(data_f',ch_table,fs,hp,lp);
%
%% make sure all bad channels are excluded from noise and MEG data again
bn = unique([noise_ch_table.name(strcmp(noise_ch_table.status, 'bad'));ch_table.name(strcmp(ch_table.status, 'bad'))])
for ii = 1:length(bn) 
    ch_table.status(find(strcmp(ch_table.name,bn{ii}))) = {'bad'};
    noise_ch_table.status(find(strcmp(noise_ch_table.name,bn{ii}))) = {'bad'};
end

%% Check final matrices are the same

writetable(ch_table,[path.data,filename,files.channels,'_new'],...
        'WriteRowNames',true,'Delimiter','tab','FileType','text')

ch_table = readtable([path.data,filename,files.channels,'_new'],...
        'Delimiter','tab','FileType','text');
    
writetable(noise_ch_table,sprintf('%s%s_%s%s_new',path.noise,files.noise,ses,files.noise_channels),...
        'WriteRowNames',true,'Delimiter','tab','FileType','text')
noise_ch_table = readtable(sprintf('%s%s_%s%s_new',path.noise,files.noise,ses,files.noise_channels),...
        'Delimiter','tab','FileType','text');
disp("Removing bad channels")
bad_chans_data = [find(startsWith(ch_table.status,'bad'))];
bad_chans_noise = [find(startsWith(noise_ch_table.status,'bad'))];

ch_table_ = ch_table;
data_f_ = data_f;
noise_ch_table_ = noise_ch_table;
noise_data_f_ = noise_data_f;


ch_table_(bad_chans_data,:) = [];
data_f_(bad_chans_data,:) = [];
noise_ch_table_(bad_chans_noise,:) = [];%clear noise_ch_table
noise_data_f_(bad_chans_noise,:) = [];

if size(data_f_,1) ~= size(noise_data_f_,1)
    error('Channel counts not the same!')
end

if ~all([ch_table_.name{:}] == [noise_ch_table_.name{:}])
    error('Noise and data channels not in same order!')
end
%%
[psd_data_all,fxx_d] = get_PSD(data_f_,fs);

[psd_data_all_n,fxx_n] = get_PSD(noise_data_f_,fs);

figure
phs = plot(fxx_n,psd_data_all_n,'r');set(gca,'YScale','Log');xlim([0,150]);hold on
for ii = 1:length(phs);phs(ii).DisplayName = noise_ch_table.name{ii};end
hold on
phs = plot(fxx_d,psd_data_all);set(gca,'YScale','Log');xlim([0,150]);
for ii = 1:length(phs);phs(ii).DisplayName = ch_table.name{ii};end
plotbrowser;