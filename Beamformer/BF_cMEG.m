%% Beamformer for cMEG data
% Before running this script extract and export your data using the 
% extract_cMEG_data app into the OPM_vars structure. Load OPM_vars into
% the workspace. Will also need coregistration information from
% apply_meshlab_trans_TRIAX function.
%
% This script uses run_beamformer_TRIAX function to beamform OPM data using
% single sphere, local spheres, BEM or single shell forward models.

clc;
addpath('R:\Matlab_files\tools\');
addpath R:\Matlab_files\fieldtrip-20190212\external\spm12
addpath('R:\Matlab_files\analysis_gui\app_functions');
if isempty(which('ft_defaults'))
    addpath('R:\Matlab_files\fieldtrip-20190212')
    ft_defaults
end

% Participant information
file_info = split(OPM_vars.filestr,'\');
date = char(file_info(4));

% MRI
path.mri = ['R:\data\MRI\' OPM_vars.Vol_nr '\'];
files.mri = [OPM_vars.Vol_nr '.mri'];
files.meshes = 'meshes.mat';

% Data
path.data = fileparts(OPM_vars.filestr);

% Results folder
path.results = [path.data,'\BF_outputs\'];
if ~exist(path.results)
    mkdir(path.results)
end

% Skanect
cd(['R:\data\skanect\' date '\' OPM_vars.Vol_nr '\']);
[files.skanect,path.skanect] = uigetfile('*.mat');
load([path.skanect,files.skanect]);

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
f1.Name = sprintf('Coverage of %s',num2str(OPM_vars.Vol_nr));
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3))
view([180 0])
fig = gcf;
fig.Color = [1,1,1];

%% Select only the sensors in use from coreg file
% need all these bits to correctly duplicate the origins for local spheres
S.sensor_info.axis = OPM_vars.QZFM_data.axis;
S.sensor_info.Nlocs = OPM_vars.Nlocs;
S.sensor_info.loc_names = OPM_vars.loc_names;
S.sensor_info.OPM_data_struct = OPM_vars.OPM_data_struct;

S.sensor_info.pos = [];
S.sensor_info.ors = [];
for kk = 1:OPM_vars.Nlocs
    Loc_ind = OPM_vars.OPM_data_struct.(['Sensor_',OPM_vars.loc_names{kk}]);
    S.sensor_info.pos = cat(1,S.sensor_info.pos,repmat(sens_info.pos(OPM_vars.pos_no(kk),:),[length(Loc_ind) 1]));
    pos_axes = OPM_vars.QZFM_data.axis(Loc_ind);
    for ii = 1:length(Loc_ind)
        if strcmpi(pos_axes(ii),'X')
            S.sensor_info.ors = cat(1,S.sensor_info.ors,sens_info.ors_X(OPM_vars.pos_no(kk),:));
        elseif strcmpi(pos_axes(ii),'Y')
            S.sensor_info.ors = cat(1,S.sensor_info.ors,sens_info.ors_Y(OPM_vars.pos_no(kk),:));
        elseif strcmpi(pos_axes(ii),'Z')
            S.sensor_info.ors = cat(1,S.sensor_info.ors,sens_info.ors_Z(OPM_vars.pos_no(kk),:));
        end
    end
end

figure(1); hold on;
plot3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),'o');
quiver3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),...
    S.sensor_info.ors(:,1),S.sensor_info.ors(:,2),S.sensor_info.ors(:,3));
text(sens_info.pos(OPM_vars.pos_no,1),sens_info.pos(OPM_vars.pos_no,2),sens_info.pos(OPM_vars.pos_no,3),OPM_vars.loc_names);

%% Make covariance matrix and regularise
% Regularisation parameter (%)
mu = 0.01;
% Data in teslas
OPM_data_T = OPM_vars.OPM_data*1e-15;
OPM_data_f_T = OPM_vars.OPM_dataf*1e-15;
OPM_data_T_mat = reshape(OPM_data_T,[OPM_vars.Nchans,OPM_vars.duration*OPM_vars.f,OPM_vars.Ntrials]);
OPM_data_f_T_mat = reshape(OPM_data_f_T,[OPM_vars.Nchans,OPM_vars.duration*OPM_vars.f,OPM_vars.Ntrials]);
% Covariance matrix
S.C = cov(OPM_data_f_T);
maxEV = max(svd(S.C));
S.Noise_Cr = min(svd(S.C)).*eye(size(S.C));
S.C = S.C + mu.*maxEV.*eye(size(S.C));
condnr = cond(S.C);
S.Cinv = inv(S.C);
figure;
imagesc(S.C)
axis square
title(sprintf('\n mu = %.3f, cond(C) = %.2f',...
    mu,condnr))

% Active and Control covariance
inputs2 = inputdlg({'ON window:','OFF window:'},'Inputs',[1 30],{'0 5','5 10'});
ONwin   = str2num(inputs2{1});
OFFwin  = str2num(inputs2{2});
win = abs(ONwin(1)-ONwin(2));
OPM_data_T_mat_A = OPM_data_f_T_mat(:,ONwin(1)*OPM_vars.f+1:ONwin(2)*OPM_vars.f,:);
OPM_data_T_mat_C = OPM_data_f_T_mat(:,OFFwin(1)*OPM_vars.f+1:OFFwin(2)*OPM_vars.f,:);
OPM_data_T_A = reshape(OPM_data_T_mat_A,OPM_vars.Nchans,OPM_vars.f*OPM_vars.Ntrials*win);
OPM_data_T_C = reshape(OPM_data_T_mat_C,OPM_vars.Nchans,OPM_vars.f*OPM_vars.Ntrials*win);
S.Ca = cov(OPM_data_T_A');
S.Cc = cov(OPM_data_T_C');
pause(0.05)

%% Organise for input to run_beamformer function

% Input window
inputs = inputdlg({'forward model: single, local, BEM, shell',...
    'do beamformer: 1 to beamform, 0 for lead fields only',...
    'max or min for T stat'},...
    'Inputs',[1 30],{'shell','1','max'});

fwd_method = inputs{1};
do_bf = str2double(inputs{2});
max_or_min = inputs{3};

% mean field correction
if isfield(OPM_vars,'MFC_multiplier')
    MFC = 1;
    S.M = OPM_vars.MFC_multiplier;
% there was a spelling mistake in the app at one point so to catch 'em all:
elseif isfield(OPM_vars,'MFC_mulitplier')
    MFC = 1;
    S.M = OPM_vars.MFC_mulitplier;
else
    MFC = 0;
end

% see help for run_beamformer_TRIAX but the following should cover
% everything you need to run any of the forward models:
S.mri_file = [path.mri files.mri];
S.skanect_file = [path.skanect files.skanect];

S.Z_inds = find(strcmpi(OPM_vars.QZFM_data.axis,'Z'));
S.Y_inds = find(strcmpi(OPM_vars.QZFM_data.axis,'Y'));
S.X_inds = find(strcmpi(OPM_vars.QZFM_data.axis,'X'));
% as per instructions in run_beamformer, if there are no tangential 
% channels make their inds equal to radial channel inds:
if isempty(S.Y_inds)
    S.Y_inds = S.Z_inds;
end
if isempty(S.X_inds)
    S.X_inds = S.Z_inds;
end

%% Run beamformer
clc; addpath('R:\Matlab_files\Beamformer\triax');
[bf_outs] = run_beamformer_TRIAX(fwd_method,sourcepos,S,1,'max',MFC);

%% Save stuff

file_path = [num2str(OPM_vars.hp) '_to_' num2str(OPM_vars.lp) 'Hz\'];
if ~exist([path.results file_path])
    mkdir([path.results file_path])
end
if ~exist([[path.results file_path],OPM_vars.Vol_nr,'_',fwd_method,'.mat'])
    save([[path.results file_path],OPM_vars.Vol_nr,'_',fwd_method,'.mat'],'mu','condnr','bf_outs')
else
    save([[path.results file_path],OPM_vars.Vol_nr,'_',fwd_method,'.mat'],'mu','condnr','bf_outs','-append')
end

%% Remember the upside down anatomical problem, fieldtrip can fix that too.
if ~exist([path.mri 'acpc_space\'])
    mkdir([path.mri 'acpc_space\'])
end

% CTF space is a total nightmare with MRI visualisation software, so we
% need to put it in a normal cordinate space, such as ACPC.
mri1 = ft_convert_coordsys(mri, 'acpc');

% Do the final flipping, and then write out hires anatomical.
mri2 = align_ijk2xyz(mri1);
if ~exist([path.mri 'acpc_space\anat.nii'])
    ft_write_mri([path.mri 'acpc_space\anat.nii'],mri2.anatomy,'dataformat',...
        'nifti','transform',mri2.transform);
end

% for the 4mm t-stat image, just doing the same thing as above never seems
% to work, so we need to downsample the ACPC space to 4mm and use all the
% transformation matrices to take the the 4mm CTF image into 4mm ACPC
cfg = [];
cfg.downsample = 4;
mri4 = ft_volumedownsample(cfg,mri1);

% This insane multiplication is required to get the 4mm ctf-space image into
% the 4mm-acpc space
kernel=mri4.transform*inv(mri1.transform)*inv(mri1.head2headOrig)*inv(brainDS.transform*inv(brain.transform));

image = brainDS;
image.anatomy(voxID) = cell2mat(bf_outs.T_stat);
image.transform=kernel*image.transform;

% Do the final flipping, and then write out lo-res functional image.
image2 = align_ijk2xyz(image);
Tstat = align_ijk2xyz(image);
ft_write_mri([[path.results file_path] [OPM_vars.Vol_nr,'_',fwd_method,'.nii']],image2.anatomy,'dataformat','nifti',...
    'transform',image2.transform);
image = brainDS;
image.anatomy(voxID) = cell2mat(bf_outs.T_stat);
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