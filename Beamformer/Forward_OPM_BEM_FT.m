function [Lead_fields] = Forward_OPM_BEM_FT(mri_file,sensor_info,skanect_file,sourcepos)
% Use FieldTrip to create a BEM lead field model
% mri_file = .mri file
% sensor_info = sensor positions (Nx3), orientations(Nx3), and labels (Nx1)
% skanect_file = where to save/load things to/from

%% Load OPM sensor positions
grad.coilpos = sensor_info.pos.*100;
grad.coilori = sensor_info.ors;
grad.label = sensor_info.label;
grad.chanpos = grad.coilpos;
grad.units = 'cm';

%% Create head model
[path.mri, filename.mri] = fileparts(mri_file);
path.mri = [path.mri filesep];
mri = ft_read_mri([path.mri filename.mri '.mri']);
if ~exist([path.mri 'segmentedmri.mat'])
    disp('Segmenting MRI...')
    cfg           = [];
    cfg.output    = {'brain','skull','scalp'};
    segmentedmri  = ft_volumesegment(cfg, mri);
    save([path.mri 'segmentedmri.mat'],'segmentedmri')
    disp('Done!')
else
    disp('Loading Segmented MRI...')
    load([path.mri 'segmentedmri.mat'])
    disp('Done!')
end

% seg_i = ft_datatype_segmentation(segmentedmri,'segmentationstyle','indexed');
% cfg = [];
% cfg.funparameter = 'seg';
% cfg.funcolormap  = gray(4); % distinct color per tissue
% cfg.location     = 'center';
% cfg.atlas        = seg_i;
% ft_sourceplot(cfg, seg_i);

cfg = [];
cfg.tissue = {'brain','skull','scalp'};
cfg.numvertices = [3000 2000 1000];
mesh_bem = ft_prepare_mesh(cfg,segmentedmri);
figure 
ft_plot_mesh(mesh_bem(1),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'skin','facealpha',0.5,'edgealpha',0.1)
ft_plot_mesh(mesh_bem(2),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'skin','facealpha',0.5,'edgealpha',0.1)
ft_plot_mesh(mesh_bem(3),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'skin','facealpha',0.5,'edgealpha',0.1)
      
cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);

%% Visualise
headmodel = ft_convert_units(vol, 'cm');

figure
ft_plot_sens(grad, 'style', '*b');
hold on
ft_plot_vol(headmodel);

%% Create sourcemodel
cfg             = [];
cfg.grad        = grad;
cfg.headmodel   = headmodel;
cfg.resolution  = 0.5;
cfg.inwardshift = -1;
sourcemodel     = ft_prepare_sourcemodel(cfg);

figure;
ft_plot_sens(grad, 'style', '*b');
ft_plot_vol(headmodel, 'edgecolor', 'none'); alpha 0.4;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));

%% Create leadfield
% cfg                  = [];
% cfg.grad             = grad;  % gradiometer distances
% cfg.headmodel        = headmodel;   % volume conduction headmodel
% cfg.sourcemodel      = sourcemodel;
% cfg.channel          = {'all'};
% cfg.singleshell.batchsize = 2000;
% lf                   = ft_prepare_leadfield(cfg);

cfg                 = [];
cfg.channel         = grad.label; % ensure that rejected sensors are not present
cfg.grad            = grad;
cfg.headmodel       = headmodel;
cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
cfg.resolution = 1;   % use a 3-D grid with a 1 cm resolution
cfg.unit       = 'cm';
cfg.tight      = 'yes';
[grid] = ft_prepare_leadfield(cfg);

end

