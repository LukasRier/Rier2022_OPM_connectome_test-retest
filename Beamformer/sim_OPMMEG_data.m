function [sim_outs] = sim_OPMMEG_data(fwd_method,S,do_fake_data)
%% Function to simulate MEG data with inputted source and sensor positions and orientations
% [sim_outs] = sim_OPMMEG_data(fwd_method,S,do_fake_data)
% [INPUTS]
% fwd_method      = 'single' (single sphere method).
%                 = 'local' (local spheres method).
%                 = 'BEM'.
%                 = 'shell' (single shell method).
% S               = parameter structure.
% do_fake_data    = Set to 1 to simulate data. Otherwise, you will get the
%                   lead fields only for a particular sesnor location and
%                   orientation.
%
% For all models
% S.mri_file      = .mri file location.
% S.source_info   = source information structure (.pos & .ors both
%                   Nsourcesx3)
% S.sensor_info   = sensor information structure (.pos & .ors both
%                   Nchansx3).
% S.Source_amp    = Source amplitude (in fT, typically set as 1).
% S.Noise_amp     = Noise amplitude (in fT, typically of order of 10).
% S.f             = Sample rate.
% S.time          = Length of signal (seconds).
%
% For Local Spheres
% S.skanect_file  = Coregistered file to save local spheres to.
% S.rad_ors       = Radial orientation associated with each position. This
%                   should be Nchansx3. No matter the channel axis, use the
%                   Z orientation.
%
% For BEM
% S.skanect_file  = Coregistered file to save volume and lead fields to.
%
%
% [OUTPUTS]
% sim_outs.LF      = leadfield matrix (Nchans x 3 x M). In XYZ format (i.e. a
%                   dipole at [1 0 0], [0 1 0], and [0 0 1]).
% sim_outs.Weights = Weights for each source location (Nchans x M).
% sim_outs.T_stat  = T stat for each source location.
% sim_outs.Z_stat  = Z stat for each source location.

[path.mri, filename.mri] = fileparts(S.mri_file);
path.mri = [path.mri filesep];
try
    mri = ft_read_mri([path.mri filename.mri '.mri']);
catch
    try
        mri = ft_read_mri([path.mri filename.mri '.nii']);
    catch
        error('MRI file is not in .mri or .nii format')
    end
end

if ~exist([path.mri 'meshes.mat'])
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
    cfg = [];
    cfg.tissue = {'brain','skull','scalp'};
    cfg.numvertices = [3000 2500 2500];
    mesh1 = ft_prepare_mesh(cfg,segmentedmri);
    for n = 1:size(mesh1,2)
        meshes(n).pnt = mesh1(n).pos;
        meshes(n).tri = mesh1(n).tri;
        meshes(n).unit = mesh1(n).unit;
        meshes(n).name = cfg.tissue{n};
        meshes(n) = ft_convert_units(meshes(n),'m');
    end
    save([path.mri,files.meshes],'meshes','segmentedmri')
else
    load([path.mri 'meshes.mat'])
end
origin = mean(meshes(1).pnt);

sourcepos = S.source_info.pos;
sourceor = S.source_info.ors;

%%%%%%%%%%% Single Sphere %%%%%%%%%%%
if strcmpi(fwd_method,'single')
    for ii = 1:size(sourcepos,1)
        % compute the lead fields
        [LF(:,:,ii)] = Forward_OPM_single_sphere_XYZ_dir(sourcepos(ii,1), sourcepos(ii,2), sourcepos(ii,3),...
            S.sensor_info.pos,S.sensor_info.ors,origin);
        fprintf('Done %u/%u \n',ii,length(sourcepos))
    end
    
    for ii = size(sourcepos,1):-1:1
        l_input(:,ii) = LF(:,:,ii)*sourceor(ii,:)';
    end
    clc
    %%%%%%%%%%% Local Spheres %%%%%%%%%%%
elseif strcmpi(fwd_method,'local')
    
    [path.skanect, filename.skanect] = fileparts(S.skanect_file);
    path.skanect = [path.skanect filesep];
    
    % Create spheres
    if exist([path.skanect,'local_spheres_info_',filename.skanect(end-5:end),'.mat'])
        load([path.skanect,'local_spheres_info_',filename.skanect(end-5:end),'.mat']);
    else
        brain_mesh.faces = meshes(1).tri;
        brain_mesh.vertices = meshes(1).pnt;
        sensors.pos = S.sensor_info.pos;
        sensors.ors = S.rad_ors;
        [so1, sr1] = create_local_spheres(brain_mesh,sensors);
        
        so = so1;
        sr = sr1;
        
        save([path.skanect,'local_spheres_info_',filename.skanect(end-5:end),'.mat'],...
            'so','sr','brain_mesh','sensors');
    end
    
    % Plot spheres
    figure
    ft_plot_mesh(meshes(1),'facecolor',[.5 .5 .5],'facealpha',0.3,'edgecolor','none')
    hold on
    [sx1,sy1,sz1] = sphere;
    axis equal
    s1 = [];
    s = [];
    ccc = jet(size(sensors.pos,1));
    for n = 1:size(sensors.pos,1)
        s = scatter3(sensors.pos(n,1),sensors.pos(n,2),sensors.pos(n,3),'MarkerFaceColor',[ccc(n,:)]);
        sx = sx1.*sr(n) + so(n,1);
        sy = sy1.*sr(n) + so(n,2);
        sz = sz1.*sr(n) + so(n,3);
        s1 = surf(sx,sy,sz,'Edgecolor','none','Facecolor',[ccc(n,:)],'FaceAlpha',0.3);
    end
    
    origins = so;
    origin_mean = mean(origins); % Mean of local sphere origins
    
    % Lead fields
    %     TR = triangulation(brain_mesh.faces,brain_mesh.vertices);
    %     vn = vertexNormal(TR);
    for ii = 1:size(sourcepos,1)
        % Find nearest surface normal if using it for tangential calculation
        source2mesh_dist = pdist2(sourcepos(ii,:),brain_mesh.vertices);
        [~,s2m_idx] = min(source2mesh_dist);
        %         nrst_norm = vn(s2m_idx,:); % Nearest surface normal to the dipole location
        %         dip_loc = sourcepos(ii,:);
        % Find lead fields and compute BF weights
        [LF(:,:,ii)] = Forward_OPM_local_spheres_XYZ_dir(sourcepos(ii,1),sourcepos(ii,2),sourcepos(ii,3),...
            S.sensor_info.pos,S.sensor_info.ors,origins);
        fprintf('Done %u/%u \n',ii,length(sourcepos))
    end
    
    for ii = size(sourcepos,1):-1:1
        l_input(:,ii) = LF(:,:,ii)*sourceor(ii,:)';
    end
    
elseif strcmpi(fwd_method,'BEM')
    % Create BEM model
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
        disp('Loaded BEM meshes')
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
    coils.p = S.sensor_info.pos;
    coils.n = S.sensor_info.ors;
    
    % Create BEM model
    if ~ismember('TBvol',who('-file',[S.skanect_file]))
        TBvol = create_BEM(bmeshes,coils,ci);
        save([S.skanect_file],'TBvol','-append')
    else
        load([S.skanect_file],'TBvol')
        disp('Loaded Volume')
    end
    
    % Create lead fields (takes a fair while...)
    if ~ismember('LF_BEM',who('-file',[S.skanect_file]))
        %         LF_BEM = zeros(size(coils.p,1),3,size(sourcepos,1));
        for n = 1:size(sourcepos,1)
            disp((n./size(sourcepos,1))*100)
            dip_loc = sourcepos(n,:);
            [LF_BEM(:,:,n)] = Forward_OPM_BEM_XYZ_dir(dip_loc,bmeshes,coils,TBvol);
            
        end
        disp('Saving...')
        save([S.skanect_file],'LF_BEM','-append')
        disp('Done!')
    else
        disp('Loading lead fields...')
        load([S.skanect_file],'LF_BEM')
        disp('Done!')
    end
    clc
    LF = LF_BEM;
    
    for ii = size(sourcepos,1):-1:1
        l_input(:,ii) = LF(:,:,ii)*sourceor(ii,:)';
    end    
    
elseif strcmpi(fwd_method,'shell')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));
    
    [LF_single_shell,L_reshaped] = Forward_OPM_shell_FT([path.mri,filename.mri],sensor_info,sourcepos);
    LF1 = L_reshaped.*1e-9;
    inside_idx = find(LF_single_shell.inside);
    LF = zeros(size(LF1));
    for n = size(sourcepos,1):-1:1
        if LF_single_shell.inside(n)
            nn = find(inside_idx == n);
            LF(:,:,n) = LF1(:,:,nn);
        else
            LF(:,:,n) = zeros(size(LF1(:,:,1)));
        end
        
    end
    clear LF1 L_reshaped
    for ii = size(sourcepos,1):-1:1
        l_input(:,ii) = LF(:,:,ii)*sourceor(ii,:)';
    end    
    
end

sim_outs.LF_XYZ = LF;
sim_outs.LF = l_input;

if do_fake_data
        time_course = S.Source_amp*randn(size(linspace(0,S.time,S.f*S.time)));
        sim_outs.source_time_course = time_course;
        fake_sens_data1 = sim_outs.LF*time_course;
        % add sensor noise
        sens_noise = ([S.Noise_amp*1e-15])*randn(size(fake_sens_data1));
        fake_sens_data = fake_sens_data1 + sens_noise;
        sim_outs.sim_data = fake_sens_data;
end

end

function [tanu, tanv] = calctangent_all_bf(RDip)
%% Same as usual calctangent but put here for ease

x=RDip(1);
y=RDip(2);
z=RDip(3);
r=sqrt(x*x+y*y+z*z);

if (x==0) && (y==0)
    tanu(1)=1.0; tanu(2)=0; tanu(3)=0;
    tanv(1)=0; tanv(2)=1.0; tanv(3)=0;
else
    RZXY= -(r-z)*x*y;
    X2Y2= 1/(x*x+y*y);
    
    tanu(1)= (z*x*x + r*y*y) * X2Y2/r;
    tanu(2)= RZXY * X2Y2/r;
    tanu(3)= -x/r;
    
    tanv(1)= RZXY * X2Y2/r;
    tanv(2)= (z*y*y + r*x*x) * X2Y2/r;
    tanv(3)= -y/r;
end
end