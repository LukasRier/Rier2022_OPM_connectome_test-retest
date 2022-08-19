function [bf_outs] = run_beamformer(fwd_method,sourcepos,S,do_bf,max_or_min,MFC)
%% Function to run whichever beamformer you choose
% [bf_outs] = run_beamformer(fwd_method,sourcepos,S,do_bf,max_or_min,MFC)
% [INPUTS]
% fwd_method      = 'single' (single sphere method).
%                 = 'local' (local spheres method).
%                 = 'BEM'.
%                 = 'shell' (single shell method).
% sourcepos       = dipole locations (M x 3).
% S               = parameter structure.
% do_bf           = set to 1 if you want to beamform, or 0 if you only want
%                   the leadfields.
% max_or_min      = look at the maximum or minimum T-Stat by setting 'Max'
%                   or 'Min' (or N/A if not used).
% MFC             = set to 1 if mean field correction is used.
%
% For all models
% S.mri_file      = .mri file location.
% S.sensor_info   = sensor information structure (.pos & .ors both
%                   Nchansx3).
% S.M             = The mean field correction matrix (used to multiply the
%                   lead fields).
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
% For beamforming (do all in Tesla)
% S.C             = Data covariance.
% S.Cinv          = Inverse of data covariance.
% S.Noise_Cr      = Noise covariance.
% S.Ca            = Active covariance. (Optional)
% S.Cc            = Control covariance. (Optional)
%
% [OUTPUTS]
% bf_outs.LF      = leadfield matrix (Nchans x 3 x M). In XYZ format (i.e. a
%                   dipole at [1 0 0], [0 1 0], and [0 0 1]).
% bf_outs.Weights = Weights for each source location (Nchans x M).
% bf_outs.T_stat  = T stat for each source location.
% bf_outs.Z_stat  = Z stat for each source location.

[path.mri, filename.mri] = fileparts(S.mri_file);
path.mri = [path.mri];
if ~exist([path.mri 'meshes.mat'])
try
    mri = ft_read_mri([path.mri filename.mri '.mri']);
catch
    try
        mri = ft_read_mri([path.mri filename.mri '.nii']);
    catch
        error('Check MRI file is  in .mri or .nii format or that FieldTrip is added')
    end
end


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

%%%%%%%%%%% Single Sphere %%%%%%%%%%%
if strcmpi(fwd_method,'single')
    for ii = 1:size(sourcepos,1)
        % compute the lead fields
        [LF(:,:,ii)] = Forward_OPM_single_sphere_XYZ_dir(sourcepos(ii,1), sourcepos(ii,2), sourcepos(ii,3),...
            S.sensor_info.pos,S.sensor_info.ors,origin);
        fprintf('Done %u/%u \n',ii,size(sourcepos,1))
    end
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end
    for n = 1:size(sourcepos,1)
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
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
%     figure
%     ft_plot_mesh(meshes(1),'facecolor',[.5 .5 .5],'facealpha',0.3,'edgecolor','none')
%     hold on
%     [sx1,sy1,sz1] = sphere;
%     axis equal
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
        fprintf('Done %u/%u \n',ii,size(sourcepos,1))
    end
    
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end

    for n = 1:size(sourcepos,1)
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
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
    if MFC % Account for mean field correction
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end
    for n = 1:size(sourcepos,1)
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end
    
elseif strcmpi(fwd_method,'shell')
    sensor_info.pos = S.sensor_info.pos;
    sensor_info.ors = S.sensor_info.ors;
    sensor_info.label = strsplit(num2str(1:size(sensor_info.pos)));
    
    [LF_single_shell,L_reshaped] = Forward_OPM_shell_FT([path.mri,filesep,filename.mri],sensor_info,sourcepos);
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
    if MFC % Account for mean field correction
       
        for n = 1:size(LF,3)
            LF(:,:,n) = S.M*LF(:,:,n);
        end
    end
    
    for n = size(sourcepos,1):-1:1
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan(:,:,n) = LF(:,:,n)*tan_ors';
    end
    
end

bf_outs.LF = LF;
bf_outs.ltan = ltan;
if do_bf   
    for n = 1:size(sourcepos,1)
        % project to two tangential fields
        R = sourcepos(n,:) - origin;
        [etheta,ephi] = calctangent_all_bf(R);
        tan_ors = [etheta;ephi];
        ltan = LF(:,:,n)*tan_ors';
        
        % get optimal combination for max SNR
        [v,d] = svd(ltan'*S.Cinv*ltan);
        [~,id] = min(diag(d));
        lopt = ltan*v(:,id); % turn to nAm amplitude
        tan_or_opt = [tan_ors'*(v(:,id))]';
        
        bf_outs.lopt(:,:,n) = lopt;
        bf_outs.or_opt(n,:) = tan_or_opt;
        W1(:,n) = (S.Cinv*lopt/(lopt'*S.Cinv*lopt));
        Z{n} = [W1(:,n)'*S.C*W1(:,n)]./[W1(:,n)'*eye(size(W1(:,n),1))*W1(:,n)];
        if isfield(S,'Ca')
            Qa = W1(:,n)'*S.Ca*W1(:,n);
            Qc = W1(:,n)'*S.Cc*W1(:,n);
            T1{n} = (Qa-Qc)./(2*Qc);
        end
        fprintf('Done %u/%u \n',n,size(sourcepos,1))
    end

    bf_outs.Weights = W1;
    bf_outs.Z_stat = Z;
    if isfield(S,'Ca')
        bf_outs.T_stat = T1;
%         max_or_min = questdlg('Do you want to see the maximum or minimum T-stat?',...
%             'Max or Min','Max','Min','Cancel','Min');
        
        if strcmpi(max_or_min,'Max')
            cax_val = [0.5*max(cell2mat(T1)) max(cell2mat(T1))];
            [~, l] = max(cell2mat(T1));
        elseif strcmpi(max_or_min,'Min')
            cax_val = [min(cell2mat(T1)), 0.5.*min(cell2mat(T1))];
            [~, l] = min(cell2mat(T1));
        end
        try
            fg = figure;
            fg.Name = fwd_method;
            subplot(2,3,[1 2 4 5])
            ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
            hold on
            scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,cell2mat(T1),'filled')
            colormap hot;caxis([cax_val]);cb = colorbar;cb.Label.String = 'Tstat';
            axis equal
            view([120 20])
            fig = gcf;
            fig.Color = 'w';
            ax = gca;
            ax.FontSize = 14;
            
            subplot(2,3,[3 6])
            % Find voxel with min/max value of Tstat
            dip_loc = sourcepos(l,:);
            dip_or = bf_outs.or_opt(l,:);
            ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
            hold on
            plot3(dip_loc(1),dip_loc(2),dip_loc(3),'.','markersize',20)
            quiver3(dip_loc(1),dip_loc(2),dip_loc(3),bf_outs.or_opt(l,1)./100,bf_outs.or_opt(l,2)./100,bf_outs.or_opt(l,3)./100)
            bf_outs.dip_loc = dip_loc;
            bf_outs.dip_or = dip_or;
        catch
            disp('No min or max T stat preference selected')
        end
    end
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