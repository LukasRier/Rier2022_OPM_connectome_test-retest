function [Sphere_orig_new,Sphere_rad_new] = create_local_spheres(mesh,sens_info)
% Create local spheres. mesh is the surface mesh to fit local spheres to
% (set as mesh.vertices and mesh.faces). sens_info is the sensor
% information using sens_info.pos and sens_info.ors.
addpath('R:\Matlab_files\tools');
vn = patchnormals(mesh);
bm = mesh.vertices;
options = optimset('MaxFunEvals',4800,'MaxIter',4800,'TolFun',1e-4,'TolX',1e-4);

ccc = jet(size(sens_info.pos,1));
for n = 1:size(sens_info.pos,1)
    disp(['Running Sensor ' num2str(n) ' of ' num2str(size(sens_info.pos,1))])
    % Finds the nearest 10% of points to each sensor, then fit a sphere to them
    % to estmate curvature
    o_idx = sens_info.ors(n,:);
    r_idx = sens_info.pos(n,:);
    sens_brain_dist = sqrt(sum([bm - o_idx].^2,2));
    [~, ord_idx] = sort(sens_brain_dist);
    nrst_pts = bm(ord_idx(1:round(0.1*size(bm,1))),:); % Find the nearest 5% of points
    [Sphere_orig(n,:),Sphere_rad(n,:)] = sphereFit(nrst_pts);
    
    % Method from Huang et al 1999, Head model comparison for magnetoencephalography
    % Overlapping spheres without BEM
    % The spheres fit in the previous section are used as the estimates before
    % getting (hopefully) improved
    X0 = [Sphere_orig(n,:), Sphere_rad(n,:)];
    F = costfn_ls(X0);
    fun = @(X0)costfn_ls(X0);
    X = fminsearch(fun,X0,options);
    Sphere_rad_new(n,:) = X(4);
    Sphere_orig_new(n,:) = X(1:3);
    
end

    function [F] = costfn_ls(X0)
        C0 = X0(1:3);
        R0 = X0(4);
        for NN = 1:size(bm,1)
            r3 = [[r_idx - bm(NN,:)]./[(sqrt(sum((r_idx - bm(NN,:)).^2))).^3]];
            rC = [[bm(NN,:) - C0]./[sqrt(sum((bm(NN,:) - C0).^2))]];
            rRC = [[r_idx - C0]./[sqrt(sum((r_idx - R0.*(rC) - C0).^2)).^3]];
            tst(NN,:) = sqrt(sum([dot(o_idx,vn(NN,:))*r3 - dot(o_idx,rC)*rRC].^2)).^2;
        end
        F = sum(tst);
        % F = sqrt(sum(EQ10.^2));
        
    end
end

