function LFMm=hbf_LFM_B_LC(bmeshes,coils,TB,spos,sdir)
% HBF_LFM_B_LC builds magnetic lead field matrix based
%
% LFM=HBF_LFM_B_LC(meshes,coils,TB,spos,sdir)
% LFM=HBF_LFM_B_LC(meshes,coils,TB,spos)
%   meshes: BEM geometry, cell array of hbf structs
%   coils:  coil description, hbf struct
%   TB:     TB matrix built with the hbf BEM solver
%   spos:   source positions, [M x 3]
%   sdir:   source orientations (unit-length), [M x 3]; optional
%
%   LFM:   lead field matrix, [Number of coils (field points) x M]
%           [l_1 ... l_M]
%           or, if 'sdir' omitted, [Number of coils x 3M]
%           [l_1x l_1y l1_z ... l_Mx l_My l_Mz]
%
% You can also compute magnetic field due to any set of directed dipoles by 
% giving the dipole moments (with amplitude) in the 'sdir' argument. If you
% omit 'sdir', LFM will be computed using xyz unit dipole triplets

%
% v160404 Matti Stenroos
if nargin==5,
    LFMm=hbf_LFM_B_LC_dir(bmeshes,coils,TB,spos,sdir);
else
    LFMm=hbf_LFM_B_LC_xyz(bmeshes,coils,TB,spos);
end