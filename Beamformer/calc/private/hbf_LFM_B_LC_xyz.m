%  HBF_LFM_B_LC_XYZ builds magnetic lead field matrix based on xyz-oriented
%    unit-current dipoles
% 
%  LFMm=HBF_LFM_B_LC_XYZ(meshes,coils,TB,spos)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    TB:     TB matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
% 
%    LFMm:   lead field matrix, [Number of coils x 3M]
%        [l_1x l_1y l1_z ... l_Mx l_My l_Mz]
% 
%  v160229 Matti Stenroos
%