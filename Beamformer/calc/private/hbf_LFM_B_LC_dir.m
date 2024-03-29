%  HBF_LFM_B_LC_DIR builds magnetic lead field matrix based on directed 
%    current dipoles.
% 
%  LFM=HBF_LFM_B_LC_DIR(meshes,coils,TB,spos,sdir)
%    meshes: BEM geometry, cell array of hbf structs
%    coils:  coil description, hbf struct
%    TB:     TB matrix built with the hbf BEM solver
%    spos:   source positions, [M x 3]
%    sdir:   source orientations (unit-length), [M x 3]
% 
%    LFM:   lead field matrix, [Number of coils (field points) x M]
%            [l_1 ... l_M]
% 
%  You can also compute magnetic field due to any set of directed dipoles by 
%  giving the dipole moments (with amplitude) in the 'sdir' argument.
% 
%  v160229 Matti Stenroos
%