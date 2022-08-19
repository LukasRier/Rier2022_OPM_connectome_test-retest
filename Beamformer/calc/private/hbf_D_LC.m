%  HBF_DI_LC makes a LC double-layer matrix between two meshes 
% 
%  D=HBF_D_LC(fieldmesh,sourcemesh,verbose)
%    fieldmesh: mesh, where double-layer potential is evaluated, hbf struct
%    sourcemesh: mesh, where the double layer is spanned, hbf struct
%    verbose (optional, default 0): give 1 for intermediate output
% 
%    D: Double-layer matrix, [N(field vertices) x N(source vertices)]
% 
%  v160229 Matti Stenroos
%