%  HBF_DB_Linear makes a DB matrix between a mesh and field points.
%    The matrix is made for an arbitrary field component and linear basis.
% 
%  DB=HBF_DB_LINEAR(fieldpoints,fielddir,sourcemesh)
%    fieldpoints: a set of points, where DB is computed, [N x 3]
%    fielddir: field directions (unit length), [N x 3]
%    sourcemesh: mesh, where the double layer is spanned, hbf struct
% 
%    DB: DB matrix, [N(field points) x N(source vertices)]
% 
%  v160229 Matti Stenroos
%