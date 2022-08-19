% HBF_SOLIDANGLES calculates solid angles spanned by triangles of a mesh in 
%    one point
% 
% omega=HBF_SOLIDANGLES(mesh,fieldpoint) 
% omega=HBF_SOLIDANGLES(elements,vertices,fieldpoint) 
%    mesh:   hbf mesh struct
%    elements: triangle description, [N(triangles) x 3]
%    vertices:   mesh vertices, [N(vertices) x 3] 
%    fieldpoint: point, in which the angles are evaluated, [1 x 3]
% 
%  Eq. 8 of vanOosterom & Strackee, IEEE TBME 1983
% 
%  v160229 Matti Stenroos
%