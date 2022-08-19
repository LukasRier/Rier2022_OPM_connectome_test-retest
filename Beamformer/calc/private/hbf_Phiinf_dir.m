%  HBF_PHIINF_DIR  computes potential due to a current dipole in infinite
%    homogeneous conductor that has unit conductivity
% 
%  phiinf=HBF_PHIINF_DIR(fp,spos,sdir)
%    fp = field points (where the field is computed), [N x 3]
%    spos  = source positions, [M x 3]
%    spdir = dipole moments (or, if unit-norm, dipole orientations), [M x 3]
% 
%    phiinf = resulting potential for each field point and dipole, [N x M]
% 
%  v160229 Matti Stenroos
%