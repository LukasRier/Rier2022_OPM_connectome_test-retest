function [TBvol] = create_BEM(bmeshes,coils,ci)
% bmeshes = meshes in order of innerskull, outerskull, scalp
% coils = sensor info as coils.p = sensor positions, coils.n = sensor
% orientations (in metres)
% ci = conductivity inside each surface --- remember the order!
p = mfilename('fullpath');
[p1,p2] = fileparts(p);
addpath([p1 '\calc'])

co = [ci(2:3) 0];

%   BEM geometry matrices
D = hbf_BEMOperatorsPhi_LC(bmeshes);
disp('Done 1/2')
DB = hbf_BEMOperatorsB_Linear(bmeshes,coils);
disp('Done 2/2')

%   (full) transfer matrices
Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);
TBvol = hbf_TM_Bvol_Linear(DB,Tphi_full,ci,co);

end

