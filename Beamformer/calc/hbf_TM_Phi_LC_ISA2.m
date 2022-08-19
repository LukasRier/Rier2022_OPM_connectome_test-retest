function [T_ISA,startinds,endinds]=hbf_TM_Phi_LC_ISA2(D,ci,co,ISAsurf)
% HBF_TM_PHI_LC_ISA2 makes a BEM transfer matrix for potential using LC 
%   BEM and isolated source approach
%
% [T,startinds,endinds]=HBF_TM_PHI_LC_ISA2(D,ci,co,ISAsurf)
%   D:  D matrix, cell [N(meshes) x N(meshes)]
%   ci: conductivity inside each boundary surface, [N(meshes) x 1]
%   co: conductivity outside each boundary surface, [N(meshes) x 1]
%   ISAsurf:  index to the mesh, on which the isolation is performed. 
%
%   T: BEM transfer matrix, N(vertices in the whole model)^2
%   startinds, endinds: indices that point to rows/columns of T that
%       correspond to each BEM boundary surface;
%       inds{1}=startinds(1):endinds(1)
% 
%   In EEG/MEG application, the isolation is typically done on the inner skull
%   boundary (mesh index 1 in three-shell model, 2 or 3 in a four-compartment
%   model). All sources must be inside the isolation boundary.
%
%   The zero level of potential is set on the isolation surface.
%
% v160303 Matti Stenroos
Nsurf=size(D,1);
zerolevel=ISAsurf;
Nop_ISA=size(D{ISAsurf,ISAsurf},1);
defvals=DefValuesD(D,zerolevel);

[T,startinds,endinds]=hbf_TM_Phi_LC(D,ci,co,zerolevel);

coISA=co(1:ISAsurf);
coISA(end)=0;
T0_ISA=hbf_TM_Phi_LC(D(1:ISAsurf,1:ISAsurf),ci(1:ISAsurf),coISA,zerolevel);
D0=zeros(endinds(end),Nop_ISA);
for I=1:Nsurf,
    D0(startinds(I):endinds(I),:)=co(ISAsurf)*(D{I,ISAsurf}+defvals(I,ISAsurf));
end
D0(startinds(ISAsurf):endinds(ISAsurf),:)=D0(startinds(ISAsurf):endinds(ISAsurf),:)-co(ISAsurf)*.5*eye(Nop_ISA);

T_ISA=T*D0*T0_ISA(startinds(ISAsurf):endinds(ISAsurf),:);
T_ISA(startinds(1):endinds(ISAsurf),:)=T_ISA(startinds(1):endinds(ISAsurf),:)+T0_ISA;

function defvals=DefValuesD(D,defsurf)
% function defvals=DefValuesD(D,defsurf);
% Calculates deflation values that set the zero level to surface "defsurf".
Nsurf=size(D,1);
N_defsurf=size(D{defsurf,defsurf},1);
defvals=zeros(size(D));
for S=1:Nsurf,
    defvals(S,defsurf)=1/N_defsurf;
end
