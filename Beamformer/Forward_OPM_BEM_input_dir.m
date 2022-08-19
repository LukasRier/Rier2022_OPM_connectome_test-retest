function [LF_matrix,UnitMDip] = Forward_OPM_BEM_input_dir(dip_loc,sourceor,bmeshes,coils,TBvol)
% dip_loc = dipole location in metres
% bmeshes = meshes in order of innerskull, outerskull, scalp
% coils = sensor info as coils.p = sensor positions, coils.n = sensor
% orientations (in metres)
% nrst_norm = nearest normal of the brain mesh to define tangent
% p = mfilename('fullpath');
% [p1,p2] = fileparts(p);
% addpath([p1 '\calc'])
% addpath R:\Matlab_files\tools\

% Find tangential components
% tic
mesh.vertices = bmeshes{1,1}.p;
mesh.faces = bmeshes{1,1}.e;
% vn = patchnormals(mesh);
% source2mesh_dist = pdist2(dip_loc,mesh.vertices);
% [~,s2m_idx] = min(source2mesh_dist);
% nrst_norm = vn(s2m_idx,:);
% % toc
% % tic
% N = 3;
% angles = (1:(180/N):180)./180.*pi;
% [etheta,ephi] = calctangent_bem(nrst_norm);
% u=cos(angles);
% v=sin(angles);
% for ix=1:3
%     MDip(:,ix)=u.*etheta(ix)+v.*ephi(ix);
% end
MDip = sourceor;
MagDip = repmat(sqrt(dot(MDip',MDip')),3,1)';
UnitMDip = MDip./MagDip;

% Convert to Am from nAm
UnitMDip = UnitMDip.*1e-9;
% toc
% Forward solutions for directed dipoles
% tic
for n = 1:size(UnitMDip,1)
    LF_matrix(:,n) = hbf_LFM_B_LC(bmeshes,coils,TBvol,dip_loc,UnitMDip(n,:));
end
% toc
end

function [tanu, tanv] = calctangent_bem(RDip)
% Same as usual but put here for ease

x=RDip(1);
y=RDip(2);
z=RDip(3);
r=sqrt(x*x+y*y+z*z);

if (x==0) && (y==0)
    tanu(1)=1.0; tanu(2)=0; tanu(3)=0;
    tanv(1)=0; tanv(2)=1.0; tanv(3)=0;
else
    RZXY= -(r-z)*x*y;
    X2Y2= 1/(x*x+y*y);
    
    tanu(1)= (z*x*x + r*y*y) * X2Y2/r;
    tanu(2)= RZXY * X2Y2/r;
    tanu(3)= -x/r;
    
    tanv(1)= RZXY * X2Y2/r;
    tanv(2)= (z*y*y + r*x*x) * X2Y2/r;
    tanv(3)= -y/r;
end
end
