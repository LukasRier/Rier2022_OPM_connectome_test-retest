function [lead_fields,UnitMDip] = Forward_OPM_single_sphere(x,y,z,slots_all,normal_unit,sphere_origin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mean_sphere_origin = sphere_origin;
%% number of channels
Nch = size(slots_all,1);
%% voxel position in CTF and in metres
vox_pos_ctf = [x y z];
%% OPM positions in CTF and in metres
slots_all = slots_all;
%%convert voxel and OPM locations to sphere regerence
vox_pos_local = vox_pos_ctf - mean_sphere_origin;
coil_position_local = slots_all - repmat(mean_sphere_origin,size(slots_all,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% number of angles and set up count matrix
N = 90;
count = 0;
for aa = 1:N
    for bb = 1:Nch
        count = count + 1;
        countmatrix(count) = aa;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compute unit dipole vectors
angles = (1:(180/N):180)./180.*pi;
[etheta, ephi] = calctangent(vox_pos_local);
u=cos(angles);
v=sin(angles);
%%dipole in tangential plane and cartesian co-ordinates
for ix=1:3
    MDip(:,ix)=u.*etheta(ix)+v.*ephi(ix);
end
%%Normalise to obtain dipole moment in for lead field calculation
MagDip = repmat(sqrt(dot(MDip',MDip')),3,1)';
UnitMDip = MDip./MagDip;
%%Convert to Am from nAm
UnitMDip = UnitMDip.*1e-9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE LEAD FIELD
%%Calculate voxel position w.r.t. local sphere origin (NB sensor -1 since sphere origins indexed like this)
voxpos_local_coils = repmat(vox_pos_local,size(slots_all,1),1);
%%vector displacements between the source and sensors
a = coil_position_local - voxpos_local_coils;
%%distances between source and sensors
am = sqrt(dot(a',a'));
%%distances between origin and sensors
rm = sqrt(dot(coil_position_local',coil_position_local'));
%%calculate f
ff = am.*(am.*rm+rm.^2-dot(voxpos_local_coils',coil_position_local'));
%%component parts of delf
K = ((am.^2)./rm + dot(a',coil_position_local')./am) + 2.*am + 2.*rm;
L = am + 2.*rm + (dot(a',coil_position_local')./am);
%%delf
delf = repmat(K,3,1).*coil_position_local'-repmat(L,3,1).*voxpos_local_coils';
%%Forward solution will hold the solution for B for each sensor
ro_matrix = repmat(voxpos_local_coils,N,1);
Q_matrix =  UnitMDip(countmatrix,:);
ff_matrix = repmat(ff',N,3);
delf_matrix = repmat(delf',N,1);
r_matrix = repmat(coil_position_local,N,1);
B_matrix = (10^-7./(ff_matrix.^2)).*(ff_matrix.*cross(Q_matrix,ro_matrix)-...
    repmat(dot(cross(Q_matrix,ro_matrix)',r_matrix'),3,1)'.*delf_matrix);
B_matrix2 = reshape(B_matrix,Nch,N,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%plot stuff
% figure
% plot3(coil_position_local(:,1),coil_position_local(:,2),coil_position_local(:,3),'r+')
% hold on
% plot3(vox_pos_local(1),vox_pos_local(2),vox_pos_local(3),'bo')
% plot3(0,0,0,'ko')
% axis equal
% angle_plot = 1;
% quiver3(vox_pos_local(1),vox_pos_local(2),vox_pos_local(3),1e7.*UnitMDip(angle_plot,1),1e7.*UnitMDip(angle_plot,2),1e7.*UnitMDip(angle_plot,3),'r','linewidth',2)
% quiver3(repmat(vox_pos_local(1),90,1),repmat(vox_pos_local(2),90,1),repmat(vox_pos_local(3),90,1),1e7.*UnitMDip(:,1),1e7.*UnitMDip(:,2),1e7.*UnitMDip(:,3))
% B_matrix2_plot = squeeze(B_matrix2(:,angle_plot,:));
% quiver3(coil_position_local(:,1),coil_position_local(:,2),coil_position_local(:,3),B_matrix2_plot(:,1),B_matrix2_plot(:,2),B_matrix2_plot(:,3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lead_fields = zeros(Nch,N);
for sensor = 1:Nch
    B_thissensor = squeeze(B_matrix2(sensor,:,:))';
    ors_mat = repmat(normal_unit(sensor,:)',1,N);
    lead_fields(sensor,:) = dot(B_thissensor,ors_mat);
end
end

function [tanu, tanv] = calctangent(RDip)
%% Same as usual calctangent but put here for ease

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