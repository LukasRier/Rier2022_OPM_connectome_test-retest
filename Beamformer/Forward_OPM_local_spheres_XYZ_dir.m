function [lead_fields,UnitMDip] = Forward_OPM_local_spheres_XYZ_dir(x,y,z,slots_all,normal_unit,sphere_origins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% number of channels
Nch = size(slots_all,1);
%% voxel position in CTF and in metres
vox_pos_ctf = [x y z];
%% OPM positions in CTF and in metres
slots_all = slots_all;
%%convert voxel and OPM locations to sphere regerence
vox_pos_local = vox_pos_ctf - sphere_origins;
coil_position_local = slots_all - sphere_origins;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% number of angles and set up count matrix
N = 3;
count = 0;
for aa = 1:N
    for bb = 1:Nch
        count = count + 1;
        countmatrix(count) = aa;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Compute unit dipole vectors
%     vox_pos_meansp = vox_pos_ctf - tan_ref;
%     angles = (1:(180/N):180)./180.*pi;
%
%     [etheta, ephi] = calctangent(vox_pos_meansp);
%
%     u=cos(angles);
%     v=sin(angles);
%     %%dipole in tangential plane and cartesian co-ordinates
%     for ix=1:3
%         MDip(:,ix)=u.*etheta(ix)+v.*ephi(ix);
%     end
MDip = [1 0 0; 0 1 0; 0 0 1];
%%Normalise to obtain dipole moment in for lead field calculation

MagDip = repmat(sqrt(dot(MDip',MDip')),3,1)';
UnitMDip = MDip./MagDip;
%%Convert to Am from nAm
UnitMDip = UnitMDip.*1e-9;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE THE LEAD FIELD
%%Calculate voxel position w.r.t. local sphere origin (NB sensor -1 since sphere origins indexed like this)  
% voxpos_local_coils = repmat(vox_pos_local,size(slots_all,1),1);
voxpos_local_coils = vox_pos_local;
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
% figure
% plot3(vox_pos_local(:,1),vox_pos_local(:,2),vox_pos_local(:,3),'.')
% axis equal
% hold all
% for n = 1:length(vox_pos_local)
%     disp(n)
%     for NN = 1:N
% %         disp(NN)
%         quiver3(vox_pos_local(n,1),vox_pos_local(n,2),vox_pos_local(n,3),...
%             UnitMDip(NN,1).*1e7,UnitMDip(NN,2).*1e7,UnitMDip(NN,3).*1e7)
%     end
% end
lead_fields = zeros(Nch,N);
for sensor = 1:Nch
    B_thissensor = squeeze(B_matrix2(sensor,:,:))';
    ors_mat = repmat(normal_unit(sensor,:)',1,N);
    lead_fields(sensor,:) = dot(B_thissensor,ors_mat);
end

