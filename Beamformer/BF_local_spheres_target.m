function [OUTS] = BF_local_spheres_target(S)
% Create lead fields and determine beamformer weights. All vars stored in S
% and output in OUTS.
%
% INPUT
% S.dip_loc        = Target dipole location
% S.pos            = Sensor positions (in m)
% S.ors            = Sensor orientations
% S.ls_Origins     = Local sphere origins (use create_local_spheres)
% S.tan_method     = Either 'mean' for mean origin of spheres or 'norm' for
%                 closest surface normal to the dipole location
% S.tan_ref        = Values for either the mean origin or the surface normal
% S.C              = Covariance matrix
% S.Cinv           = Inverse of covariance
% S.Noise_Cr       = Noise covariance
%
% OUTPUT
% OUTS.lead_fields = Lead fields for all angles;
% OUTS.Zang        = Z value for all angles;
% OUTS.Weights     = Weights;
% OUTS.leads_vec   = Lead fields for best angle;
% OUTS.best_or     = Best orientation;
% OUTS.best_ang    = Best angle;

% number of channels
slots_all = S.pos;
normal_unit = S.ors;
Nch = size(slots_all,1);

% voxel position in metres
vox_pos_ctf = S.dip_loc;

% convert voxel and OPM locations to sphere regerence
sphere_origins = S.ls_origins;
vox_pos_local = vox_pos_ctf - sphere_origins;
coil_position_local = slots_all - sphere_origins;

% number of angles and set up count matrix
N = 90;
count = 0;
for aa = 1:N
    for bb = 1:Nch
        count = count + 1;
        countmatrix(count) = aa;
    end
end

% Compute unit dipole vectors
if strcmpi(S.tan_method,'mean')
    vox_pos_meansp = vox_pos_ctf - S.tan_ref;
    angles = (1:(180/N):180)./180.*pi;
    [etheta, ephi] = calctangent_ls(vox_pos_meansp);
    u=cos(angles);
    v=sin(angles);
    
    % dipole in tangential plane and cartesian co-ordinates
    for ix=1:3
        MDip(:,ix)=u.*etheta(ix)+v.*ephi(ix);
    end
    
    % Normalise to obtain dipole moment in for lead field calculation
    MagDip = repmat(sqrt(dot(MDip',MDip')),3,1)';
    UnitMDip = MDip./MagDip;
    
    % Convert to Am from nAm
    UnitMDip = UnitMDip.*1e-9;
end

if strcmpi(S.tan_method,'norm')
    angles = (1:(180/N):180)./180.*pi;
    [etheta,ephi] = calctangent_ls(S.tan_ref);
    u=cos(angles);
    v=sin(angles);
    for ix=1:3
        MDip(:,ix)=u.*etheta(ix)+v.*ephi(ix);
    end
    MagDip = repmat(sqrt(dot(MDip',MDip')),3,1)';
    UnitMDip = MDip./MagDip;
    
    %Convert to Am from nAm
    UnitMDip = UnitMDip.*1e-9;
    %     figure
    %     quiver3(1,1,1,n(1),n(2),n(3),'b')
    %     axis equal
    %     hold on
    %     tst = ones(90,3);
    %     quiver3(tst(:,1),tst(:,2),tst(:,3),UnitMDip(:,1),UnitMDip(:,2),UnitMDip(:,3),'r')
    
end

% CALCULATE THE LEAD FIELD
% Calculate voxel position w.r.t. local sphere origin (NB sensor -1 since sphere origins indexed like this)
% voxpos_local_coils = repmat(vox_pos_local,size(slots_all,1),1);
voxpos_local_coils = vox_pos_local;

% vector displacements between the source and sensors
a = coil_position_local - voxpos_local_coils;

% distances between source and sensors
am = sqrt(dot(a',a'));

% distances between origin and sensors
rm = sqrt(dot(coil_position_local',coil_position_local'));

% calculate f
ff = am.*(am.*rm+rm.^2-dot(voxpos_local_coils',coil_position_local'));

% component parts of delf
K = ((am.^2)./rm + dot(a',coil_position_local')./am) + 2.*am + 2.*rm;
L = am + 2.*rm + (dot(a',coil_position_local')./am);

% delf
delf = repmat(K,3,1).*coil_position_local'-repmat(L,3,1).*voxpos_local_coils';

% Forward solution will hold the solution for B for each sensor
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

% Compute beamformer weights
C = S.C;
Cinv = S.Cinv;
Noise_Cr = S.Noise_Cr;
for angle = 1:size(lead_fields,2)
    Aa(:,angle) = (Cinv*lead_fields(:,angle))/(lead_fields(:,angle)'*...
        Cinv*lead_fields(:,angle));
    Zang_a(angle) = (Aa(:,angle)'*C*Aa(:,angle))/(Aa(:,angle)'*Noise_Cr*Aa(:,angle));
    X(angle) = (Aa(:,angle)'*C*Aa(:,angle))/(Aa(:,angle)'*(1e-30*Noise_Cr)*Aa(:,angle));
end
I = find(Zang_a==max(Zang_a));
W1 = Aa(:,I);
leads_vec = lead_fields(:,I);
Best_OR = UnitMDip(I,:);
best_ang = (180/size(lead_fields,2))*I;

% Outputs
OUTS.lead_fields = lead_fields;
OUTS.Zang = Zang_a;
OUTS.Weights = W1;
OUTS.leads_vec = leads_vec;
OUTS.best_or = Best_OR;
OUTS.best_ang = best_ang;

    function [tanu, tanv] = calctangent_ls(RDip)
        % Same as calctangent but put here to save adding paths
        
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
end

