function Conn_analysis_MJB(sub,ses)
%
%clear all
close all
clc
restoredefaultpath
%sub = '010';
%ses = '001';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\ppzmjb2\OneDrive - The University of Nottingham\Matlab_files\fieldtrip-20190212')
addpath('C:\Users\ppzmjb2\The University of Nottingham\OPM Big Data PC - Elena\dsmovie_allfiles\Matt')

ft_defaults
datadir = 'C:\Users\ppzmjb2\The University of Nottingham\OPM Big Data PC - Elena\dsmovie_allfiles\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = ['sub-',sub,'_','ses-',ses];
path.main = [datadir,'sub-',sub,'\'];
path.data = [path.main,'ses-',ses,'\'];
path.mri = [path.main,'mri\'];
files.mri = ['sub-',sub,'.nii'];
files.meshes = 'meshes.mat';                                                                                                            
global files
cd(path.data)
load([path.data,filename,'_meg.mat'])
basefilename = [path.main,'ses-',ses,'\Matt\'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check sensors are co-registered to MRI
channels = importdata([path.data,filename,'_channels_coreg.tsv']);
helmet.sens_pos = channels.data(:,1:3);
helmet.sens_ors_X = channels.data(1:3:end,4:6);
helmet.sens_ors_Y = channels.data(2:3:end,4:6);
helmet.sens_ors_Z = channels.data(3:3:end,4:6);
mri = ft_read_mri([path.mri files.mri]);
cfg = [];
cfg.output    = {'scalp'; 'skull' ; 'brain'};
cd(path.mri)
if exist('segmentedmri.mat')
    load('segmentedmri.mat')
else
    segmentedmri  = ft_volumesegment(cfg, mri);
    BETbrain = ft_read_mri([path.mri 'sub-' sub '_brain.nii'])
    BETlog = BETbrain.anatomy > 1000;
    segmentedmri.brain = BETlog;
    save segmentedmri.mat segmentedmri
end

cfg = [];
cfg.tissue = {'scalp'; 'brain'};
cfg.numvertices = [5000];
mesh2 = ft_prepare_mesh(cfg,segmentedmri);
mesh1 = ft_convert_units(mesh2,'m');
for n = 1:size(mesh1,2)
    meshes(n).pnt = mesh1(n).pos;
    meshes(n).tri = mesh1(n).tri;
    meshes(n).unit = mesh1(n).unit;
    meshes(n).name = cfg.tissue{n};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the source locations
load([path.mri 'voxlox.mat']);
XX = mri.anatomy;
%% this is a hack!
if strmatch(sub,'009')
    voxlox(3,71) = 158;
end
if strmatch(sub,'007')
    voxlox(3,73) = 125;
end

if strmatch(sub,'010')
    voxlox(3,71) = 170;
end

for n = 1:78
    XX(voxlox(1,n),voxlox(2,n),voxlox(3,n)) = 32767;
    sourcepos1 = ft_warp_apply(mri.transform,[voxlox(1,n) voxlox(2,n) voxlox(3,n)]);
    sourcepos(:,n) = sourcepos1'/1000; % convert to metres
    sourcepos(:,n)'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the lead fields in cartesian axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'C:\Users\ppzmjb2\OneDrive - The University of Nottingham\Matlab_files\Beamformer'
S.mri_file = [path.mri files.mri];
S.sensor_info.pos = [helmet.sens_pos];
S.sensor_info.ors = zeros(size(S.sensor_info.pos));
S.sensor_info.ors(1:3:end,:) = helmet.sens_ors_X;
S.sensor_info.ors(2:3:end,:) = helmet.sens_ors_Y;
S.sensor_info.ors(3:3:end,:) = helmet.sens_ors_Z;
dip_loc = sourcepos';                                            %Dipole locations
[bf_outs_shell] = run_beamformer2('shell',dip_loc,S,0); 
lead_fields_shell_xyz = bf_outs_shell.LF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% convert orientation of sources to polar
X = meshes.pnt; Origin = mean(X,1);
Ndips = length(sourcepos);
for n = 1:Ndips;
    thispos = sourcepos(:,n);
    [phi,theta1,r] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,thispos(3) - Origin(3));
    theta = pi/2 - theta1;
    Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
    Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
    Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
    Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make a plot of the geometry...
figure(1);
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
scatter3(sourcepos(1,:),sourcepos(2,:),sourcepos(3,:),'ro','linewidth',3)
view([130,30])
fig = gcf;
fig.Color = [1,1,1];
plot3(helmet.sens_pos(:,1),helmet.sens_pos(:,2),helmet.sens_pos(:,3),'o')
quiver3(helmet.sens_pos(1:3:end,1),helmet.sens_pos(1:3:end,2),helmet.sens_pos(1:3:end,3),...
    S.sensor_info.ors(1:3:end,1),S.sensor_info.ors(1:3:end,2),S.sensor_info.ors(1:3:end,3),'r','linewidth',2)
quiver3(helmet.sens_pos(1:3:end,1),helmet.sens_pos(1:3:end,2),helmet.sens_pos(1:3:end,3),...
    S.sensor_info.ors(2:3:end,1),S.sensor_info.ors(2:3:end,2),S.sensor_info.ors(2:3:end,3),'g','linewidth',2)
quiver3(helmet.sens_pos(1:3:end,1),helmet.sens_pos(1:3:end,2),helmet.sens_pos(1:3:end,3),...
    S.sensor_info.ors(3:3:end,1),S.sensor_info.ors(3:3:end,2),S.sensor_info.ors(3:3:end,3),'b','linewidth',2)
plot3(Origin(1),Origin(2),Origin(3),'bo','linewidth',4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% take a random lead field and plot it...
figure(2);
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
ft_plot_topo3d(double(helmet.sens_pos(1:3:end,:)),Lead_fields(3:3:end,2,16))
alpha(gca,0.5)
plot3(helmet.sens_pos(1:3:end,1),helmet.sens_pos(1:3:end,2),helmet.sens_pos(1:3:end,3),'go','linewidth',3)
scatter3(sourcepos(1,16),sourcepos(2,16),sourcepos(3,16),'r','linewidth',4)
quiver3(helmet.sens_pos(1:3:end,1),helmet.sens_pos(1:3:end,2),helmet.sens_pos(1:3:end,3),...
    S.sensor_info.ors(1:3:end,1).*Lead_fields(1:3:end,2,16) + S.sensor_info.ors(2:3:end,1).*Lead_fields(2:3:end,2,16) + S.sensor_info.ors(3:3:end,1).*Lead_fields(3:3:end,2,16),...
    S.sensor_info.ors(1:3:end,2).*Lead_fields(1:3:end,2,16) + S.sensor_info.ors(2:3:end,2).*Lead_fields(2:3:end,2,16) + S.sensor_info.ors(3:3:end,2).*Lead_fields(3:3:end,2,16),...
    S.sensor_info.ors(1:3:end,3).*Lead_fields(1:3:end,2,16) + S.sensor_info.ors(2:3:end,3).*Lead_fields(2:3:end,2,16) + S.sensor_info.ors(3:3:end,3).*Lead_fields(3:3:end,2,16),'r','linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% conditino the data and use the beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% segment the data to when the movie was playing
tr = diff(triggers(:,1)>2);
start_sample1 = find(tr == 1);
start_sample = start_sample1(1);
end_sample = find(tr == -1);
duration = floor((end_sample - start_sample)/fs);
data1 = data(start_sample+1:start_sample+duration*fs,:);
%% filter the OPM data 13-30 Hz
hp = 13;
lp = 30;
[b,a] = butter(4,2*[hp lp]/fs);
data_f = [filtfilt(b,a,data1)]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strmatch(sub,'010') & strmatch(ses,'001')
    epoch_length = 5;
    time_epoch = linspace(0, epoch_length, epoch_length*fs);
%     for n = 1:duration/epoch_length
%         figure(567);
%         clf
%         dat_epoch = data_f(:, (n-1)*(epoch_length*fs)+1:(n*epoch_length*fs));
%         imagesc(time_epoch,1:size(dat_epoch,1),dat_epoch);shading flat
%         colorbar
%         ylabel('channel')
%         xlabel('time (s)')
%         caxis([-500 500])
%         title(sprintf('Epoch %d',n));
%         pause(1)
%     end
    %% remove bad trials...
    data_f_mat = reshape(data_f,[size(data_f,1),epoch_length*fs,duration/epoch_length]);
    badtrials = [71 79 82 83 91]
    trials1 = ones(1,duration/epoch_length); trials1(badtrials) = 0;
    trials_to_keep = find(trials1 == 1);
    data_f_mat2 = data_f_mat(:,:,trials_to_keep);
    clear data_f
    clear duration
    duration = epoch_length*length(trials_to_keep)
    data_f = reshape(data_f_mat2,size(data_f_mat2,1),duration*fs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
C = cov(data_f');
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);
for n = 1:length(sourcepos)
     this_L = Lead_fields(:,:,n);
     W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
     iPower_v = this_L'*Cr_inv*this_L;
     [v,d] = svd(iPower_v);
     [~,id] = min(diag(d));
     lopt = this_L*v(:,id); % turn to nAm amplitude
     w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
     VE(:,n) = w*data_f;
end
if exist(basefilename);
    cd(basefilename)
else
    mkdir(basefilename)
    cd(basefilename)
end
save VE.mat VE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% down_f = 10;
% Nlocs = 78;
% VE_orth = symmetric_orthogonalise(VE, 1);
% for reg = 1:78
%     H_VE(:,reg) = abs(hilbert(VE_orth(:,reg)));
%     H_VE_d(:,reg) = mean(reshape(H_VE(:,reg),fs/down_f,duration*down_f,1));
% end
% AEC_b = corrcoef(H_VE_d);
% figure(1)
% subplot(121)
% imagesc(AEC_b);colorbar;
% subplot(122)
% go_netviewer_Matt(AEC_b,0.7)
% save AEC_b.mat AEC_b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
down_f = 10;
Nlocs = 78;
for seed = 1:78;
    Xsig = squeeze(VE(:,seed));
    for test = 1:78;
        if seed == test
            AEC(seed,test) = NaN;
        else
            Ysig = squeeze(VE(:,test));
            X_win = (Xsig - mean(Xsig));
            Y_win = (Ysig - mean(Ysig));
            %%regress leakage
            beta_leak = (pinv(X_win)*Y_win);
            Y_win_cor = Y_win - X_win*beta_leak;
            %%calculate envelopes
            H_X = abs(hilbert(X_win));
            H_X_d = mean(reshape(H_X,fs/down_f,duration*down_f,1));
            %%calculate envelopes
            H_Y = abs(hilbert(Y_win_cor));
            H_Y_d = mean(reshape(H_Y,fs/down_f,duration*down_f,1));
            AEC(seed,test) = corr(H_X_d',H_Y_d');
        end
    end
    clc;disp(sprintf('running beta connectivity region %d',seed));
end
AEC_b = 0.5*(AEC + AEC');
figure(1)
subplot(121)
imagesc(AEC_b);colorbar;
subplot(122)
go_netviewer_Matt(AEC_b,0.7)
save AEC_b.mat AEC_b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hp = 8;
lp = 12;
[b,a] = butter(4,2*[hp lp]/fs);
data_f = [filtfilt(b,a,data1)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strmatch(sub,'010') & strmatch(ses,'001')
    duration = 600
    epoch_length = 5;
    time_epoch = linspace(0, epoch_length, epoch_length*fs);
%     for n = 1:duration/epoch_length
%         figure(567);
%         clf
%         dat_epoch = data_f(:, (n-1)*(epoch_length*fs)+1:(n*epoch_length*fs));
%         imagesc(time_epoch,1:size(dat_epoch,1),dat_epoch);shading flat
%         colorbar
%         ylabel('channel')
%         xlabel('time (s)')
%         caxis([-500 500])
%         title(sprintf('Epoch %d',n));
%         pause(1)
%     end
    %% remove bad trials...
    data_f_mat = reshape(data_f,[size(data_f,1),epoch_length*fs,duration/epoch_length]);
    badtrials = [71 79 82 83 91]
    trials1 = ones(1,duration/epoch_length); trials1(badtrials) = 0;
    trials_to_keep = find(trials1 == 1);
    data_f_mat2 = data_f_mat(:,:,trials_to_keep);
    clear data_f
    clear duration
    duration = epoch_length*length(trials_to_keep)
    data_f = reshape(data_f_mat2,size(data_f_mat2,1),duration*fs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



C = cov(data_f');
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);
for n = 1:length(sourcepos)
     this_L = Lead_fields(:,:,n);
     W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
     iPower_v = this_L'*Cr_inv*this_L;
     [v,d] = svd(iPower_v);
     [~,id] = min(diag(d));
     lopt = this_L*v(:,id); % turn to nAm amplitude
     w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt));
     VEa(:,n) = w*data_f;
end
if exist(basefilename);
    cd(basefilename)
else
    mkdir(basefilename)
    cd(basefilename)
end
save VEa.mat VEa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
down_f = 10;
Nlocs = 78;
for seed = 1:78;
    Xsig = squeeze(VEa(:,seed));
    for test = 1:78;
        if seed == test
            AEC(seed,test) = NaN;
        else
            Ysig = squeeze(VEa(:,test));
            X_win = (Xsig - mean(Xsig));
            Y_win = (Ysig - mean(Ysig));
            %%regress leakage
            beta_leak = (pinv(X_win)*Y_win);
            Y_win_cor = Y_win - X_win*beta_leak;
            %%calculate envelopes
            H_X = abs(hilbert(X_win));
            H_X_d = mean(reshape(H_X,fs/down_f,duration*down_f,1));
            %%calculate envelopes
            H_Y = abs(hilbert(Y_win_cor));
            H_Y_d = mean(reshape(H_Y,fs/down_f,duration*down_f,1));
            AEC(seed,test) = corr(H_X_d',H_Y_d');
        end
    end
    fprintf('running alpha connectivity region %2d\n',seed);
end
AEC_a = 0.5*(AEC + AEC');
figure(242)
subplot(121)
imagesc(AEC_a);colorbar;
subplot(122)
go_netviewer_Matt(AEC_a,0.7)
save AEC_a.mat AEC_a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

