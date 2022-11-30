clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_dir = 'F:\Rdrive\movie\';
% project_dir = '/net/cador/data_local/Lukas/movie/';

results_dir = [project_dir,'results',filesep,'individual_connectomes',filesep];
mkdir(results_dir);
datadir = [project_dir,'data',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpfs = [4, 8,13,30,35,40];
lpfs = [8,12,30,40,45,48];

if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end
band_name = {'\theta','\alpha','\beta','\gamma_1','\gamma_2','\gamma_3'};
mean_corr = 1;dosave = 1;
pctl = 0.95;
predef_lims = [2.1489, 2.5776, 2.9694, 3.8135, 2.3530, 2.2000];

mask=triu(ones(78),1)==1;

for f_ind = [2,1,3:6]
    if ~exist([results_dir,band_name{f_ind}(2:end),filesep],'dir')
        mkdir([results_dir,band_name{f_ind}(2:end),filesep]);
    end
    hp = hpfs(f_ind);
    lp = lpfs(f_ind);
    
    ses = '001';
    for n = 1:10
        sub = sprintf('%3d',n);sub(sub == ' ') = '0';
        path.AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
        files.AEC = ['sub-',sub,'_','run-',ses,'_AEC'];
        load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.AEC,files.AEC,hp,lp),'AEC')
        X = isnan(AEC);
        AEC(X) = 0;
        if mean_corr
            AEC_1_all(:,:,n) = AEC./sqrt(mean(AEC(mask).^2));
        else
            AEC_1_all(:,:,n) = AEC;
        end
    end
    
    ses = '002';
    for n = 1:10
        sub = sprintf('%3d',n);sub(sub == ' ') = '0';
        path.AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
        files.AEC = ['sub-',sub,'_','run-',ses,'_AEC'];
        load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.AEC,files.AEC,hp,lp),'AEC')
        X = isnan(AEC);
        AEC(X) = 0;
        if mean_corr
            AEC_2_all(:,:,n) = AEC./sqrt(mean(AEC(mask).^2));
        else
            AEC_2_all(:,:,n) = AEC;
        end
    end
    
    
    
    %first run
    for subi = 1:10
        AECmat = AEC_1_all(:,:,subi);
        %%
        plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,band_name)
        %
        sub = sprintf('%3d',subi);sub(sub == ' ') = '0';
        f = gcf;f.Renderer = 'painters';
        f.Children(3).Children.String = ['Subject ',sub];%f.Children(3).Children.FontSize = 25;
         %%
        if dosave
            if mean_corr
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz_meancorr.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz_meancorr.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
            else
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
            end
        end
    end
    %second run
    for subi = 1:10
        
        AECmat = AEC_2_all(:,:,subi);
        plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,band_name)
        
        sub = sprintf('%3d',subi);sub(sub == ' ') = '0';
        f = gcf;f.Renderer = 'painters';
        f.Children(3).Children.String = ['Subject ',sub];%f.Children(3).Children.FontSize = 25;
        if dosave
            if mean_corr
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz_meancorr.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz_meancorr.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
            else
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
            end
        end
    end
end

function plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,band_name)
figure;
set(gcf,'Position',[680 95 550 883],'Color','w','Renderer','painters')
subplot(4,2,1:4)
imagesc(AECmat);colorbar;
%     cLims = [-max(abs(mean_AEC_b(:))) max(abs(mean_AEC_b(:)))];
cLims = [-1,1]*predef_lims(f_ind);
caxis(cLims);
axis square; yticks([5, 14, 25, 37, 44, 53, 64, 76]);
yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xticks([5, 14, 25, 37, 44, 53, 64, 76]);
xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xtickangle(45)
ax=gca;
ax.FontSize=15; ax.Position([1,2]) = [0.15,0.67];

drawnow
ax3 = subplot(4,2,7);        ax3.Position = [0.05,0.005,0.45,0.24];
go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
view([-180,0]) %FRONT VIEW
ax2 = subplot(4,2,5);        ax2.Position = [0.05,0.24,0.42,0.24];
go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
view([-90,0]) %side view
ax1 = subplot(4,2,[6,8]);ax1.Position = [0.40,0,0.55,0.55];
go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind)) %top view
st = suptitle(band_name{f_ind});
st.FontSize = 25;
drawnow
ax1.Position = [0.40,0,0.55,0.55];
ax2.Position = [0.05,0.20,0.42,0.24];
ax3.Position = [0.05,0.005,0.45,0.24];
drawnow
set(gcf,'Renderer','painters')

end
