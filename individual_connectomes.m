clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_dir = 'R:\OPMMEG\Projects\movie\';
project_dir = '/net/cador/data_local/Lukas/movie/';

results_dir = [project_dir,'results',filesep];
datadir = [project_dir,'data',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpfs = [4, 8,13,30,35,40];
lpfs = [8,12,30,40,45,48];

if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end
band_name = {'\theta','\alpha','\beta','\gamma_1','\gamma_2','\gamma_3'};
mean_corr = 1;dosave = 1;
pctl = 0.95;
predef_lims = [2.2,2.8,2.8,5.5,5.5,10];


for f_ind = [1,3:length(hpfs)]
    
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
            AEC_1_all(:,:,n) = AEC./mean(AEC(:));
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
            AEC_2_all(:,:,n) = AEC./mean(AEC(:));
        else
            AEC_2_all(:,:,n) = AEC;
        end
    end
    
    
    
    %first run
    for subi = 1:10
        figure;
        set(gcf,'Position',[680 95 550 883],'Color','w')
        AECmat = AEC_1_all(:,:,subi);
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
        ax3 = subplot(4,2,7);ax3.Position = [0.05,0.05,0.4,0.20];
        go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
        view([-180,0])
        ax2 = subplot(4,2,5);ax2.Position = [0.05,0.25,0.4,0.20];
        go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
        view([-90,0])
        ax1 = subplot(4,2,[6,8]);ax1.Position = [0.45,0.05,0.5,0.4];
        go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
        
        st = suptitle(band_name{f_ind});
        st.FontSize = 35;
        drawnow
        sub = sprintf('%3d',subi);sub(sub == ' ') = '0';

        if dosave
            if mean_corr
                saveas(gcf,sprintf('%ssub-%s-indiv_conn1_%d_%d_Hz_meancorr.png',results_dir,sub,hp,lp));
            else
                saveas(gcf,sprintf('%ssub-%s-indiv_conn1_%d_%d_Hz.png',results_dir,sub,hp,lp));
            end
        end
    end
    %second run
    for subi = 1:10
        figure;
        set(gcf,'Position',[680 95 550 883],'Color','w')
        AECmat = AEC_2_all(:,:,subi);
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
        ax3 = subplot(4,2,7);ax3.Position = [0.05,0.05,0.4,0.20];
        go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
        view([-180,0])
        ax2 = subplot(4,2,5);ax2.Position = [0.05,0.25,0.4,0.20];
        go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
        view([-90,0])
        ax1 = subplot(4,2,[6,8]);ax1.Position = [0.45,0.05,0.5,0.4];
        go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
        
        st = suptitle(band_name{f_ind});
        st.FontSize = 35;
        drawnow
        sub = sprintf('%3d',subi);sub(sub == ' ') = '0';

        if dosave
            if mean_corr
                saveas(gcf,sprintf('%ssub-%s-indiv_conn2_%d_%d_Hz_meancorr.png',results_dir,sub,hp,lp));
            else
                saveas(gcf,sprintf('%ssub-%s-indiv_conn2_%d_%d_Hz.png',results_dir,sub,hp,lp));
            end
        end
    end
end
