clear all
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set project dir to directory containing data from https://doi.org/10.5072/zenodo.1134455
project_dir = '/path/to/data/folder';
project_dir = 'F:\Rdrive\movie\';

if ~exist(project_dir,'dir')
    error('Set project directory!')
end
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
        sub = sprintf('%3d',subi);sub(sub == ' ') = '0';
        title_str = ['Subject ',sub];
        plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,title_str)
        
         %%
        if dosave
            if mean_corr
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz_meancorr.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                print(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz_meancorr.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp),'-dpng','-r600');
            else
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                print(gcf,sprintf('%s%s%ssub-%s-indiv_conn1_%d_%d_Hz.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp),'-dpng','-r600');
            end
        end
    end
    %second run
    for subi = 1:10
        
        AECmat = AEC_2_all(:,:,subi);
        sub = sprintf('%3d',subi);sub(sub == ' ') = '0';
        title_str = ['Subject ',sub,'_2'];
        plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,title_str)
        
        if dosave
            if mean_corr
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz_meancorr.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
                print(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz_meancorr.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp),'-dpng','-r600');
            else
                print(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz.png',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp),'-dpng','-r600');
                saveas(gcf,sprintf('%s%s%ssub-%s-indiv_conn2_%d_%d_Hz.svg',results_dir,band_name{f_ind}(2:end),filesep,sub,hp,lp));
            end
        end
    end
end

function plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,title_str)
fh = figure;
set(fh,'Units','centimeters','Color','w','Renderer','painters');
fwidth = 3.8;
fheight = 5.9;
fh.Position([3,4]) = [fwidth,fheight];
ax_mat = axes;
imagesc(AECmat);cb = colorbar('Location','southoutside','FontSize',8);
%     cLims = [-max(abs(mean_AEC_b(:))) max(abs(mean_AEC_b(:)))];
cLims = [-1,1]*predef_lims(f_ind);
caxis(cLims);
axis square; 
yticks([5, 14, 25, 37, 44, 53, 64, 76]);
yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xticks([5, 14, 25, 37, 44, 53, 64, 76]);
xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
% xtickangle(45)
ax_mat.Position = [0.46,0.58,0.52,0.52];
ax_mat.FontSize=6;
% = [0.15,0.67];

drawnow
ax3 = axes;        ax3.Position = [-0.1,0.0,0.6,0.3];
go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
view([-180,0]) %FRONT VIEW
ax2 = axes;
ax2.Position = [0.4,0.0,0.6,0.4];
go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind))
view([-90,0]) %side view
ax1 = axes;ax1.Position = [-0.05,0.25,0.55,0.4];
go_netviewer_perctl_lim(AECmat,pctl,predef_lims(f_ind)) %top view
drawnow

set(fh.Children(2:end),'Units','centimeters')
fh.Position(4) = fh.Position(4)+0.5;

st = annotation('textbox',[0,0.95,1,0.05],'String',title_str,...
        'FontSize',9,...
        'EdgeColor','w',...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'Margin',0);
set(fh,'Renderer','painters')

end