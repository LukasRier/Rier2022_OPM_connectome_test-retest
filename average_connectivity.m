clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_dir = 'R:\OPMMEG\Projects\movie\';
project_dir = '/net/cador/data_local/Lukas/movie/';
results_dir = [project_dir,'results',filesep,'average_connectomes',filesep];

datadir = [project_dir,'data',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpfs = [4, 8,13,30,35,40];
lpfs = [8,12,30,40,45,48];

if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end
band_name = {'\theta','\alpha','\beta','\gamma_1','\gamma_2','\gamma_3'};
mean_corr = 0;dosave = 0;
pctl = 0.95;

if mean_corr
    predef_lims = [2.2,2.8,2.8,5.5,5.5,10];
else
    predef_lims = [0.2,0.23,0.23,0.23,0.07,0.06];
end
for f_ind = 1:length(hpfs)
    clearvars -except project_dir predef_lims mean_corr pctl dosave results_dir datadir hpfs lpfs f_ind pval all_real_diffs all_mean_within_corrs band_name AEC_b_all mean_AEC_b_all_frq mean_AEC_b2_all_frq
    
    
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
            AEC_b_all(:,:,n) = AEC./mean(AEC(:));
        else
            AEC_b_all(:,:,n) = AEC;
        end
    end
    mean_AEC_b = mean(AEC_b_all,3);
    mean_AEC_b_all_frq{f_ind} = mean_AEC_b;
    
    %mean_AEC_a = mean(AEC_a_all,3);
    
    figure;
    set(gcf,'Position',[680 95 550 883])
    subplot(4,2,1:4)
    
    limit = prctile(mean_AEC_b(triu(ones(size(mean_AEC_b)),1)==1),pctl*100);
    mask = abs(mean_AEC_b) >= limit;
    mean_AEC_b(~mask)=0;
    imagesc(mean_AEC_b);colorbar;
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
    go_netviewer_perctl_lim(mean_AEC_b,pctl,predef_lims(f_ind))
    view([-180,0])
    ax2 = subplot(4,2,5);ax2.Position = [0.05,0.25,0.4,0.20];
    go_netviewer_perctl_lim(mean_AEC_b,pctl,predef_lims(f_ind))
    view([-90,0])
    ax1 = subplot(4,2,[6,8]);ax1.Position = [0.45,0.05,0.5,0.4];
    go_netviewer_perctl_lim(mean_AEC_b,pctl,predef_lims(f_ind))
    
    st = suptitle(band_name{f_ind});
    st.FontSize = 35;
    drawnow
    if dosave
    if mean_corr
        saveas(gcf,sprintf('%saverage_conn1_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
    else
        saveas(gcf,sprintf('%saverage_conn1_%d_%d_Hz.png',results_dir,hp,lp));
    end
    end
    
    %%
    
%     figure
%     C(1:8,:) = mean_AEC_b;
%     ax3 = subplot(4,2,7);ax3.Position = [0.05,0.05,0.4,0.20];
%     go_netviewer_perctl_lim(C,pctl,predef_lims(f_ind))
%     view([-180,0])
%     ax2 = subplot(4,2,5);ax2.Position = [0.05,0.25,0.4,0.20];
%     go_netviewer_perctl_lim(C,pctl,predef_lims(f_ind))
%     view([-90,0])
%     ax1 = subplot(4,2,[6,8]);ax1.Position = [0.45,0.05,0.5,0.4];
%     go_netviewer_perctl_lim(C,pctl,predef_lims(f_ind))

end
%
for f_ind = 1:length(hpfs)
    clearvars -except project_dir results_dir predef_lims pctl dosave mean_corr datadir hpfs lpfs f_ind pval all_real_diffs all_mean_within_corrs band_name mean_AEC_b_all_frq mean_AEC_b2_all_frq AEC_b_all2 mean_AEC_b AEC_b_all
    hp = hpfs(f_ind);
    lp = lpfs(f_ind);
    
    ses = '002';
    for n = 1:10
        sub = sprintf('%3d',n);sub(sub == ' ') = '0';
        path.AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
        files.AEC = ['sub-',sub,'_','run-',ses,'_AEC'];
        load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.AEC,files.AEC,hp,lp),'AEC')
        X = isnan(AEC);
        AEC(X) = 0;
        if mean_corr
            AEC_b_all2(:,:,n) = AEC./mean(AEC(:));
        else
            AEC_b_all2(:,:,n) = AEC;
        end
    end
    mean_AEC_b2 = mean(AEC_b_all2,3);
    mean_AEC_b2_all_frq{f_ind} = mean_AEC_b2;
    figure
    set(gcf,'Position',[680 95 550 883])
    subplot(4,2,1:4)
    
    limit = prctile(mean_AEC_b2(triu(ones(size(mean_AEC_b2)),1)==1),pctl*100);
    mask = abs(mean_AEC_b2) >= limit;
    mean_AEC_b2(~mask)=0;
    
    imagesc(mean_AEC_b2);colorbar;
%     cLims = [-max(abs(mean_AEC_b2(:))) max(abs(mean_AEC_b2(:)))];
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
    go_netviewer_perctl_lim(mean_AEC_b2,pctl,predef_lims(f_ind))
    view([-180,0])
    ax2 = subplot(4,2,5);ax2.Position = [0.05,0.25,0.4,0.20];
    go_netviewer_perctl_lim(mean_AEC_b2,pctl,predef_lims(f_ind))
    view([-90,0])
    ax1 = subplot(4,2,[6,8]);ax1.Position = [0.45,0.05,0.5,0.4];
    go_netviewer_perctl_lim(mean_AEC_b2,pctl,predef_lims(f_ind))
    
    
    st = suptitle(band_name{f_ind});
    st.FontSize = 35;
    drawnow
    
    if dosave
    if mean_corr
        saveas(gcf,sprintf('%saverage_conn2_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
    else
        saveas(gcf,sprintf('%saverage_conn2_%d_%d_Hz.png',results_dir,hp,lp));
    end
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for xx = 1:78
    for yy = 1:78
        a = squeeze(AEC_b_all(xx,yy,:));
        b = squeeze(AEC_b_all2(xx,yy,:));
        SR(xx,yy) = signrank(a-b);
    end
end
if mean_corr;mean_corr_str = ' meancorr';else;mean_corr_str = '';end

figure
set(gcf,'Position',[680 95 507 883],'Name',['sign rank test',mean_corr_str])
subplot(211)
imagesc(SR);colorbar;


axis square; yticks([5, 14, 25, 37, 44, 53, 64, 76]);
yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xticks([5, 14, 25, 37, 44, 53, 64, 76]);
xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
    'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
xtickangle(45)
subplot(212)
go_netviewer_perctl(SR,pctl)

drawnow
st = suptitle(band_name{f_ind});
st.FontSize = 35;
%%

for  f_ind = 1:length(hpfs)
    figure
    set(gcf,'Color','w')
    hp = hpfs(f_ind);
    lp = lpfs(f_ind);
    subplot(1,3,1)
    imagesc(mean_AEC_b_all_frq{f_ind}-mean_AEC_b2_all_frq{f_ind})
    colorbar;
    caxis([ -0.05    0.05])
    
    
    axis square; yticks([5, 14, 25, 37, 44, 53, 64, 76]);
    yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
        'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
    xticks([5, 14, 25, 37, 44, 53, 64, 76]);
    xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
        'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
    xtickangle(45)
    title('Difference')
    
    subplot(1,3,2)
    imagesc(mean_AEC_b_all_frq{f_ind})
    colorbar;
    cax = caxis;
    
    axis square; yticks([5, 14, 25, 37, 44, 53, 64, 76]);
    yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
        'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
    xticks([5, 14, 25, 37, 44, 53, 64, 76]);
    xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
        'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
    xtickangle(45)
    title('Run 1')
    
    subplot(1,3,3)
    imagesc(mean_AEC_b2_all_frq{f_ind})
    colorbar;
    
    
    axis square; yticks([5, 14, 25, 37, 44, 53, 64, 76]);
    yticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
        'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
    xticks([5, 14, 25, 37, 44, 53, 64, 76]);
    xticklabels({'L.M.Frontal', 'L.P.Motor', 'L.Calcarine',...
        'L.Cingulum', 'R.M.Frontal', 'R.P.Motor', 'R.Calcarine', 'R.Cingulum'});
    xtickangle(45)
    caxis(cax)
    set(gcf,'Position',[0,479,1983,520])
    title('Run 2')

    drawnow
    st = suptitle(band_name{f_ind});
    drawnow
    % st.FontSize = 355
    if dosave
    if mean_corr
        saveas(gcf,sprintf('%saverage_conn_and_diff_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
    else
        saveas(gcf,sprintf('%saverage_conn_and_diff_%d_%d_Hz.png',results_dir,hp,lp));
    end
    end
end
%%
f = figure(9999);clf
f.Color = 'w';
ax = axes;
hold on
c1 = [0.1216    0.4706    0.7059];
c2 = [  0.2000    0.6275    0.1725];


for  f_ind = 1:length(hpfs)
    mask=triu(ones(78),1)==1;
    bar(ax,f_ind+[-.2],mean(mean_AEC_b_all_frq{f_ind}(mask)),0.3,'FaceColor',c1)
    bar(ax,f_ind+[+.2],mean(mean_AEC_b2_all_frq{f_ind}(mask)),0.3,'FaceColor',c2)
    errorbar(ax,f_ind+[-.2],mean(mean_AEC_b_all_frq{f_ind}(mask)),std(mean_AEC_b_all_frq{f_ind}(mask)),'k')
    errorbar(ax,f_ind+[+.2],mean(mean_AEC_b2_all_frq{f_ind}(mask)),std(mean_AEC_b2_all_frq{f_ind}(mask)),'k')
end
ax.XTick=[1:length(hpfs)];
ax.XTickLabel = band_name;ax.FontSize=15
title(['Mean connectivity per band',mean_corr_str])
ylabel(['Mean global connectivity'])
legend('Run1','Run2','std. dev.')
if dosave
if mean_corr
    saveas(f,sprintf('%saverage_conn_per_run_meancorr.png',results_dir));
else
    saveas(f,sprintf('%saverage_conn_per_run.png',results_dir));
end
end
%%
f = figure(1111);clf
f.Color = 'w';
ax = axes;
hold on
frq_cols = [237,248,177
199,233,180
127,205,187
65,182,196
29,145,192
34,94,168]./256;frq_cols = flipdim(frq_cols,1);

for  f_ind = 1:length(hpfs)
    diffs = mean_AEC_b_all_frq{f_ind}-mean_AEC_b2_all_frq{f_ind};
    mask=triu(ones(78),1)==1;
    bar(ax,f_ind,mean(diffs(mask)),'FaceColor',frq_cols(f_ind,:))
    errorbar(ax,f_ind,mean(diffs(mask)),std(diffs(mask)),'k')
end
ax.XTick=[1:length(hpfs)];
ax.XTickLabel = band_name;ax.FontSize=15;
title(['Mean connectivity differences',mean_corr_str])
ylabel(sprintf('Mean Connectivity'))
% legend('Run1-Run2','std. dev.','Location','best')
if dosave
if mean_corr
    saveas(f,sprintf('%saverage_conn_differences_meancorr.png',results_dir));
else
    saveas(f,sprintf('%saverage_conn_differences.png',results_dir));
end
end
%%

fcorr = figure();clf
fcorr.Color = 'w';
fcorr.Position = [133         326        1747         669];

ax1 = subplot(2,length(hpfs),[length(hpfs)+1 : length(hpfs)*2]);
hold on
xlabel('Frequency (Hz)')
figure(fcorr);
for f_ind = 1:length(hpfs)
    [rho_,pval_] = corr(mean_AEC_b_all_frq{f_ind}(triu(ones(78),1)==1),mean_AEC_b2_all_frq{f_ind}(triu(ones(78),1)==1))
    pval(f_ind) = pval_;
    corr_vals(f_ind) = rho_;
    subplot(2,length(hpfs),f_ind)
    plot(mean_AEC_b_all_frq{f_ind}(triu(ones(78),1)==1),mean_AEC_b2_all_frq{f_ind}(triu(ones(78),1)==1),'k.')
    xlabel('Run 1')
    ylabel('Run 2')
    hold on
    xlims = xlim;
    title(band_name{f_ind});
    xlim(xlims);
    plot(xlim,xlim,'r')
    axis square
    
end
[~,sortid]= sort(pval);
crit_pval_benjamini = 0.05./(length(hpfs):-1:1);
mincolor = 0.3;
colorstep = 0.015;
for f_ind = 1:length(hpfs)


center = (lpfs(f_ind) + hpfs(f_ind))./2;
width = (lpfs(f_ind) - hpfs(f_ind))*0.975;
val = corr_vals(f_ind);
bh1(f_ind) = bar(ax1,center,val,width,'FaceColor',frq_cols(f_ind,:),'FaceAlpha',0.8);

if pval(f_ind) < crit_pval_benjamini(sortid(f_ind))
    th(f_ind) = text(ax1, center, val + 0.1, '*','FontSize',15,'HorizontalAlignment','center');
    bh1(f_ind).FaceAlpha = 1;
end
end
% bh1(5).Visible = 'off';
% bh2(5).Visible = 'off';
% th(5).Visible = 'off';
ax1.YLabel.String = sprintf('Between run correlation');
set([ax1],'FontSize',15,'TickLength',[0,0],'XTickLabelRotation',70,'XTick',unique([hpfs,lpfs]));
drawnow
if dosave
if mean_corr
    saveas(fcorr,sprintf('%saverage_conn_betw_run_correlations_meancorr.png',results_dir));
else
    saveas(fcorr,sprintf('%saverage_conn_betw_run_correlations.png',results_dir));
end
end