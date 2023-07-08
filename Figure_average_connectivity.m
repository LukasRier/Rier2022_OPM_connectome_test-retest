clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set project dir to directory containing data from https://doi.org/10.5072/zenodo.1134455
project_dir = 'F:\Rdrive\movie\';
if ~exist(project_dir,'dir')
    error('Set project directory!')
end
results_dir = [project_dir,'results',filesep,'average_connectomes',filesep];
mkdir(results_dir)
datadir = [project_dir,'data',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpfs = [4, 8,13,30,35,40];
lpfs = [8,12,30,40,45,48];

if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end
band_name = {'\theta','\alpha','\beta','\gamma_1','\gamma_2','\gamma_3'};
dosave = 1;
pctl = 0.95;
if dosave
resf = fopen([results_dir,'Average_connectivity_results.csv'],'w');
fprintf(resf,'Current plot included edges percentile = %1.2f\n',100*(1-pctl));
fprintf(resf,'Measure,mean,std.dev,p\n');
end
mean_corr = 1;
% if mean_corr
%     predef_lims = [2.3,2.8,3,5.5,5.5,15];
% else
%     predef_lims = [0.2,0.23,0.23,0.23,0.07,0.06];
% end
for f_ind = 1:length(hpfs)
    
    
    hp = hpfs(f_ind);
    lp = lpfs(f_ind);
    mask=triu(ones(78),1)==1;

    ses = '001';
    for n = 1:10
        sub = sprintf('%3d',n);sub(sub == ' ') = '0';
        path.AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
        files.AEC = ['sub-',sub,'_','run-',ses,'_AEC'];
        load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.AEC,files.AEC,hp,lp),'AEC')
        X = isnan(AEC);
        AEC(X) = 0;
        
        AEC_b_all_corr(:,:,n) = AEC./sqrt(mean(AEC(mask).^2));
        
        AEC_b_all(:,:,n) = AEC;
        mask=triu(ones(78),1)==1;
        global_conn_1(f_ind,n) = mean(AEC(mask));
    end
    mean_AEC_b_corr = mean(AEC_b_all_corr,3);
    mean_AEC_b_all_frq_corr{f_ind} = mean_AEC_b_corr;
    AEC_b_all_frq{f_ind} = AEC_b_all;
        
%     plot_mat_n_brains(mean_AEC_b_corr,predef_lims,f_ind,pctl,band_name)
%     
%     if dosave
%         saveas(gcf,sprintf('%saverage_conn1_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
%     end
%     
end

%  
mask=triu(ones(78),1)==1;

for f_ind = 1:length(hpfs)
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
        AEC_b_all2_corr(:,:,n) = AEC./sqrt(mean(AEC(mask).^2));
        AEC_b_all2(:,:,n) = AEC;
        
        
        global_conn_2(f_ind,n) = mean(AEC(mask));
    end
    mean_AEC_b2_corr = mean(AEC_b_all2_corr,3);
    mean_AEC_b2_all_frq_corr{f_ind} = mean_AEC_b2_corr;
    AEC_b2_all_frq{f_ind} = AEC_b_all2;
    
%     plot_mat_n_brains(mean_AEC_b2_corr,predef_lims,f_ind,pctl,band_name)
%     
%     drawnow
%     
%     if dosave
%         saveas(gcf,sprintf('%saverage_conn2_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
%     end
end


%% find limits 
for f_ind = 1:length(hpfs)
    mean_AEC_b_corr = mean_AEC_b_all_frq_corr{f_ind};
    mean_AEC_b2_corr = mean_AEC_b2_all_frq_corr{f_ind};
    predef_lims(f_ind)= max(abs([min([mean_AEC_b_corr(mask);mean_AEC_b2_corr(mask)]),max([mean_AEC_b_corr(mask);mean_AEC_b2_corr(mask)])]));
end
for f_ind = 1:length(hpfs)
    hp = hpfs(f_ind);
    lp = lpfs(f_ind);
    
    mean_AEC_b_corr = mean_AEC_b_all_frq_corr{f_ind};

    plot_mat_n_brains(mean_AEC_b_corr,predef_lims,f_ind,pctl,band_name)
    
    if dosave
        saveas(gcf,sprintf('%saverage_conn1_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
        saveas(gcf,sprintf('%saverage_conn1_%d_%d_Hz_meancorr.svg',results_dir,hp,lp));
    end
    %%
    
    mean_AEC_b2_corr = mean_AEC_b2_all_frq_corr{f_ind};

    plot_mat_n_brains(mean_AEC_b2_corr,predef_lims,f_ind,pctl,band_name)
       
    if dosave
        saveas(gcf,sprintf('%saverage_conn2_%d_%d_Hz_meancorr.png',results_dir,hp,lp));
        saveas(gcf,sprintf('%saverage_conn2_%d_%d_Hz_meancorr.svg',results_dir,hp,lp));
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% global mean connectivity figure
f = figure(9999);clf;f.Renderer = 'painters';
f.Color = 'w';
f.Units = 'centimeters';
fwidth = 8;
fheight = 4;
f.Position([3,4]) = [fwidth,fheight];
ax = axes;
hold on
c1 = [0.1216    0.4706    0.7059];
c2 = [  0.2000    0.6275    0.1725];


for  f_ind = 1:length(hpfs)
    mask=triu(ones(78),1)==1;
    mean_run1 = mean(global_conn_1(f_ind,:));
    mean_run2 = mean(global_conn_2(f_ind,:));
    std_run1 = std(global_conn_1(f_ind,:));
    std_run2 = std(global_conn_2(f_ind,:));
  
    bar(ax,f_ind+[-.2],mean_run1,0.3,'FaceColor',[1,1,1].*0.43)
    bar(ax,f_ind+[+.2],mean_run2,0.3,'FaceColor',[1,1,1].*0.86)
    errorbar(ax,f_ind+[-.2],mean_run1,std_run1,'k')
    errorbar(ax,f_ind+[+.2],mean_run2,std_run2,'k')
    
    if dosave
    fprintf(resf,'Global conn. run 1 (%d-%dHz),%1.8f,%1.8f,\n',hpfs(f_ind),lpfs(f_ind),mean_run1,std_run1);
    fprintf(resf,'Global conn. run 2 (%d-%dHz),%1.8f,%1.8f,\n',hpfs(f_ind),lpfs(f_ind),mean_run2,std_run2);
    [~,p,~,stats] = ttest(global_conn_1(f_ind,:),global_conn_2(f_ind,:));
    fprintf(resf,'T-val global diff (%d-%dHz),%1.8f,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),stats.tstat,p);
    fprintf('T-val global diff (%d-%dHz),%1.8f,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),stats.tstat,p);
    
    p = signrank(global_conn_1(f_ind,:),global_conn_2(f_ind,:));
    fprintf(resf,'signrank global diff (%d-%dHz),,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),p);
    fprintf('signrank global diff (%d-%dHz),,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),p);
    end
end
ax.XTick=[1:length(hpfs)];
ax.XTickLabel = band_name;

% title(['Mean connectivity per band'])
ylabel(sprintf('Global mean\nconnectivity'))
legend('Run1','Run2','std. dev.')
f.Renderer = 'painters';
set(f.Children,'FontSize',9);

if dosave
    saveas(f,sprintf('%saverage_conn_per_run.png',results_dir));
    saveas(f,sprintf('%saverage_conn_per_run.svg',results_dir));
end
%% paired differences plot
f = figure(1111);clf
f.Color = 'w';f.Renderer = 'painters';
f.Units = 'centimeters';
fwidth = 8;
fheight = 4;
f.Position([3,4]) = [fwidth,fheight];
ax = axes;
hold on
frq_cols = [237,248,177
199,233,180
127,205,187
65,182,196
29,145,192
34,94,168]./256;frq_cols = flipdim(frq_cols,1);

for  f_ind = 1:length(hpfs)  
    diffs =  AEC_b_all_frq{f_ind} - AEC_b2_all_frq{f_ind};
    for subi = 1:size(diffs,3)
        sub_diff = diffs(:,:,subi);
        mean_diffs(subi) = mean(sub_diff(mask));
    end
    bar(ax,f_ind,mean(mean_diffs),'FaceColor',frq_cols(f_ind,:))
    errorbar(ax,f_ind,mean(mean_diffs),std(mean_diffs),'k')
    p = signrank(mean_diffs);
    if dosave; fprintf(resf,'signrank mean of paired diffs (%d-%dHz),,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),p);end;
    fprintf('signrank mean of paired diffs (%d-%dHz),,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),p);

end
ax.XTick=[1:length(hpfs)];
ax.XTickLabel = band_name;
ylabel(sprintf(['Mean Connectivity\ndifference']))
f.Renderer = 'painters';
set(f.Children,'FontSize',9);
% legend('Run1-Run2','std. dev.','Location','best')
if dosave
    saveas(f,sprintf('%saverage_conn_differences.png',results_dir));
    saveas(f,sprintf('%saverage_conn_differences.svg',results_dir));
end
%%

fcorr = figure();clf
fcorr.Color = 'w';fcorr.Renderer = 'painters';fcorr.Units = 'centimeters';
tit_fs = 10.8
ax1 = subplot(3,length(hpfs),[length(hpfs)+1 : length(hpfs)*2]);
hold on
xlabel('Frequency (Hz)')
figure(fcorr);
for f_ind = 1:length(hpfs)
    [rho_,pval_] = corr(mean_AEC_b_all_frq_corr{f_ind}(triu(ones(78),1)==1),mean_AEC_b2_all_frq_corr{f_ind}(triu(ones(78),1)==1));
    pval(f_ind) = pval_;
    corr_vals(f_ind) = rho_;
    subplot(3,length(hpfs),f_ind)
    plot(mean_AEC_b_all_frq_corr{f_ind}(triu(ones(78),1)==1),mean_AEC_b2_all_frq_corr{f_ind}(triu(ones(78),1)==1),'k.','MarkerSize',0.2)
    xlabel('Run 1')
    ylabel('Run 2')
    hold on
    xlims = ylim;
    title(band_name{f_ind},'FontSize',tit_fs);
    xlim(xlims);ylim(xlims)
    plot(xlim,xlim,'r')
    axis square
%     set(gca,'FontSize',11)

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
%     th(f_ind) = text(ax1, center, val + 0.1, '*','FontSize',15,'HorizontalAlignment','center');
    bh1(f_ind).FaceAlpha = 1;
end
if dosave;fprintf(resf,'Average conn correlation (%d-%dHz),%1.8f,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),corr_vals(f_ind),pval(f_ind));end
fprintf('Average conn correlation (%d-%dHz),%1.8f,p=,%1.8f\n',hpfs(f_ind),lpfs(f_ind),corr_vals(f_ind),pval(f_ind));
end
% bh1(5).Visible = 'off';
% bh2(5).Visible = 'off';
% th(5).Visible = 'off';
ax1.YLabel.String = sprintf('Between-run correlation');
% set([ax1],'FontSize',11,'TickLength',[0,0],'XTickLabelRotation',70,'XTick',unique([hpfs,lpfs]));
set([ax1],'TickLength',[0,0],'XTickLabelRotation',70,'XTick',unique([hpfs,lpfs]));

bs_Ns = 2:10;
% Impact of N on group average, between-run correlation
corrs=repmat({},length(hpfs),length(bs_Ns));
mean_corrs = zeros(length(hpfs),length(bs_Ns));
std_corrs = zeros(length(hpfs),length(bs_Ns));


for f_ind = 1:length(hpfs)
    for bs_ind = 1:length(bs_Ns)
        n_bs = bs_Ns(bs_ind);
       [choices] = nchoosek(1:10,n_bs);
        for ch_i = 1:size(choices,1)
            g1 = mean(AEC_b_all_frq{f_ind}(:,:,choices(ch_i,:)),3);
            g2 = mean(AEC_b2_all_frq{f_ind}(:,:,choices(ch_i,:)),3);
            corrs{f_ind,bs_ind}(ch_i) = corr(g1(mask),g2(mask));
        end
        mean_corrs(f_ind,bs_ind) = mean(corrs{f_ind,bs_ind});
        std_corrs(f_ind,bs_ind) = std(corrs{f_ind,bs_ind});
    end
    subplot(3,length(hpfs),2*length(hpfs)+f_ind)
    plot(bs_Ns,mean_corrs(f_ind,:),'k-x')
    hold on
    errorbar(bs_Ns(1:end-1),mean_corrs(f_ind,1:end-1),std_corrs(f_ind,1:end-1),'k')
    title(band_name(f_ind),'FontSize',tit_fs)
    xticks(bs_Ns)
    xlabel('Group size')
    ylabel('Average corr.')
    ylim([0,1])
    xlim([1.5,10.5])
%     set(gca,'FontSize',11)
end
fcorr.Position([3,4]) = [25 19];

drawnow
fcorr.Renderer = 'painters';

if dosave
    saveas(fcorr,sprintf('%saverage_conn_betw_run_correlations_meancorr.png',results_dir));
    saveas(fcorr,sprintf('%saverage_conn_betw_run_correlations_meancorr.svg',results_dir));
end


%%
function plot_mat_n_brains(AECmat,predef_lims,f_ind,pctl,band_name)
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
% st = sgtitle(band_name{f_ind});
drawnow
set(fh,'Renderer','painters')

end