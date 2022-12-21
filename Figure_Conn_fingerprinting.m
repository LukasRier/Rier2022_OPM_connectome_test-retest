clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set project dir to directory containing data from https://doi.org/10.5072/zenodo.1134455
project_dir = '/path/to/data/folder';
if ~exist(project_dir,'dir')
    error('Set project directory!')
end

results_dir = [project_dir,'results',filesep,'connectome_fingerprints',filesep];
mkdir(results_dir);
datadir = [project_dir,'data',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hpfs = [4, 8,13,30,35,40];
lpfs = [8,12,30,40,45,48];
band_name = {'\theta','\alpha','\beta','\gamma_1','\gamma_2','\gamma_3'};

if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end
frq_cols = [237,248,177
199,233,180
127,205,187
65,182,196
29,145,192
34,94,168]./256;frq_cols = flipdim(frq_cols,1);
resf = fopen([results_dir,'FingerPrinting_results.csv'],'w')
fprintf(resf,'Measure,mean,std.dev.\n');
for f_ind = 1:length(hpfs)
close all
clearvars -except resf min_corr max_corr  project_dir results_dir datadir hpfs lpfs f_ind frq_cols pval all_real_diffs all_mean_within_corrs band_name


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
    AEC_b_all(:,:,n) = AEC./sqrt(mean(AEC(mask).^2));
end
%%
ses = '002';
for n = 1:10
    sub = sprintf('%3d',n);sub(sub == ' ') = '0';
    path.AEC = [datadir,'derivatives',filesep,'AEC',filesep,'sub-',sub,filesep];
    files.AEC = ['sub-',sub,'_','run-',ses,'_AEC'];
    load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.AEC,files.AEC,hp,lp),'AEC')
    X = isnan(AEC);
    AEC(X) = 0;
    AEC_b_all2(:,:,n) = AEC./sqrt(mean(AEC(mask).^2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% neural fingerprinting

%% calculate the correlations
count = 0;
for n = 1:10
    for m = 1:10
        count = count + 1;
        run1sa = AEC_b_all(:,:,n);
        run2sb = AEC_b_all2(:,:,m);
        run1sa_v_t = run1sa(mask);
        run2sb_v_t = run2sb(mask);
        sub_correl(n,m) = corr(run1sa_v_t,run2sb_v_t);
    end
end
min_corr = -0.15 ;if min(sub_correl(:))<-0.15;error('change lims');end
max_corr = 1; %max(sub_correl(:))
%% Fingerprinting matrix
figure(32767);
pcolor(1:10,1:10,sub_correl);shading flat
xticks([1 2 3 4 5 6 7 8 9 10]);
xticklabels({'Subject 1 Run 1', 'Subject 2 Run 1', 'Subject 3 Run 1', 'Subject 4 Run 1', 'Subject 5 Run 1', 'Subject 6 Run 1', 'Subject 7 Run 1', 'Subject 8 Run 1', 'Subject 9 Run 1', 'Subject 10 Run 1'});
xtickangle(45)
yticklabels({'Subject 1 Run 2', 'Subject 2 Run 2', 'Subject 3 Run 2', 'Subject 4 Run 2', 'Subject 5 Run 2', 'Subject 6 Run 2', 'Subject 7 Run 2', 'Subject 8 Run 2', 'Subject 9 Run 2', 'Subject 10 Run 2'});
ytickangle(45)

[~,max_inds] = max(sub_correl,[],1);
correctly_identified = sum(max_inds == 1:10);
fprintf(resf,'Fingerprinting success (%d-%dHz),%d/10,\n',hp,lp,correctly_identified);

set(gcf,'color',[1 1 1])
saveas(gcf,sprintf('%sFingerprinting_matrix_%d_%d_Hz.png',results_dir,hp,lp));

%% subject 1
S1_val = sub_correl(1,1);
S1_bet = [sub_correl(1,2:10) sub_correl(2:10,1)']
%% subject 2
S2_val = sub_correl(2,2);
S2_bet = [sub_correl(2,1) sub_correl(2,3:10) sub_correl(1,2)' sub_correl(3:10,2)']
%% subject 3
S3_val = sub_correl(3,3);
S3_bet = [sub_correl(3,1:2) sub_correl(3,4:10) sub_correl(1:2,3)' sub_correl(4:10,3)']
%% subject 4
S4_val = sub_correl(4,4);
S4_bet = [sub_correl(4,1:3) sub_correl(4,5:10) sub_correl(1:3,4)' sub_correl(5:10,4)']
%% subject 5
S5_val = sub_correl(5,5);
S5_bet = [sub_correl(5,1:4) sub_correl(5,6:10) sub_correl(1:4,5)' sub_correl(6:10,5)']
%% subject 6
S6_val = sub_correl(6,6);
S6_bet = [sub_correl(6,1:5) sub_correl(6,7:10) sub_correl(1:5,6)' sub_correl(7:10,6)']
%% subject 7
S7_val = sub_correl(7,7);
S7_bet = [sub_correl(7,1:6) sub_correl(7,8:10) sub_correl(1:6,7)' sub_correl(8:10,7)']
%% subject 8
S8_val = sub_correl(8,8);
S8_bet = [sub_correl(8,1:7) sub_correl(8,9:10) sub_correl(1:7,8)' sub_correl(9:10,8)']
%% subject 9
S9_val = sub_correl(9,9);
S9_bet = [sub_correl(9,1:8) sub_correl(9,10) sub_correl(1:8,9)' sub_correl(10,9)']
%% subject 10
S10_val = sub_correl(10,10);
S10_bet = [sub_correl(10,1:9) sub_correl(1:9,10)']

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103)
set(gcf,'position',[10,10,500,700])
%% subject 1
plot(ones(size(S1_bet))+0.3*(rand(size(S1_bet))-0.5),S1_bet,'k.','linewidth',1)
hold on
plot(1,S1_val,'r+','linewidth',4)
%% subject 2
plot(2*ones(size(S2_bet))+0.3*(rand(size(S2_bet))-0.5),S2_bet,'k.','linewidth',1)
hold on
plot(2,S2_val,'r+','linewidth',4)
%% subject 3
plot(3*ones(size(S3_bet))+0.3*(rand(size(S3_bet))-0.5),S3_bet,'k.','linewidth',1)
hold on
plot(3,S3_val,'r+','linewidth',4)
%% subject 4
plot(4*ones(size(S4_bet))+0.3*(rand(size(S4_bet))-0.5),S4_bet,'k.','linewidth',1)
hold on
plot(4,S4_val,'r+','linewidth',4)
%% subject 5
plot(5*ones(size(S5_bet))+0.3*(rand(size(S5_bet))-0.5),S5_bet,'k.','linewidth',1)
hold on
plot(5,S5_val,'r+','linewidth',4)
%% subject 6
plot(6*ones(size(S6_bet))+0.3*(rand(size(S6_bet))-0.5),S6_bet,'k.','linewidth',1)
hold on
plot(6,S6_val,'r+','linewidth',4)
%% subject 7
plot(7*ones(size(S7_bet))+0.3*(rand(size(S7_bet))-0.5),S7_bet,'k.','linewidth',1)
hold on
plot(7,S7_val,'r+','linewidth',4)
%% subject 8
plot(8*ones(size(S8_bet))+0.3*(rand(size(S8_bet))-0.5),S8_bet,'k.','linewidth',1)
hold on
plot(8,S8_val,'r+','linewidth',4)
%% subject 9
plot(9*ones(size(S9_bet))+0.3*(rand(size(S9_bet))-0.5),S9_bet,'k.','linewidth',1)
hold on
plot(9,S9_val,'r+','linewidth',4)
%% subject 10
plot(10*ones(size(S10_bet))+0.3*(rand(size(S10_bet))-0.5),S10_bet,'k.','linewidth',1)
hold on
plot(10,S10_val,'r+','linewidth',4)

axis([0 11 min_corr max_corr])
grid on
ylabel('Correlation')
xticks([1 2 3 4 5 6 7 8 9 10]);
xticklabels({'Subject 1', 'Subject 2', 'Subject 3', 'Subject 4', 'Subject 5', 'Subject 6', 'Subject 7', 'Subject 8', 'Subject 9', 'Subject 10'});
xtickangle(45)
set(gcf,'color',[1 1 1])
legend('between subject','Within subject','location','southwest')
set(gca,'fontsize',14)
saveas(gcf,sprintf('%sFingerprinting_corr_scatter_%d_%d_Hz.png',results_dir,hp,lp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(104);clf
set(gcf,'position',[10,10,500,400])
allwithin = diag(sub_correl);
all_mean_within_corrs(f_ind) = mean(allwithin);
all_between = [S1_bet S2_bet S3_bet S4_bet S5_bet S6_bet S7_bet S8_bet];
bar([1:2],[mean(allwithin) mean(all_between)],'FaceColor',frq_cols(f_ind,:));
hold on
plot(1*ones(size(allwithin))+0.5*(rand(size(allwithin))-0.5),allwithin,'r+','linewidth',4)
plot(2*ones(size(all_between))+0.5*(rand(size(all_between))-0.5),all_between,'k.','linewidth',1)
ylabel('Correlation')
xticks([1 2]);
xticklabels({'Within subject', 'Between subject'});
set(gca,'fontsize',14)
set(gcf,'color',[1 1 1])
xtickangle(0)
axis([0.5 2.5 min_corr max_corr])
tx = text(gca,1.5,0.9,sprintf('%s',band_name{f_ind}),'FontSize',25,...
    'HorizontalAlignment','center','VerticalAlignment','middle');
saveas(gcf,sprintf('%sFingerprinting_within_between_bar_chart_%d_%d_Hz.png',results_dir,hp,lp));
saveas(gcf,sprintf('%sFingerprinting_within_between_bar_chart_%d_%d_Hz.svg',results_dir,hp,lp));
fprintf(resf,'Within-sub corr.(%d-%dHz),%1.8f,%1.8f\n',hp,lp,mean(allwithin),std(allwithin));
fprintf(resf,'Between-sub corr.(%d-%dHz),%1.8f,%1.8f\n',hp,lp,mean(all_between),std(all_between));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% statistical test 
realdiff = mean(allwithin) - mean(all_between)
all_real_diffs(f_ind) = realdiff;
alldata = [allwithin' all_between];
for n = 1:100000
    alldata_shuff = alldata(randperm(length(alldata)));
    allwithin_sham = alldata_shuff(1:10);
    all_between_sham = alldata_shuff(11:end);
    shamdiff(n) = mean(allwithin_sham) - mean(all_between_sham);
end
figure(105)
set(gcf,'position',[10,10,500,400])
[null,bins] = hist(shamdiff,100);
% pval = sum(abs(shamdiff) > abs(realdiff))./length(shamdiff);
pval(f_ind) = sum(shamdiff > realdiff)./length(shamdiff);%onesided

hold on
plot(bins,null);
hold on
plot(realdiff,0,'r+','linewidth',4);
xlabel('Correlation difference')
ylabel('Count')
set(gcf,'color',[1 1 1])
legend('null distribution','real value','location','northeast')
title(sprintf('P = %1.5f',pval(f_ind)))
set(gca,'fontsize',14)
saveas(gcf,sprintf('%sFingerprinting_permutation_test_%d_%d_Hz.png',results_dir,hp,lp));
end
fclose(resf);

%% Identifiability plot

f1 = figure();
f1.Position = [500,500,1100,300];
ax1 = axes;
hold on

xlabel('Frequency (Hz)')
set(f1,'Color','w')
set(ax1,'FontSize',13,'TickLength',[0,0],'XTickLabelRotation',70,'XTick',unique([hpfs,lpfs]));
[~,sortid]= sort(pval);
crit_pval_benjamini = 0.05./(20:-1:1);
mincolor = 0.3;
colorstep = 0.015;
for f_ind = 1:length(hpfs)
    
center = (lpfs(f_ind) + hpfs(f_ind))./2;
width = (lpfs(f_ind) - hpfs(f_ind))*0.975;
val = all_real_diffs(f_ind);
bh1(f_ind) = bar(ax1,center,val,width,'FaceColor',frq_cols(f_ind,:),'FaceAlpha',0.8);

if pval(f_ind) < crit_pval_benjamini(sortid(f_ind))
    th(f_ind) = text(ax1, center, val + 0.06, '*','FontSize',25,'HorizontalAlignment','center','VerticalAlignment','top');
    bh1(f_ind).FaceAlpha = 1;
end
end
% bh1(5).Visible = 'off';
% bh2(5).Visible = 'off';
% th(5).Visible = 'off';
ax1.Position(1) = 0.1;ax2.Position(1) = 0.1;
ax1.Position(3) = 0.85;ax2.Position(3) = 0.85;
ax1.YLabel.String = sprintf('Identifiability');
ax1.YLabel.FontSize = 13;
saveas(gcf,sprintf('%sFingerprinting_across_freqs.svg',results_dir));
saveas(gcf,sprintf('%sFingerprinting_across_freqs.png',results_dir));
save(sprintf('%sFingerprinting_across_freqs.mat',results_dir),'hpfs','lpfs','pval','all_real_diffs','all_mean_within_corrs');