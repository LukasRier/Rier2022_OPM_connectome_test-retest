if~exist('RESTOREDEFAULTPATH_EXECUTED');
    restoredefaultpath;
end
clearvars
clc
addpath ./BrainPlots/
addpath ./gifti-1.8/
project_dir = '/net/cador/data_local/Lukas/movie/';
results_dir = [project_dir,'results',filesep,'average_spectra',filesep];
if ~exist(results_dir,'dir'); mkdir(results_dir);end
datadir = [project_dir,'data',filesep];
exp_type = 'task-movie';
epoch_length = 5;

%
% hpfs = [4, 8,13,30, 50,30,35,40,52,55,60,65,70,75,80,85, 90,95,102,108];
% lpfs = [8,12,30,50,100,40,45,48,60,65,70,75,80,85,90,95,100,98,110,112];

hpfs = [4, 8,13,30,35,40];
lpfs = [8,12,30,40,45,48];



if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end

n_subs = 10;
n_ses = 2;
band_vars = nan(78,length(hpfs),n_subs*n_ses);
for sub_i = 1:n_subs
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0';
    for ses_i = 1:n_ses
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0';
        % sub = '001'
        % ses = '001'
        fprintf(repmat(' ',1,41))
        
        filename = ['sub-',sub,'_',exp_type,'_','run-',ses];
        path.main = [datadir,'sub-',sub,filesep];
        path.data = [path.main,'meg',filesep];
        path.pow = [datadir,'derivatives',filesep,'POW',filesep,'sub-',sub,filesep];
        
        files.pow = ['sub-',sub,'_','run-',ses,'_POW'];
        
        load([path.data,filename,'_meg.mat'],'fs');
        for f_i = 1:length(hpfs)
            hp = hpfs(f_i);
            lp = lpfs(f_i);
            fname = sprintf('%s%s_%d_%d_Hz_Z_standard.mat',path.pow,files.pow,hp,lp);
            load(fname,'filt_vars','ENV_std_over_mean');
            band_vars(:,f_i,sub_i+n_subs*(ses_i-1)) = filt_vars;
            fprintf([repmat('\b',1,41),'Sub: %s | Ses: %s | %3d-%3d Hz (%2d/%2d)\n'],sub,ses,hp,lp,f_i,length(hpfs));
        end
    end
end
%
drawArrow = @(ax,x,y) quiver(ax, x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'k',...
    'ShowArrowHead','off');

mean_pow = mean(band_vars,3); % across subjects
fcents = (hpfs + lpfs)/2;fcents = fcents./(max(fcents));
picpos = linspace(0.01,1,length(hpfs)+1);picpos(end)=[];
[~,inds] = sort(fcents); picpos = picpos(inds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variance in filter band
f = figure;
set(f,'Color','w','WindowState','Maximized','Position',get(0, 'Screensize'))
% set(f,'outerposition',[0,0,1,1]);
% set(f,'Position',[0   -0.4583    1.5000    1.3500]);
faxis = axes;
faxis.Position = [0.005,0.4,0.98,0.2];
faxis.XLim=[0,max(lpfs)];faxis.YLim=[0,.5];
faxis.XLim=[-0.05*max(lpfs),max(lpfs)*1.1];faxis.YLim=[0,.5];
faxis.NextPlot = 'add';faxis.YAxis.Visible='off';
set(faxis,'FontSize',8,'TickLength',[0,0],'XTickLabelRotation',70,'XTick',unique([hpfs,lpfs]));
% mincolor = 0.3;
% colorstep = 0.015;
ax_wdth = faxis.XLim(2) - faxis.XLim(1);
xlabel(faxis,'Frequency (Hz)');

% faxis.YAxis.Visible='off';
for f_i = 1:length(hpfs)
    col_range = [0.5,1].*max(mean_pow(:,f_i));
    %         col_range = [0,max(mean_pow(:,f_i))];
    %     col_range = [min(mean_pow(:)),max(mean_pow(:))];
    ax(f_i) = axes;
    PaintBrodmannAreas_1view(mean_pow(:,f_i),78, 256, col_range,[], [] )
    ax(f_i).Position = [picpos(f_i),0.6,mean(diff(picpos)),1.5*mean(diff(picpos))];
    cb(f_i) = colorbar;set(cb(f_i),'Location','southoutside');
    cb(f_i).Label.FontSize = 6;cb(f_i).FontSize=6;
    ax(f_i).Position = [picpos(f_i),0.6,mean(diff(picpos)),1.5*mean(diff(picpos))];
    cb(f_i).Position([2,3]) = [cb(f_i).Position(2)+0.002,cb(f_i).Position(3)-0.004];
    
    center = (lpfs(f_i) + hpfs(f_i))./2;
    width = (lpfs(f_i) - hpfs(f_i))*0.975;
    val = 0.15+0.01*mod(f_i,2);
    %     bh1(f_i) = bar(faxis,center,val,width,'FaceColor',[mincolor+colorstep*f_i,mincolor+colorstep*f_i,1],'FaceAlpha',0.8);
    bh1(f_i) = bar(faxis,center,val,width,'FaceColor',[1,1,1],'FaceAlpha',0.4);
    
    drawArrow(faxis,[center,...
        -((ax_wdth./faxis.Position(3))*faxis.Position(1)-faxis.XLim(1)) + (ax_wdth./faxis.Position(3))*(picpos(f_i)+mean(diff(picpos))./2)],...
        [val,.6*faxis.YLim(2)])
    drawnow
    % pause
end
suptitle('Movie experiment | Mean variance in bands across all recordings')

%%
band_name = {'\theta','\alpha','\beta','\gamma_1','\gamma_2','\gamma_3'};

for f_i = 1:length(hpfs)
    %     col_range = [0,max(mean_pow(:))];
    %     col_range = [min(mean_pow(:)),max(mean_pow(:))];
    col_range = [0,1].*max(mean_pow(:,f_i));
    f1 = figure;
    set(f1,'Color','w','Position',[680    71   386   907])
    ax(1) = axes;
    PaintBrodmannAreas_chooseview(mean_pow(:,f_i),78, 256, col_range,[], [], [-90,0])
    ax(1).Position = [0.01,0.1,0.9,0.3];
    cb = colorbar;set(cb,'Location','southoutside');
    cb.Label.FontSize = 6;cb.FontSize=15;
    cb.Position([1,2]) = cb.Position([1,2])-[-0.05,0.07];
    %     cb.Position([2,3]) = [cb.Position(2)-0.002,cb.Position(3)-0.004];
    drawnow
    
    ax(2) = axes;
    PaintBrodmannAreas_chooseview(mean_pow(:,f_i),78, 256, col_range,[], [] )
    ax(2).Position = [0.08,0.70,0.9,0.3];
    drawnow
    
    ax(3) = axes;
    PaintBrodmannAreas_chooseview(mean_pow(:,f_i),78, 256, col_range,[], [], [90,0] )
    ax(3).Position = [0.1,0.38,0.9,0.3];
    
    st = suptitle(band_name{f_i})
    % pause
    st.FontSize = 35;
    drawnow
    saveas(f1,sprintf('%s%s_average_brains.png',results_dir,band_name{f_i}(2:end)))
end


%% % Example spectra
hp = 1;
lp = 150;
n_subs = 10;
n_ses = 2;
PSDs_std = [];PSDs_raw = [];noisePSDs_std = [];noisePSDs_raw = [];

for sub_i = 1:n_subs
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:n_ses
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        
        path.VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
        files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
        filename = ['sub-',sub,'_',exp_type,'_','run-',ses];
        path.main = [datadir,'sub-',sub,filesep];
        path.data = [path.main,'meg',filesep];
        load([path.data,filename,'_meg.mat'],'fs')
        load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.VEs,files.VEs,hp,lp),'VE')
        VE_raw = VE;
        [psd_raw,fxx] = pwelch(VE_raw,fs*5,[],[],fs);
        PSDs_raw = cat(3,PSDs_raw,psd_raw);
        load(sprintf('%s%s_%d_%d_Hz_Z_standard.mat',path.VEs,files.VEs,hp,lp),'VE')
        VE_std = VE;
        [psd_std,fxx] = pwelch(VE_std,fs*5,[],[],fs);
        PSDs_std = cat(3,PSDs_std,psd_std);
        
        
        path.noiseVEs = [datadir,'derivatives',filesep,'noiseVEs',filesep,'sub-',sub,filesep];
        files.noiseVEs = ['sub-',sub,'_','run-',ses,'_noiseVE'];
        load(sprintf('%s%s_%d_%d_Hz_Z.mat',path.noiseVEs,files.noiseVEs,hp,lp),'VE_noise')
        VE_raw = VE_noise;
        [psd_raw,fxx] = pwelch(VE_raw,fs*5,[],[],fs);
        noisePSDs_raw = cat(3,noisePSDs_raw,psd_raw);
        load(sprintf('%s%s_%d_%d_Hz_Z_standard.mat',path.noiseVEs,files.noiseVEs,hp,lp),'VE_noise')
        VE_std = VE_noise;
        [psd_std,fxx] = pwelch(VE_std,fs*5,[],[],fs);
        noisePSDs_std = cat(3,noisePSDs_std,psd_std);
    end
end

mean_PSD_std = mean(PSDs_std(:,:,1:2:n_subs*n_ses),3);
mean_PSD_std2 = mean(PSDs_std(:,:,2:2:n_subs*n_ses),3);
err_PSD_std =  std(PSDs_std(:,:,1:2:n_subs*n_ses),[],3)./sqrt(n_subs);
err_PSD_std2 = std(PSDs_std(:,:,2:2:n_subs*n_ses),[],3)./sqrt(n_subs);

mean_PSD_raw = mean(PSDs_raw(:,:,1:2:n_subs*n_ses),3);
mean_PSD_raw2 = mean(PSDs_raw(:,:,2:2:n_subs*n_ses),3);
err_PSD_raw =  std(PSDs_raw(:,:,1:2:n_subs*n_ses),[],3)./sqrt(n_subs);
err_PSD_raw2 = std(PSDs_raw(:,:,2:2:n_subs*n_ses),[],3)./sqrt(n_subs);

mean_noisePSD_std = mean(noisePSDs_std(:,:,1:2:n_subs*n_ses),3);
mean_noisePSD_std2 = mean(noisePSDs_std(:,:,2:2:n_subs*n_ses),3);
err_noisePSD_std =  std(noisePSDs_std(:,:,1:2:n_subs*n_ses),[],3)./sqrt(n_subs);
err_noisePSD_std2 = std(noisePSDs_std(:,:,2:2:n_subs*n_ses),[],3)./sqrt(n_subs);

mean_noisePSD_raw = mean(noisePSDs_raw(:,:,1:2:n_subs*n_ses),3);
mean_noisePSD_raw2 = mean(noisePSDs_raw(:,:,2:2:n_subs*n_ses),3);
err_noisePSD_raw =  std(noisePSDs_raw(:,:,1:2:n_subs*n_ses),[],3)./sqrt(n_subs);
err_noisePSD_raw2 = std(noisePSDs_raw(:,:,2:2:n_subs*n_ses),[],3)./sqrt(n_subs);
%
% mean_PSD_raw = sqrt(mean(PSDs_raw(:,:,1:2:n_subs*n_ses),3));
% mean_PSD_raw2 = sqrt(mean(PSDs_raw(:,:,2:2:n_subs*n_ses),3));
% err_PSD_raw =  sqrt(std(PSDs_raw(:,:,1:2:n_subs*n_ses),[],3)./sqrt(n_subs));
% err_PSD_raw2 = sqrt(std(PSDs_raw(:,:,2:2:n_subs*n_ses),[],3)./sqrt(n_subs));
% mean_noisePSD_raw = sqrt(mean(noisePSDs_raw(:,:,1:2:n_subs*n_ses),3));
% mean_noisePSD_raw2 = sqrt(mean(noisePSDs_raw(:,:,2:2:n_subs*n_ses),3));
% err_noisePSD_raw =  sqrt(std(noisePSDs_raw(:,:,1:2:n_subs*n_ses),[],3)./sqrt(n_subs));
% err_noisePSD_raw2 = sqrt(std(noisePSDs_raw(:,:,2:2:n_subs*n_ses),[],3)./sqrt(n_subs));


%%

regs = [16,25,7];
xlims = [0,50];
c1 = [0.1216    0.4706    0.7059];
c2 = [  0.2000    0.6275    0.1725];
for x = regs
    figure;set(gcf,'Color','w','Position',[63   507   815   420]);
    s1 = subplot(2,1,1);
    plot(fxx,mean_PSD_raw(:,x),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_raw(:,x)-err_PSD_raw(:,x),mean_PSD_raw(:,x)+err_PSD_raw(:,x),fxx,c1)
    
    plot(fxx,mean_PSD_raw2(:,x),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_raw2(:,x)-err_PSD_raw2(:,x),mean_PSD_raw2(:,x)+err_PSD_raw2(:,x),fxx,c2)
    
    xlim(xlims)
    xlabel('Frequency (Hz)')
    legend('Run1','std. err.','Run2','std. err.','Location','southeast')
%     legend('Run1','Run2','Location','southeast')
    ylabel('A.U.')
    
    s2 = subplot(2,1,2);
    plot(fxx,mean_PSD_raw(:,x+39),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_raw(:,x+39)-err_PSD_raw(:,x+39),mean_PSD_raw(:,x+39)+err_PSD_raw(:,x+39),fxx,c1)
    
    plot(fxx,mean_PSD_raw2(:,x+39),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_raw2(:,x+39)-err_PSD_raw2(:,x+39),mean_PSD_raw2(:,x+39)+err_PSD_raw2(:,x+39),fxx,c2)

    xlim(xlims)
    xlabel('Frequency (Hz)')
    ylabel('A.U.')
    
    ax = axes;
    rrr=.1.*ones(1,78);rrr(x) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax.Position=[0.8,0.75,0.1,0.15];
    
    ax2=axes;
    rrr=.1.*ones(1,78);rrr(x+39) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax2.Position=[0.8,0.27,0.1,0.15];
    saveas(gcf,sprintf('%sregion_%d_average_specs_nonorm.png',results_dir,x))

end

%
regs = [16,25,7];
xlims = [0,50];
c1 = [0.1216    0.4706    0.7059];
c2 = [  0.2000    0.6275    0.1725];
for x = regs
    figure;set(gcf,'Color','w','Position',[63   507   815   420]);
    s1 = subplot(2,1,1);
    plot(fxx,mean_PSD_std(:,x),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_std(:,x)-err_PSD_std(:,x),mean_PSD_std(:,x)+err_PSD_std(:,x),fxx,c1)
    
    plot(fxx,mean_PSD_std2(:,x),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_std2(:,x)-err_PSD_std2(:,x),mean_PSD_std2(:,x)+err_PSD_std2(:,x),fxx,c2)
    
    xlim(xlims)
    xlabel('Frequency (Hz)')
    legend('Run1','std. err.','Run2','std. err.','Location','southeast')
%     legend('Run1','Run2','Location','southeast')
    ylabel('A.U.')
    
    s2 = subplot(2,1,2);
    plot(fxx,mean_PSD_std(:,x+39),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_std(:,x+39)-err_PSD_std(:,x+39),mean_PSD_std(:,x+39)+err_PSD_std(:,x+39),fxx,c1)
    
    plot(fxx,mean_PSD_std2(:,x+39),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_std2(:,x+39)-err_PSD_std2(:,x+39),mean_PSD_std2(:,x+39)+err_PSD_std2(:,x+39),fxx,c2)

    xlim(xlims)
    xlabel('Frequency (Hz)')
    ylabel('A.U.')
    
    ax = axes;
    rrr=.1.*ones(1,78);rrr(x) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax.Position=[0.8,0.75,0.1,0.15];
    
    ax2=axes;
    rrr=.1.*ones(1,78);rrr(x+39) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax2.Position=[0.8,0.27,0.1,0.15];
    saveas(gcf,sprintf('%sregion_%d_average_specs.png',results_dir,x))

end

%% with noise

regs = [16,25,7];
xlims = [0,50];
c1 = [0.1216    0.4706    0.7059];
c2 = [  0.2000    0.6275    0.1725];
for x = regs
    figure;set(gcf,'Color','w','Position',[63   507   815   420]);
    s1 = subplot(2,1,1);
    plot(fxx,mean_PSD_std(:,x),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_std(:,x)-err_PSD_std(:,x),mean_PSD_std(:,x)+err_PSD_std(:,x),fxx,c1)
    
    plot(fxx,mean_PSD_std2(:,x),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_std2(:,x)-err_PSD_std2(:,x),mean_PSD_std2(:,x)+err_PSD_std2(:,x),fxx,c2)
    
    plot(fxx,mean_noisePSD_std(:,x),'Color','r','Linewidth',2);
    ciplot(mean_noisePSD_std(:,x)-err_noisePSD_std(:,x),mean_noisePSD_std(:,x)+err_noisePSD_std(:,x),fxx,[1,0,0])
    
    
    xlim(xlims)
    xlabel('Frequency (Hz)')
    legend('Run1','std. err.','Run2','std. err.','Empty Room','std. err.','Location','north')
%     legend('Run1','Run2','Location','southeast')
    ylabel('A.U.')
    
    s2 = subplot(2,1,2);
    plot(fxx,mean_PSD_std(:,x+39),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_std(:,x+39)-err_PSD_std(:,x+39),mean_PSD_std(:,x+39)+err_PSD_std(:,x+39),fxx,c1)
    
    plot(fxx,mean_PSD_std2(:,x+39),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_std2(:,x+39)-err_PSD_std2(:,x+39),mean_PSD_std2(:,x+39)+err_PSD_std2(:,x+39),fxx,c2)
    
    plot(fxx,mean_noisePSD_std(:,x),'Color','r','Linewidth',2);
    ciplot(mean_noisePSD_std(:,x)-err_noisePSD_std(:,x),mean_noisePSD_std(:,x)+err_noisePSD_std(:,x),fxx,[1,0,0])
    
    xlim(xlims)
    xlabel('Frequency (Hz)')
    ylabel('A.U.')
    
    ax = axes;
    rrr=.1.*ones(1,78);rrr(x) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax.Position=[0.8,0.75,0.1,0.15];
    
    ax2=axes;
    rrr=.1.*ones(1,78);rrr(x+39) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax2.Position=[0.8,0.27,0.1,0.15];
    saveas(gcf,sprintf('%sregion_%d_average_specs_w_noise_std.png',results_dir,x))

end

%%

regs = [16,25,7];
xlims = [0,50];
c1 = [0.1216    0.4706    0.7059];
c2 = [  0.2000    0.6275    0.1725];
for x = regs
    figure;set(gcf,'Color','w','Position',[63   507   815   420]);
    s1 = subplot(2,1,1);
    plot(fxx,mean_PSD_raw(:,x),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_raw(:,x)-err_PSD_raw(:,x),mean_PSD_raw(:,x)+err_PSD_raw(:,x),fxx,c1)
    
    plot(fxx,mean_PSD_raw2(:,x),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_raw2(:,x)-err_PSD_raw2(:,x),mean_PSD_raw2(:,x)+err_PSD_raw2(:,x),fxx,c2)
    
    plot(fxx,mean_noisePSD_raw(:,x),'Color','r','Linewidth',2);
    ciplot(mean_noisePSD_raw(:,x)-err_noisePSD_raw(:,x),mean_noisePSD_raw(:,x)+err_noisePSD_raw(:,x),fxx,[1,0,0])
    
    
    xlim(xlims)
    xlabel('Frequency (Hz)')
    legend('Run1','std. err.','Run2','std. err.','Empty Room','std. err.','Location','north')
%     legend('Run1','Run2','Location','southeast')
    ylabel('fT^2')
    
    s2 = subplot(2,1,2);
    plot(fxx,mean_PSD_raw(:,x+39),'Color',c1,'Linewidth',2);
    hold on
    ciplot(mean_PSD_raw(:,x+39)-err_PSD_raw(:,x+39),mean_PSD_raw(:,x+39)+err_PSD_raw(:,x+39),fxx,c1)
    
    plot(fxx,mean_PSD_raw2(:,x+39),'Color',c2,'Linewidth',2);
    ciplot(mean_PSD_raw2(:,x+39)-err_PSD_raw2(:,x+39),mean_PSD_raw2(:,x+39)+err_PSD_raw2(:,x+39),fxx,c2)
    
    plot(fxx,mean_noisePSD_raw(:,x),'Color','r','Linewidth',2);
    ciplot(mean_noisePSD_raw(:,x)-err_noisePSD_raw(:,x),mean_noisePSD_raw(:,x)+err_noisePSD_raw(:,x),fxx,[1,0,0])
    
    xlim(xlims)
    xlabel('Frequency (Hz)')
    ylabel('fT^2')
    
    ax = axes;
    rrr=.1.*ones(1,78);rrr(x) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax.Position=[0.8,0.75,0.1,0.15];
    
    ax2=axes;
    rrr=.1.*ones(1,78);rrr(x+39) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    ax2.Position=[0.8,0.27,0.1,0.15];
    saveas(gcf,sprintf('%sregion_%d_average_specs_w_noise.png',results_dir,x))

end