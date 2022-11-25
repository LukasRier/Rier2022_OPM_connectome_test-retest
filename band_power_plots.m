if~exist('RESTOREDEFAULTPATH_EXECUTED');
    restoredefaultpath;
end
clearvars
clc
addpath ./BrainPlots/
addpath ./gifti-1.8/
project_dir = '/net/cador/data_local/Lukas/movie/';
datadir = [project_dir,'data',filesep];
exp_type = 'task-movie';
epoch_length = 5;

hpfs = [4, 8,13,30, 50,30,35,40,52,55,60,65,70,75,80,85, 90,95,102,108];
lpfs = [8,12,30,50,100,40,45,48,60,65,70,75,80,85,90,95,100,98,110,112];

hpfs = [4, 8,13,30,35,40,52,55,60,65,70,75,80,85, 90,95,102,108];
lpfs = [8,12,30,40,45,48,60,65,70,75,80,85,90,95,100,98,110,112];

if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_subs = 10;
n_ses = 2;
vars = nan(78,length(hpfs),n_subs*n_ses);
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
        load([path.data,filename,'_meg.mat'],'fs')
        
        
        path.VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
        if ~exist(path.VEs,'dir'); error("Run VEs first");end
        files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
        for f_i = 1:length(hpfs)
            hp = hpfs(f_i);
            lp = lpfs(f_i);
            fname = sprintf('%s%s_%d_%d_Hz_Z.mat',path.VEs,files.VEs,hp,lp);
            if exist(fname,'file')
                load(fname,'VE')
                VE = VE';
                
                % hilbert env
                
                %
                VE_mat = reshape(VE,[size(VE,1),epoch_length*fs,(round(size(VE,2)/fs))/epoch_length]);
                vars(:,f_i,sub_i+n_subs*(ses_i-1)) = mean(squeeze(var(VE_mat,[],2)),2);
                fprintf([repmat('\b',1,41),'Sub: %s | Ses: %s | %3d-%3d Hz (%2d/%2d)\n'],sub,ses,hp,lp,f_i,length(hpfs));
            end
        end
        
    end
end
%%
drawArrow = @(ax,x,y) quiver(ax, x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'k',...
    'ShowArrowHead','off');

mean_pow = mean(vars,3); % across subjects
fcents = (hpfs + lpfs)/2;fcents = fcents./(max(fcents));
picpos = linspace(0.01,1,length(hpfs)+1);picpos(end)=[];
[~,inds] = sort(fcents); picpos = picpos(inds);
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
    range = [];
    % range = [min(mean_pow(:)),max(mean_pow(:))];
    ax(f_i) = axes;
    PaintBrodmannAreas_1view(mean_pow(:,f_i),78, 256, range,[], [] )
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
suptitle('Movie experiment | Mean variance of narrow-band pseudo-Z VE across all recordings')

%% second version using broadband data and getting relative power
clearvars
addpath ./BrainPlots/
addpath ./gifti-1.8/
project_dir = '/net/cador/data_local/Lukas/movie/';
datadir = [project_dir,'data',filesep];
exp_type = 'task-movie';
epoch_length = 5;
%
% hpfs = [4, 8,13,30, 50,30,35,40,52,55,60,65,70,75,80,85, 90,95,102,108];
% lpfs = [8,12,30,50,100,40,45,48,60,65,70,75,80,85,90,95,100,98,110,112];

hpfs = [4, 8,13,30,35,40,52,55,60,65,70,75,80,85];
lpfs = [8,12,30,40,45,48,60,65,70,75,80,85,90,95];
%
% hpfs = [4, 8,13,30,35,40];
% lpfs = [8,12,30,40,45,48];


if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end
n_subs = 10;
n_ses = 2;
vars_psd = nan(78,length(hpfs),n_subs*n_ses);
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
        load([path.data,filename,'_meg.mat'],'fs')
        
        
        path.VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
        if ~exist(path.VEs,'dir'); error("Run VEs first");end
        
        files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
        fname = sprintf('%s%s_%d_%d_Hz.mat',path.VEs,files.VEs,1,150);
        load(fname,'VE')
        VE = VE';VE = VE - mean(VE,2);VE = VE./std(VE,[],2);
        [PSDs,fxx] = pwelch(VE',fs*epoch_length,0,[],fs);
        for f_i = 1:length(hpfs)
            hp = hpfs(f_i);
            lp = lpfs(f_i);
            band_inds = fxx >= hp & fxx <= lp;
            vars_psd(:,f_i,sub_i+n_subs*(ses_i-1)) = trapz(fxx(band_inds),PSDs(band_inds,:));
            fprintf([repmat('\b',1,41),'Sub: %s | Ses: %s | %3d-%3d Hz (%2d/%2d)\n'],sub,ses,hp,lp,f_i,length(hpfs));
        end
        
    end
end
%%
drawArrow = @(ax,x,y) quiver(ax, x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'k',...
    'ShowArrowHead','off');

mean_pow = mean(vars_psd,3); % across subjects
fcents = (hpfs + lpfs)/2;fcents = fcents./(max(fcents));
picpos = linspace(0.01,1,length(hpfs)+1);picpos(end)=[];
[~,inds] = sort(fcents); picpos = picpos(inds);
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
for f_i = 1:6;%length(hpfs)
    range = [];
    % range = [min(mean_pow(:)),max(mean_pow(:))];
    ax(f_i) = axes;
    PaintBrodmannAreas_1view(mean_pow(:,f_i),78, 256, range,[], [] )
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
suptitle('Movie experiment | Mean relative PSD across all recordings')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Method #3
% filt 4-150, beamform including sqrt(w^Tw) normalisation
% filt VE to bands of interest and look at variance and std dev of hilb env
% over mean of hilb env


clearvars
addpath ./BrainPlots/
addpath ./gifti-1.8/
project_dir = '/net/cador/data_local/Lukas/movie/';
datadir = [project_dir,'data',filesep];
exp_type = 'task-movie';
epoch_length = 5;
%
% hpfs = [4, 8,13,30, 50,30,35,40,52,55,60,65,70,75,80,85, 90,95,102,108];
% lpfs = [8,12,30,50,100,40,45,48,60,65,70,75,80,85,90,95,100,98,110,112];

hpfs = [4, 8,13,30,35,40,52,55,60,65,70,75,80,85];
lpfs = [8,12,30,40,45,48,60,65,70,75,80,85,90,95];


hpfs = [1,4, 8,13,30,35,40,52,55,60,65,70,75,80,85, 90,95,102,105,...
    110,115,120,125,130,135,140];
lpfs = [4,8,12,30,40,45,48,60,65,70,75,80,85,90,95,100,98,110,115,...
    120,125,130,135,140,145,150];
%
% hpfs = [4, 8,13,30,35,40];
% lpfs = [8,12,30,40,45,48];


if length(hpfs) ~= length(lpfs);error("Freq bands don't match");end
if ~all((hpfs-lpfs)<0);error("Not all hipass smaller than lowpass");end

n_subs = 10;
n_ses = 2;
band_vars = nan(78,length(hpfs),n_subs*n_ses);
env_fracs = nan(78,length(hpfs),n_subs*n_ses);
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
        path.VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
        
        files.pow = ['sub-',sub,'_','run-',ses,'_POW'];
        files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];

        load([path.data,filename,'_meg.mat'],'fs');
        for f_i = 1:length(hpfs)
            hp = hpfs(f_i);
            lp = lpfs(f_i);
            fname = sprintf('%s%s_%d_%d_Hz_Z.mat',path.pow,files.pow,hp,lp);
            load(fname,'filt_vars','ENV_std_over_mean');
            band_vars(:,f_i,sub_i+n_subs*(ses_i-1)) = filt_vars;
            env_fracs(:,f_i,sub_i+n_subs*(ses_i-1)) = ENV_std_over_mean;
            fprintf([repmat('\b',1,41),'Sub: %s | Ses: %s | %3d-%3d Hz (%2d/%2d)\n'],sub,ses,hp,lp,f_i,length(hpfs));
        end
    end
end
%
drawArrow = @(ax,x,y) quiver(ax, x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'k',...
    'ShowArrowHead','off');

mean_pow = mean(band_vars,3); % across subjects
mean_env_fracs = mean(env_fracs,3); % across subjects
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
    range = [];
    % range = [min(mean_pow(:)),max(mean_pow(:))];
    ax(f_i) = axes;
    PaintBrodmannAreas_1view(mean_pow(:,f_i),78, 256, range,[], [] )
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
    range = [];
    % range = [min(mean_pow(:)),max(mean_pow(:))];
    ax(f_i) = axes;
    PaintBrodmannAreas_1view(mean_env_fracs(:,f_i),78, 256, range,[], [] )
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
suptitle('Movie experiment | Mean \sigma_{env}/\mu_{env} fraction in bands across all recordings')

%%
close all
x=VE(:,24);
hp = 4;
lp = 100;
[b,a] = butter(4,2*[hp lp]/fs);
x = filtfilt(b,a,x);

figure
cdfplot((x - mean(x))./std(x))
hold on
x_values = linspace(min(x),max(x));
plot(x_values,normcdf(x_values,0,1),'r-')
legend('Empirical CDF','Standard Normal CDF','Location','best')
figure
histogram(x)
kstest((x-mean(x))/std(x))

h_x = abs(hilbert(x));
figure
subplot(2,2,[1,2])
plot(x,'k')
hold on
plot(h_x,'r')
xlim([0,5*600])

subplot(2,2,3)
histogram(x)
title('VE distribution (normal ish)')
subplot(2,2,4)
histogram(h_x)
title('Hilbert envelope distribution (Raleigh ish)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % Example spectra
hp = 4;
lp = 150;
n_subs = 10;
n_ses = 2;
PSDs_raw = [];
PSDs_std = [];
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
        load(sprintf('%s%s_%d_%d_Hz_Z_standard.mat',path.VEs,files.VEs,hp,lp),'VE')
        VE_std = VE;
        
        [psd,fxx] = pwelch(VE_raw(:,16),fs*5,[],[],fs);
        [psd_std,~] = pwelch(VE_std(:,16),fs*5,[],[],fs);
        
        PSDs_raw = cat(2,PSDs_raw,psd);
        PSDs_std = cat(2,PSDs_std,psd_std);
            
        figure
        subplot(2,1,1)
        plot(fxx,psd);
        xlim([0,50])
        subplot(2,1,2)
        plot(fxx,psd_std);
        xlim([0,50])
        suptitle(['sub ',sub,' session ',ses])
        drawnow
    end
end

mean_PSD_raw = mean(PSDs_raw,2);
mean_PSD_std = mean(PSDs_std,2);
figure
subplot(2,1,1)
plot(fxx,mean_PSD_raw);
xlim([0,50])
subplot(2,1,2)
plot(fxx,mean_PSD_std);
xlim([0,50])
drawnow
regs = [16,25];

%%
xlims = [0,50];
for x = regs
    figure;set(gcf,'Color','w');
    subplot(2,3,1)
    rrr=nan(1,78);rrr(x) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end
    subplot(2,3,[2,3])
    plot(fxx,mean_PSD_raw(:,x));
    title("Raw psd")    
    xlim(xlims)


    subplot(2,3,4)
    rrr=nan(1,78);rrr(x+39) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
        if x == 25;view([0,0]),camlight;end

    subplot(2,3,[5,6])
    plot(fxx,mean_PSD_raw(:,x+39));
    xlim(xlims)

    %%%%%%%%%%
    
    figure;set(gcf,'Color','w');
    subplot(2,3,1)
    rrr=nan(1,78);rrr(x) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end

    subplot(2,3,[2,3])
    plot(fxx,mean_PSD_std(:,x));
    title("Standardised psd")
    xlim(xlims)

    subplot(2,3,4)
    rrr=nan(1,78);rrr(x+39) = 1;
    PaintBrodmannAreas_1view(rrr,78, 256, [0,1],[], [] )
    if x == 25;view([0,0]),camlight;end

    subplot(2,3,[5,6])
    plot(fxx,mean_PSD_std(:,x+39));
    xlim(xlims)

end

%% Run1 run2 comparison


mean_PSD_raw = mean(PSDs_raw,3);
mean_PSD_std = mean(PSDs_std,3);

[psd,fx] = pwelch(VE_raw,fs*5,[],[],fs);
[psd_std,~] = pwelch(VE_std,fs*5,[],[],fs);
%%

figure
subplot(2,1,1)
plot(fx,psd(:,16));
xlim([0,50])
subplot(2,1,2)
plot(fx,psd_std(:,16));
xlim([0,50])