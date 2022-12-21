function [xn, yn, cn, xi, yi,p] = plot_fieldmaps_LR(vals, lay, ax)
%% PLOT_FIELDMAPS_LR
% Generate a fieldmap for OPM-MEG data
%
% vals...values to be plotted (one per channel in lay
% lay....layout in FieldTrip format containing fields pos with x/y coordinates and
% outline of the head
% ax.....Figure axes handle

layout = lay;

% define colormap
ncolor = 100;
cmap =[linspace(1,0,ncolor/4);linspace(1,0,ncolor/4); linspace(1,1,ncolor/4)]';
cmap =[cmap;[linspace(0,0,ncolor/4); linspace(0,0,ncolor/4);linspace(1,0,ncolor/4)]'];
cmap =[cmap;[linspace(0,1,ncolor/4);linspace(0,0,ncolor/4); linspace(0,0,ncolor/4)]'];
cmap =[cmap;[linspace(1,1,ncolor/4); linspace(0,1,ncolor/4);linspace(0,1,ncolor/4)]'];
% white: maximal field/potential
% red:   postive field/potential
% black: zero field/potential
% blue:  negative field/potential
% white: minimal field/potential

c2 = vals;
z_idx = find(isnan(c2)); 
c2(z_idx) = [];
layout.pos(z_idx,:) = [];

% create grid

ngrid  = 500;
x = layout.pos(:,1)';
y = layout.pos(:,2)';

xmin = min(layout.outline{1,1}(:,1));
xmax = max(layout.outline{1,1}(:,1));
ymin = min(layout.outline{1,1}(:,2));
ymax = max(layout.outline{1,1}(:,2));
[xi, yi] = meshgrid(linspace(xmin, xmax, ngrid), linspace(ymin, ymax, ngrid));

if size(c2,1) ~= length(x)
    disp('ERROR: size(cData) should be Nx1 (where N = number of channels)')
    return
end

cn = [c2',zeros(1,size(layout.outline{1,1},1))];
xn = [x,layout.outline{1,1}(:,1)'];
yn = [y,layout.outline{1,1}(:,2)'];
[xi, yi, ci] = griddata(xn, yn, cn, xi, yi, 'cubic');


% plot
% figure
ax;
colormap(cmap)
p = pcolor(ax,xi,yi,ci);
p.EdgeColor = 'none';
hold on
ft_plot_lay(layout,'label','no','box','no')
% cb = colorbar;
% caxis([-max(abs(c)),max(abs(c))])
% fig = gcf;fig.Color = 'w';