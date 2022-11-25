function PaintBrodmannAreas_1view(data, ...
    nr_ROIs_to_plot, nr_colors_to_plot, colour_range, cbar_label, colourbar_threshold)
% get labels from Arjans selection
mesh_type=[];
load('ROI_MNI_V4_List.mat','ROI');
ROI_indices = select_ROIs_from_full_AAL;
labels_temp = char(ROI.Nom_L);
labels = labels_temp(ROI_indices,:);

addpath(genpath('/home/ppzpkt/matlab/spm12'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LR hack for thresholding
if isempty(colour_range)
    colour_range = minmax(data');
end
data(data>colour_range(2)) = colour_range(2);
data(data<colour_range(1)) = colour_range(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%
% Plots data on cortical surface mesh from SPM (for AAL atlas).
% The data en labels should be in the same order (i.e. data(i) corresponds to labels{i})

% make all the areas grey
%braincolor = [0 0 0]; % black
braincolor = [0.6 0.6 0.6];% grey

lighting_type = 'phong'; % very shiny - metalic like
%lighting_type = 'gouraud'; % a bit shiny
%lighting_type = 'flat';

%lcolor=[219,112,147]/255; % pink light
lcolor='w'; % white light
%brain_material([.55,.6,.4,10]);
brain_material = 'dull';
%brain_material = 'shiny';

patch_handles = [];

if nargin<4 | isempty(colour_range)
    colour_range = [min(data), max(data)]; % set the colour_range based on the provided data
    if colour_range(1)== colour_range(2)
        colour_range(1) = -1*colour_range(1); colour_range = sort(colour_range);
        disp(sprintf('setting colour_range to [%5.2f %5.2f]',colour_range))
    end
elseif length(colour_range)~=2
    error('Please provide a maximum and a minimum value for the colour range')
end
if nargin<6
    colourbar_threshold=[]; % for thresholding of the images
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the mesh and the labels for the lh mesh
tmp = load('cortex_20484.surf.gii_labelled_lh.mat');
% %             remove 2 lines below for Helen:
test = gifti(tmp.lhmesh); % had to do this in order to get compiled version working
clear test
%
tmp.positions = tmp.lhmesh.vertices;
tmp.polygons = tmp.lhmesh.faces;

positions=tmp.positions;
polygons=tmp.polygons;
approximate_meshlabels=char(tmp.approximate_meshlabels);
clear tmp

colormap_index3_lh=[];
% make all the areas grey
% colormap_index3_lh(1:size(positions,1),:) = 0.*ones(size(positions,1),3);
colormap_index3_lh(1:size(positions,1),:) = repmat(braincolor, size(positions,1),1);

map=colormap(turbo(nr_colors_to_plot));
% map=colormap(parula(nr_colors_to_plot));
%map=colormap(winter(nr_colors_to_plot));
% map=colormap(summer(nr_colors_to_plot));
%map=colormap(bone(nr_colors_to_plot));
map=colormap(flipud(gray(nr_colors_to_plot)));
%map=colormap(autumn(nr_colors_to_plot));
%my_autumn = flipud(colormap(autumn(nr_colors_to_plot)));
%map=colormap(my_autumn);
%%map=flipud(map);
colormap(map);

%% Main figure, left hemisphere, top view
% fh1 = figure(gcf);
% fh1.Position = [26 548 560 420];
% subplot(2,2,[1 3])
patches=patch('faces', polygons, 'vertices', positions);
set(patches, 'CDataMapping', 'direct' );
set(patches, 'FaceColor', 'interp');
set(patches, 'EdgeColor', 'none');
set(patches, 'FaceLighting',lighting_type);

[colormap_index3_lh] = local_colourin_ROIs(colormap_index3_lh, ...
    approximate_meshlabels, labels, data, nr_ROIs_to_plot, nr_colors_to_plot, ...
    map, colour_range, braincolor, mesh_type);

% color the patches
set(patches, 'FaceVertexCData',colormap_index3_lh);

% change appearance
caxis([min(colour_range) max(colour_range)]);
% AH, March 2013 -- repeat the colormap statement here in order to update the colorbar
colormap(map);

axis off;
%legh=legend(alllh, legend_labels_lh,'Location','BestOutside');
% set(gcf,'Position', [0         165        1225         780]);
% axis equal
lh=camlight;
set(lh, 'Color', lcolor)
material(brain_material)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the mesh and colours for the right hemisphere
% load the mesh and the labels for the mesh

tmp = load('cortex_20484.surf.gii_labelled_rh.mat');

% remove 2 lines below for Helen:
test = gifti(tmp.rhmesh); % had to do this in order to get compiled version working
clear test

tmp.positions = tmp.rhmesh.vertices;
tmp.polygons = tmp.rhmesh.faces;

positions_rh=tmp.positions;
polygons_rh=tmp.polygons;
approximate_meshlabels_rh=char(tmp.approximate_meshlabels);
clear tmp

colormap_index3_rh=[];
% make all the areas grey
%colormap_index3_rh(1:size(positions_rh,1),:) = 0.*ones(size(positions_rh,1),3);
colormap_index3_rh(1:size(positions_rh,1),:) = repmat(braincolor, size(positions_rh,1),1);

% find the colours for the patches
[colormap_index3_rh] = local_colourin_ROIs(colormap_index3_rh, ...
    approximate_meshlabels_rh, labels, data, nr_ROIs_to_plot, ...
    nr_colors_to_plot, map, colour_range, braincolor, mesh_type);


%% Main figure - right hemisphere; top view
% subplot(2,2,[1 3])
patches6=patch('faces',polygons_rh, 'vertices', positions_rh);
set(patches6, 'CDataMapping', 'direct');
set(patches6, 'FaceColor', 'interp');
set(patches6, 'EdgeColor', 'none');
set(patches6, 'FaceVertexCData',colormap_index3_rh);
set(patches6, 'FaceLighting',lighting_type);
material(brain_material)
% whitebg(gcf,'k')
% set(gcf,'color','k')
% cbh = colorbar;
% set(cbh,'Visible','on', 'location', 'southoutside');
% cbh.Label.String = cbar_label;
%%
if ~isempty(colourbar_threshold)
    colorindex1 = 1 + fix((colourbar_threshold(1)-colour_range(1))/(colour_range(2)-colour_range(1))*(nr_colors_to_plot-1));
    colorindex2 = 1 + fix((colourbar_threshold(2)-colour_range(1))/(colour_range(2)-colour_range(1))*(nr_colors_to_plot-1));
    map(colorindex1:colorindex2,:)=repmat([0 0 0],colorindex2-colorindex1+1,1);
    %     figure(123)
    colormap(map)
    %     figure(124)
    %     colormap(map)
end

function [colormap_index3] = local_colourin_ROIs(colormap_index3, ...
    approximate_meshlabels, labels, data, nr_ROIs_to_plot, ...
    nr_colors_to_plot, map, colour_range, braincolor, mesh_type)

for j=1:nr_ROIs_to_plot
    searchlabel = deblank(char(labels(j,:)));
    colorindex = 1 + fix((data(j)-colour_range(1))/(colour_range(end)-colour_range(1))*(nr_colors_to_plot-1));
    if ~isnan(colorindex)
        surfacecolor = map(colorindex,:);
    else
        % color is invisible
        %surfacecolor = [0 0 0];
        surfacecolor = braincolor;
    end
    labelind = strmatch(searchlabel,approximate_meshlabels(:,1:length(searchlabel)), 'exact');
    if ~isempty(labelind) % color all relevant triangles
        for i=1:length(labelind)
            colormap_index3(labelind(i),:) = surfacecolor;
        end
    end
end


