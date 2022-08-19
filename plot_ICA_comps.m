function [bad_comps] = plot_ICA_comps(comp150,ch_table,lay,bad_comps)

%% Lukas Rier July 2022
% NOT QUITE DONE

comps_n_page = 7;
currpage = 1;
n_comps = size(comp150.topo,2);

comps_table = cell(n_comps,3);
table_color = zeros(n_comps,3);
for comp_i = 1:n_comps
    comps_table{comp_i,1} = comp_i;
    comps_table{comp_i,2} = 'Good';
    comps_table{comp_i,3} = true;
    table_color(comp_i,:) = [0.6 1 0.6];
    if ismember(comp_i,bad_comps)
        comps_table{comp_i,2} = 'Bad';
        comps_table{comp_i,3} = false;
        table_color(comp_i,:) = [1 0.6 0.6];
    end
end
gap = 0.12;
comps_good = ones(n_comps,1);comps_good(bad_comps)=0;
% n_pages = ceil(n_comps./comps_n_page);
curr_comps = 1 + mod(((currpage - 1)*comps_n_page) : currpage*comps_n_page-1,n_comps);
compdata = [comp150.trial{1,:}];

f = figure;
f.Units = 'Normalized';f.Position = [0,0,1,1];
f.Color = 'w';

for compi = 1:comps_n_page
    ax1(compi) = axes('Position',[0.05,0.85 - gap.*(compi-1),0.1,0.1]);
    % y comps
    ax2(compi) = axes('Position',[0.13,0.85 - gap.*(compi-1),0.1,0.1]);
    % z comps
    ax3(compi) = axes('Position',[0.21,0.85 - gap.*(compi-1),0.1,0.1]);
    % comp time course
    ax4(compi) = axes('Position',[0.35,0.85 - gap.*(compi-1),0.3,0.1]);
end


t = uitable('Data',comps_table,'ColumnName',{'Comp Number','Status','Toggle'},...
    'Units', 'Normalized', 'Position',[0.7, 0.1, 0.2, 0.8],'BackgroundColor',table_color);
set(t,'ColumnEditable',[false false true])

update_button = uicontrol('Style','PushButton','Units', 'Normalized',...
    'Position',[0.7 0.05 0.05 0.05],'String','Update','Callback',@update_bad_ch);
done_button = uicontrol('Style','PushButton','Units', 'Normalized',...
    'Position',[0.75 0.05 0.1 0.05],'String','Save bad components','Callback',@done_bad_comps);

next_page_button = uicontrol('Style','PushButton','Units', 'Normalized',...
    'Position',[0.01 0.1 0.05 0.05],'String','Next','Callback',@next_page);
prev_page_button = uicontrol('Style','PushButton','Units', 'Normalized',...
    'Position',[0.01 0.05 0.05 0.05],'String','Previous','Callback',@prev_page);

win_pos_slider = uicontrol('Style','slider','Units', 'Normalized',...
    'Position',[0.35,0.85 - gap.*(comps_n_page),0.3,0.05],'Value',1,'Callback',@update_slider);

xlimits = [0,size([comp150.trial{1,:}],2)./comp150.fsample];

win_width = 25;
win_pos = 1;
update_slider
[px,py,pz,p,map_resource] = create_plot

uiwait

    function update_bad_ch(src,event)
        for nn = 1:size(t.Data,1)
            if t.Data{nn,3}
                comps_good(nn) = 1;
                t.Data{nn,2} = 'Good';
                table_color(nn,:) = [0.6 1 0.6];
            else
                comps_good(nn) = 0;
                t.Data{nn,2} = 'Bad';
                table_color(nn,:) = [1 0.6 0.6];
            end
        end
        set(t,'BackgroundColor',table_color)
        for ind = 1:comps_n_page
            if comps_good(curr_comps(ind)) == 0
                ax4(ind).Color = [1 0.6 0.6];
            else
                ax4(ind).Color = 'w';
            end
        end
    end

    function done_bad_comps(src,event)
        for nn = 1:size(t.Data,1)
            if strcmpi(t.Data{nn,2},'Bad')
                comps_good(nn) = 0;
            else
                comps_good(nn) = 1;
            end
        end
        idx2 = find(comps_good == 0);
        bad_comps = idx2';
        close(f)
    end

    function next_page(src,event)
        currpage = currpage + 1;
        curr_comps = 1 + mod(((currpage - 1)*comps_n_page) : currpage*comps_n_page-1,n_comps);
        curr_comps
        update_plot
    end

    function prev_page(src,event)
        currpage = currpage - 1;
        curr_comps = 1 + mod(((currpage - 1)*comps_n_page) : currpage*comps_n_page-1,n_comps);
        curr_comps
        update_plot
    end

    function [px,py,pz,p,map_resource] = create_plot
        for compi = 1:comps_n_page
            % x comps
            subplot(ax1(compi))
            x_topo = nan(size(lay.label));
            x_topo(ch_table.slot_no(ch_table.isx==1),1) = comp150.topo(ch_table.isx==1,...
                curr_comps(compi));
            [xn, yn, cn, xi, yi,px(compi)] = plot_fieldmaps_LR(x_topo,lay,ax1(compi));
            map_resource = {xn, yn, cn, xi, yi};

            % y comps
            subplot(ax2(compi))
            y_topo = nan(size(lay.label));
            y_topo(ch_table.slot_no(ch_table.isy==1),1) = comp150.topo(ch_table.isy==1,...
                curr_comps(compi));
            [~,~,~,~,~,py(compi)] = plot_fieldmaps_LR(y_topo,lay,ax2(compi));
            
            % z comps
            subplot(ax3(compi))
            z_topo = nan(size(lay.label));
            z_topo(ch_table.slot_no(ch_table.isz==1),1) = comp150.topo(ch_table.isz==1,...
                curr_comps(compi));
            [~,~,~,~,~,pz(compi)] = plot_fieldmaps_LR(z_topo,lay,ax3(compi));
            
            % time course
            subplot(ax4(compi))         
            time = (0:size(compdata,2)-1)./comp150.fsample;
            p(compi) = plot(ax4(compi),time,compdata(curr_comps(compi),:),'k');
            ylabel(sprintf('Comp %d',curr_comps(compi)));
            ax4(compi).Box='off';
            if comps_good(curr_comps(compi))==0;ax4(compi).Color = [1,0.6,0.6];end
            
            drawnow
        end
        update_slider
    end

    function update_plot(src,event)
        for compi = 1:comps_n_page
            x_topo = nan(size(lay.label));
            x_topo(ch_table.slot_no(ch_table.isx==1),1) = comp150.topo(ch_table.isx==1,...
                curr_comps(compi));
            cn = [x_topo(~isnan(x_topo))',zeros(1,size(lay.outline{1,1},1))];
            xn = [lay.pos(~isnan(x_topo),1)',lay.outline{1,1}(:,1)'];
            yn = [lay.pos(~isnan(x_topo),2)',lay.outline{1,1}(:,2)'];
            [~,~,ci] = griddata(xn, yn, cn, map_resource{4}, map_resource{5},'cubic');
            px(compi).CData = ci;
            
            % y comps
            y_topo = nan(size(lay.label));
            y_topo(ch_table.slot_no(ch_table.isy==1),1) = comp150.topo(ch_table.isy==1,...
                curr_comps(compi));
            cn = [y_topo(~isnan(y_topo))',zeros(1,size(lay.outline{1,1},1))];
            xn = [lay.pos(~isnan(y_topo),1)',lay.outline{1,1}(:,1)'];
            yn = [lay.pos(~isnan(y_topo),2)',lay.outline{1,1}(:,2)'];
            [~,~,ci] = griddata(xn, yn, cn, map_resource{4}, map_resource{5},'cubic');
            py(compi).CData = ci;
            % z comps
            z_topo = nan(size(lay.label));
            z_topo(ch_table.slot_no(ch_table.isz==1),1) = comp150.topo(ch_table.isz==1,...
                curr_comps(compi));
            cn = [z_topo(~isnan(z_topo))',zeros(1,size(lay.outline{1,1},1))]
            xn = [lay.pos(~isnan(z_topo),1)',lay.outline{1,1}(:,1)'];
            yn = [lay.pos(~isnan(z_topo),2)',lay.outline{1,1}(:,2)'];
            [~,~,ci] = griddata(xn, yn, cn, map_resource{4}, map_resource{5},'cubic');
            pz(compi).CData = ci;
            
            compdata = [comp150.trial{1,:}];
            p(compi).YData = compdata(curr_comps(compi),:);
            ax4(compi).YLabel.String = sprintf('Comp %d',curr_comps(compi));
            if comps_good(curr_comps(compi))==0;ax4(compi).Color = [1,0.6,0.6];end
            drawnow 
        end
        update_slider
    end

    function update_slider(src,event)
        win_pos = get(win_pos_slider,'Value');
        x_top = xlimits(2) - (win_pos - 1)*(win_width - xlimits(2));
        x_bottom = x_top - win_width;
        set(ax4(1:comps_n_page),'XLim',[x_bottom,x_top])
    end
end

