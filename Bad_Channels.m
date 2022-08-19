function [ch_table_new] = Bad_Channels(data_f,ch_table,fs,hp,lp)
%Plots channels and provides ui to remove bad channels 
% data_f... nsamps x nchans matrix

% time...time vector

% ch_table...table containing columns name, status, with channel names and
% good/bad status

% hp/lp...filter cut off values for band pass filter 
if size(data_f,2) > size(data_f,1)
    error("Make sure the data matrix is in (N_sample x N_channel) form!")
end
Nchans = size(data_f,2);
fall_ch = figure;
fall_ch.Units = 'normalized';fall_ch.Position=[0,0,1,1];
fall_ch.Name = 'All channels';
sep = max(max(data_f))./1;
chan_space = sep*1.1;

time = linspace(0,size(data_f,1)./fs,size(data_f,1));
p1 = plot(time,data_f+(chan_space*[0:1:Nchans-1]));
ax = gca;
xlim([time(1), time(end)])

fall_ch.CurrentAxes.Position = [0.08 0.05 0.45 0.95];
fall_ch.CurrentAxes.FontSize = 10;
%             fall_ch.CurrentAxes.FontSizeMode = 'auto';
fall_ch.CurrentAxes.Units = 'normalized';
fall_ch.CurrentAxes.YTick = 0:chan_space:chan_space*Nchans-1;
for nc = 1:Nchans
    fall_ch.CurrentAxes.YTickLabel{nc} = [ch_table.name{nc},sprintf(' - Ch %u',nc)];
end
fig = gcf;
fig.Color = [1,1,1];
xlabel('Time (s)');ylabel('Channels')
title(sprintf('All channels - Filtered data [%u - %u] Hz',hp,lp))

% Load previous bad channels
bad_ch = find(startsWith(ch_table.status,'bad'));
for n = 1:Nchans
    table_data{n,1} = ch_table.name{n};
    table_data{n,2} = ch_table.status{n};
    table_data{n,3} = true;
    table_color(n,:) = [0.6 1 0.6];
end

for b_i = 1:size(bad_ch,1)
table_data{bad_ch(b_i),3} = false;
table_color(bad_ch(b_i),:) = [1 0.6 0.6];
p1(bad_ch(b_i)).Color = 'r';
% p1(bad_ch(b_i)).LineStyle = '--';
p1(bad_ch(b_i)).LineWidth=0.1;
end
t = uitable('Data',table_data,'ColumnName',{'Channel Name','Status','Toggle'},...
    'Units', 'Normalized', 'Position',[0.7, 0.1, 0.2, 0.8],'BackgroundColor',table_color);
set(t,'ColumnEditable',[false false true])

update_button = uicontrol('Style','PushButton','Units', 'Normalized',...
    'Position',[0.7 0.05 0.05 0.05],'String','Update','Callback',@update_bad_ch);
done_button = uicontrol('Style','PushButton','Units', 'Normalized',...
    'Position',[0.75 0.05 0.1 0.05],'String','Remove bad channels','Callback',@Done_bad_ch);

win_pos_slider = uicontrol('Style','slider','Units', 'Normalized',...
    'Position',[0.6 0.05 0.01 0.6],'Value',1,'Callback',@update_slider);

ylimits = ([-chan_space,(chan_space*Nchans)+chan_space]);
win_width = (ylimits(2) - ylimits(1))/5;
win_pos = 1;
update_slider

uiwait

ch_table_new = ch_table;
for b_i = 1:size(bad_ch,1)
    ch_table_new.status{bad_ch(b_i)} = 'bad';
end


function update_bad_ch(src,event)
    fprintf('t data size during update: %d',size(t.Data,1));
for nn = 1:size(t.Data,1)
    if t.Data{nn,3}
        idx(nn) = 1;
        t.Data{nn,2} = 'Good';
        p1(nn).Visible = 'on';
        table_color(nn,:) = [0.6 1 0.6];
    else
        idx(nn) = 0;
        t.Data{nn,2} = 'Bad';
        p1(nn).Visible = 'off';
        table_color(nn,:) = [1 0.6 0.6];
    end
end
set(t,'BackgroundColor',table_color)
end

function Done_bad_ch(src,event)
for nn = 1:size(t.Data,1)
    if strcmpi(t.Data{nn,2},'Bad')
        idx(nn) = 0;
    else
        idx(nn) = 1;
    end
end
idx2 = find(idx == 0);

bad_ch = idx2';
close(fig)
end

function update_slider(src,event)
    win_pos = get(win_pos_slider,'Value');
    y_top = ylimits(2) - (win_pos - 1)*(win_width - ylimits(2));
    y_bottom = y_top - win_width;
    ylim(ax,[y_bottom,y_top])
end
end