% figure

function slid_test_delme
x = (0:20*100)./100;
plot(x,sin(x))
ylimits = xlim
win_width = (ylimits(2) - ylimits(1))/5;
win_pos = 1;
y_top = ylimits(2) - (win_pos - 1)*(win_width - ylimits(2))
y_bottom = y_top - win_width
xlim([y_bottom,y_top])
ax = gca;
win_pos_slider = uicontrol('Style','slider','Units', 'Normalized',...
    'Position',[0.6 0.05 0.01 0.6],'Value',1,'Callback',@update_slider);
uiwait



function update_slider(src,event)
    win_pos = get(win_pos_slider,'Value');
    y_top = ylimits(2) - (win_pos - 1)*(win_width - ylimits(2));
    y_bottom = y_top - win_width;
    xlim(ax,[y_bottom,y_top])
end
end