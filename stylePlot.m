function stylePlot(ax)

% Manage input arguments
if nargin == 0
    ax = gca;
end

%% Paramters definition
fontName = "Times New Roman";

%% Style figure
f = gcf;

%  Background color
f.Color = 'black';
f.Position = [0, 0, 1000, 700];

%% Style axes
ax = gca;

% Background color
ax.Color = 'black';

% Title
ax.Title.FontName = fontName;
ax.Title.FontWeight = 'bold';
ax.Title.Color = 'white';
ax.Title.FontSize = 20;

% Labels
ax.XColor = 'white';
ax.YColor = 'white';
ax.ZColor = 'white';

ax.XLabel.FontName = fontName;
ax.XLabel.FontSize = 18;

ax.YLabel.FontName = fontName;
ax.YLabel.FontSize = 18;
ax.YLabel.Rotation = 0;

ax.ZLabel.FontName = fontName;
ax.ZLabel.FontSize = 18;

% Line format
ax.LineWidth = 1.05;

% Grid
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';

ax.GridColor = 'white';
ax.GridLineWidth = 1.5;

% Legend
leg = ax.Legend;
if ~isempty(leg)
    leg.FontName = fontName;
    leg.TextColor = 'white';
    leg.FontSize = 15;
    leg.EdgeColor = 'white';
    leg.Color = 'black';
end

% Text
txt = ax.findobj('type', 'text');
if ~isempty(txt)
    for i = 1:size(txt, 1)
        txt(i).FontName = fontName;
        txt(i).FontSize = 15;
        txt(i).FontWeight = 'bold';
        txt(i).HorizontalAlignment = 'center';
        txt(i).VerticalAlignment = 'middle';
        % txt(i).BackgroundColor = 'white';
        pos = txt(i).Position;
        txt(i).Position = [pos(1).*1.05, pos(2).*1.05, pos(3)];
    end
end

end