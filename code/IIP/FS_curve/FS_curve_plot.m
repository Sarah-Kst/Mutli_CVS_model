%% FS curve alone

figure

load T_FS_curve

ax(1) = plot(T_FS_curve.Preload,T_FS_curve.S_V,'--o','Color',[220/255 0/255 0/255],'LineWidth',1.5,'MarkerSize',10);
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Stroke volume';'(ml)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'Preload (µm)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[0.99 1.15];
%ax(1).YTick = [20:20:140];
%ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off

FigHandle = figure(1);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(1);
%export_fig -painters -r200 -transparent FS_curve.png;












