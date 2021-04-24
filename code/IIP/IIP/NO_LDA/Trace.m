load SBV_NO_FS1

ax(1) = plot(a,S_V,'--o','Color',[220/255 0 0],'LineWidth',1.5,'MarkerSize',10);
%hold on
%ax(2)= plot(a(4),S_V(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'Stroke volume (ml)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Times New Roman';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Times New Roman';

ax(1).FontName = 'Times New Roman';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[500 1500];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[500,750,1000,1250,1500];


box off

hold on

load SBV_FS

ax(2) = plot(a,S_V,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
%hold on
%ax(2)= plot(a(4),S_V(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
hold on


ax(2) = gca; 
ax(2).XColor = 'k';
ax(2).YColor = 'k';
ax(2).YLabel.String = 'Stroke volume (ml)';
ax(2).YLabel.FontSize = 18;
ax(2).YLabel.FontName='Times New Roman';
ax(2).XLabel.String = 'Stressed blood volume (ml)';
ax(2).XLabel.FontSize = 18;
ax(2).XLabel.FontName='Times New Roman';

ax(2).FontName = 'Times New Roman';
ax(2).FontSize = 18;

%ax(1).YLim=[50 130];
ax(2).XLim=[500 1500];
%ax(1).YTick = [20:20:140];
ax(2).XTick=[500,750,1000,1250,1500];
