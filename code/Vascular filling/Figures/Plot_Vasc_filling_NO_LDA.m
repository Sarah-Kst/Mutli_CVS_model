%% Figure FS/NO FS

figure

load SBV_2018_bas.mat

t = 0:0.8:800;

ax(1) = plot(a,S_V,'--o','Color',[53/255 131/255 141/255],'LineWidth',1.5,'MarkerSize',10);
hold on
ax(2)= plot(a(7),S_V(7),'ko','MarkerFaceColor',[220/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[220/255 0 0]);
hold on

load SBV_both_2018_no_fs.mat

ax(3) = plot(a,S_V,'--o','Color',[220/255 0 0],'LineWidth',1.5,'MarkerSize',10);
hold on
ax(4)= plot(a(7),S_V(7),'ko','MarkerFaceColor',[220/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[220/255 0 0]);
hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Stroke volume';'(ml)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Times New Roman';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Times New Roman';

ax(1).FontName = 'Times New Roman';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off
