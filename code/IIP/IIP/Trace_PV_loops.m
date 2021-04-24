%% Protocol IIP
% ----------------------------------------------------------------------

figure % active force

load Variables_FS7d

time = t - min(t);

ax(1) = plot(time,Fm_lv(:,5)-Fparall_lv(:,5),'Color',[220/255 0 0],'LineWidth',1.5);
hold on

load Variables_BL_dtoutput_fzero
time = t - min(t);

ax(1) = plot(time,Fm_lv(:,10)-Fparall_lv(:,10),'Color',[53/255 131/255 141/255],'LineWidth',1.5);
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'LV active force (mN/mm²)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'Time (ms)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

box off

hold on

load Variables_NO_FS7d

time = t - min(t);

ax(3) = plot(time,Fm_lv(:,5)-Fparall_lv(:,5),'Color',[220/255 0 0],'LineWidth',1.5, 'LineStyle','--');
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'LV active force (mN/mm²)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'Time (ms)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

% ----------------------------------------------------------------------

figure % PV loop

load Variables_BL_dtoutput_fzero

%time = t - min(t);

ax(1) = plot(Vlv(:,10),Plv(:,10),'Color',[53/255 131/255 141/255],'LineWidth',1.5,'MarkerSize',10);
hold on

load Variables_FS7d

%time = t - min(t);

ax(1) = plot(Vlv(:,5),Plv(:,5),'Color',[220/255 0 0],'LineWidth',1.5,'MarkerSize',10);
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'LV Pressure (mmHg)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'LV Volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

box off

hold on

load Variables_NO_FS7d

%time = t - min(t);

ax(2) = plot(Vlv(:,5),Plv(:,5),'Color',[220/255 0 0],'LineWidth',1.5,'MarkerSize',10,'LineStyle','--');
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'LV Pressure (mmHg)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'LV Volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

%% Aortic flow 

figure 

load Variables_BL_dtoutput_fzero

time = t - min(t);

ax(1) = plot(time,Qvao(:,10),'Color',[53/255 131/255 141/255],'LineWidth',1.5,'MarkerSize',10);

box off

hold on

load Variables_NO_FS7d


time = t - min(t);

ax(3) = plot(time,Qvao(:,5),'Color',[220/255 0 0],'LineWidth',1.5,'MarkerSize',10,'LineStyle','--');
hold on

load Variables_FS7d

time = t - min(t);

ax(2) = plot(time,Qvao(:,5),'Color',[220/255 0 0],'LineWidth',1.5,'MarkerSize',10);
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'Aortic flow (ml/ms)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'Time (ms)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;


%% Export fig

FigHandle = figure(1);set(FigHandle, 'Position', [100, 100, 600, 400]);%hFigure = figure(7);
% export_fig -painters -r400 -transparent aortic_FS.png;
%  
FigHandle = figure(2);set(FigHandle, 'Position', [100, 100, 600, 400]);%hFigure = figure(8);
% export_fig -painters -r400 -transparent aortic_FS2.png;

%% Aortic and ventricular pressure

figure % active force

load Variables_FS7d

time = t - min(t);

ax(1) = plot(time,Plv(:,5),'Color',[220/255 0 0],'LineWidth',1.5);
hold on
ax(4) = plot(time,Pao(:,5),'Color',[220/255 0 0],'LineWidth',1.5);

hold on

load Variables_BL_dtoutput_fzero
time = t - min(t);

ax(1) = plot(time,Plv(:,10),'Color',[53/255 131/255 141/255],'LineWidth',1.5);
hold on
ax(7) = plot(time,Pao(:,10),'Color',[53/255 131/255 141/255],'LineWidth',1.5);
hold on

ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'LV active force (mN/mm²)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'Time (ms)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

box off

hold on

load Variables_NO_FS7d

time = t - min(t);

ax(3) = plot(time,Plv(:,5),'Color',[220/255 0 0],'LineWidth',1.5, 'LineStyle','--');
hold on
ax(5) = plot(time,Pao(:,5),'Color',[220/255 0 0],'LineWidth',1.5, 'LineStyle','--');

hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = 'Pressure (mmHg)';
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Calibri';
ax(1).XLabel.String = 'Time (ms)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Calibri';

ax(1).FontName = 'Calibri';
ax(1).FontSize = 18;

