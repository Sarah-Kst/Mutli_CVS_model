%%%%%%% SBV/SV %%%%%%%


load SBV_both_2018_no_fs

t = 0:0.8:800;

ax(1) = plot(a,S_V,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
ax(2)= plot(a(7),S_V(7),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Stroke volume';'(ml)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off

%%% Preload

load SBV_both_2018_no_fs

figure

ax(1) = plot(a,LM,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
ax(2)= plot(a(7),LM(7),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
hold on

ax(1) = gca; % current axes
ax(1).XColor = 'k';
ax(1).YColor = 'k';
%h1.Color = 'k';
ax(1).YLabel.String = {'Preload';'(µm)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

% ax(1).YLim=[80 150];
ax(1).XLim=[250 1750];
% ax(1).YTick = [80:10:120];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off

%%%%%% Preload/SV %%%%%%

load SBV_both_2018_no_fs

figure

t = 0:0.8:800;

ax(1) = plot(LM,S_V,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
ax(2)= plot(LM(7),S_V(7),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Stroke volume';'(ml)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Preload (µm)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
%%ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
%ax(1).XTick=[250,500,750,1000,1250,1500,1750];

ax(1).XLim=[0.98 1.1];


box off

hold off

%%%%% Pmax/SBV

% load Test.mat
% 
% figure
% 
% t = 0:0.8:800;
% 
% ax(1) = plot(a,LM,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
% hold on
% % ax(2)= plot(LM(4),S_V(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% % hold on
% 
% 
% ax(1) = gca; 
% ax(1).XColor = 'k';
% ax(1).YColor = 'k';
% ax(1).YLabel.String = 'Preload (µm)';
% ax(1).YLabel.FontSize = 18;
% ax(1).YLabel.FontName='Candara';
% ax(1).XLabel.String = 'Stress Blood Volume (ml)';
% ax(1).XLabel.FontSize = 18;
% ax(1).XLabel.FontName='Candara';
% 
% ax(1).FontName = 'Candara';
% ax(1).FontSize = 18;
% 
% %ax(1).YLim=[7.5 10];
% %ax(1).XLim=[250 1750];
% %ax(1).YTick = [20:20:140];
% ax(1).XTick=[250,500,750,1000,1250,1500,1750];
% 
% 
% box off
% 
% hold off


%%%%%% Preload/SV %%%%%%

% load SBV_both_2018_no_fs
% 
% figure
% 
% t = 0:0.8:800;
% 
% ax(1) = plot(LM,S_V,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
% hold on
% ax(2)= plot(LM(4),S_V(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on
% 
% 
% ax(1) = gca; 
% ax(1).XColor = 'k';
% ax(1).YColor = 'k';
% ax(1).YLabel.String = 'Stroke volume (ml)';
% ax(1).YLabel.FontSize = 18;
% ax(1).YLabel.FontName='Candara';
% ax(1).XLabel.String = 'Preload (µm)';
% ax(1).XLabel.FontSize = 18;
% ax(1).XLabel.FontName='Candara';
% 
% ax(1).FontName = 'Candara';
% ax(1).FontSize = 18;
% 
% %ax(1).YLim=[50 130];
% %%ax(1).XLim=[250 1750];
% %ax(1).YTick = [20:20:140];
% %ax(1).XTick=[250,500,750,1000,1250,1500,1750];
% 
% 
% box off
% 
% hold off

%%%%% Preload/SV

% figure
% 
% load SBV_both_2018_no_fs
% 
% t = 0:0.8:800;
% 
% ax(3) = plot(EDV,S_V,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
% hold on
% ax(4)= plot(EDV(4),S_V(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on
% 
% 
% ax(3) = gca; 
% ax(3).XColor = 'k';
% ax(3).YColor = 'k';
% ax(3).YLabel.String = 'Stroke volume (ml)';
% ax(3).YLabel.FontSize = 18;
% ax(3).YLabel.FontName='Candara';
% ax(3).XLabel.String = 'End-diastolic volume (ml)';
% ax(3).XLabel.FontSize = 18;
% ax(3).XLabel.FontName='Candara';
% 
% ax(3).FontName = 'Candara';
% ax(3).FontSize = 18;
% 
% %ax(1).YLim=[50 130];
% %ax(3).XLim=[500 1500];
% %ax(1).YTick = [20:20:140];
% %ax(3).XTick=[500,750,1000,1250,1500];
% 
% 
% box off
% 
% hold off

%%%% Preload

% load SBV_both_2018_no_fs
% 
% figure
% 
% ax(1) = plot(a,LM,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
% hold on
% ax(2)= plot(a(4),LM(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on
% 
% ax(1) = gca; % current axes
% ax(1).XColor = 'k';
% ax(1).YColor = 'k';
% %h1.Color = 'k';
% ax(1).YLabel.String = 'Preload (µm)';
% ax(1).YLabel.FontSize = 18;
% ax(1).YLabel.FontName='Candara';
% ax(1).XLabel.String = 'Stressed blood volume (ml)';
% ax(1).XLabel.FontSize = 18;
% ax(1).XLabel.FontName='Candara';
% 
% ax(1).FontName = 'Candara';
% ax(1).FontSize = 18;
% 
% % ax(1).YLim=[80 150];
% %ax(1).XLim=[250 1750];
% % ax(1).YTick = [80:10:120];
% ax(1).XTick=[250,500,750,1000,1250,1500,1750];
% 
% box off
% 
% 
%%% Pmax & Afterload 

load SBV_both_2018_no_fs

figure

ax(1) = plot(a,Pmax,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
%ax(2)= plot(a(4),Pmax(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
%hold on
ax(1) = gca; % current axes
ax(1).XColor = 'k';
ax(1).YColor = 'k';
%h1.Color = 'k';
ax(1).YLabel.String = {'Maximal left ventricular';'pressure (mmHg)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';
%handle=title('Maximal left ventricular pressure','FontName','Candara','FontSize',16,'FontWeight','normal');
%set(handle,'Position',[900,141,0]);


ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

% ax(1).YLim=[80 150];
ax(1).XLim=[250 1750];
% ax(1).YTick = [80:10:120];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off

figure

ax(1) = plot(a,Afterload,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
%hold on
%ax(2)= plot(a(4),Afterload(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
hold on
ax(1) = gca; % current axes
ax(1).XColor = 'k';
ax(1).YColor = 'k';
%h1.Color = 'k';
ax(1).YLabel.String = {'Aortic pressure at aortic';'valve opening (mmHg)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

% handle=title({'Aortic pressure at aortic',' valve opening (afterload)'},'FontName','Candara','FontSize',16,'FontWeight','normal');
% set(handle,'Position',[830,89,0]);

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

% ax(1).YLim=[80 150];
ax(1).XLim=[250 1750];
% ax(1).YTick = [80:10:120];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off
% 
% 
% %%% Pmax - Afterload
% 
% load SBV_both_2018_no_fs
% 
% figure
% 
% ax(1) = plot(a,Pmax-Afterload,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
% %hold on
% %ax(2)= plot(a(4),Pmax(4)-Afterload(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on
% ax(1) = gca; % current axes
% ax(1).XColor = 'k';
% ax(1).YColor = 'k';
% %h1.Color = 'k';
% ax(1).YLabel.String = 'Plv,max - Pao,avo (mmHg)';
% ax(1).YLabel.FontSize = 18;
% ax(1).YLabel.FontName='Candara';
% ax(1).XLabel.String = 'Stressed blood volume (ml)';
% ax(1).XLabel.FontSize = 18;
% ax(1).XLabel.FontName='Candara';
% 
% ax(1).FontName = 'Candara';
% ax(1).FontSize = 18;
% 
% % ax(1).YLim=[32.5 34];
% % %ax(1).XLim=[250 1750];
% ax(1).YTick = [48:1:51];
% % ax(1).XTick=[250,500,750,1000,1250,1500,1750];
% 
% box off


%%% Sepsis

% clear
% 
% figure
% 
% load Total_Sepsis
% 
% ax(1) = plot(a,S_V,'--ko','LineWidth',1.5,'MarkerSize',10);
% hold on
% ax(2)= plot(a(4),S_V(4),'ko','MarkerFaceColor',[0,0,0],'MarkerSize',10);
% 
% hold on
% 
% load Total_BL
% 
% ax(1) = plot(a,S_V,'ko','LineWidth',1.5,'MarkerSize',10);
% hold on
% ax(2)= plot(a(4),S_V(4),'ko','MarkerFaceColor',[0,0,0],'MarkerSize',10);
% 
% ax(1) = gca; % current axes
% ax(1).XColor = 'k';
% ax(1).YColor = 'k';
% %h1.Color = 'k';
% ax(1).YLabel.String = 'Stroke volume (ml)';
% ax(1).YLabel.FontSize = 18;
% ax(1).YLabel.FontName='Candara';
% ax(1).XLabel.String = 'Stressed blood volume (ml)';
% ax(1).XLabel.FontSize = 18;
% ax(1).XLabel.FontName='Candara';
% 
% ax(1).FontName = 'Candara';
% ax(1).FontSize = 18;
% 
% ax(1).YLim=[40 65];
% ax(1).XLim=[350 1250];
% %ax(1).YTick = [20:20:140];
% ax(1).XTick=[250,500,750,1000,1250,1500,1750];
% 
% box off

hold off

%%%%%%% t_ejc %%%%%%%


load SBV_both_2018_no_fs
load T_both_NO_FS

t = 0:0.8:800;

figure

ax(1)= plot(a,T.t_ejc,'o','MarkerFaceColor',[150/255 185/255 220/255],'MarkerSize',8,'MarkerEdgeColor',[150/255 185/255 220/255]);
ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Blood ejection';'duration (ms)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
ax(1).YTick = [236:4:248];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off

%%%%%%% t_dur %%%%%%%


load SBV_both_2018_no_fs
load T_both_NO_FS

t = 0:0.8:800;

figure


ax(1)= plot(a,T.t_dur,'o','MarkerFaceColor',[150/255 185/255 220/255],'MarkerSize',8,'MarkerEdgeColor',[150/255 185/255 220/255]);
ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Isovolumic contraction';'duration (ms)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [236:4:248];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off

%%%%%%% Mitral closing %%%%%%%


load SBV_both_2018_no_fs
load T_both_NO_FS

t = 0:0.8:800;

figure


ax(1)= plot(a,T.t_init,'o','MarkerFaceColor',[150/255 185/255 220/255],'MarkerSize',8,'MarkerEdgeColor',[150/255 185/255 220/255]);
ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Mitral closing time';'(ms)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [236:4:248];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off

%%%%%%% Aortic opening %%%%%%%


load SBV_both_2018_no_fs
load T_both_NO_FS

t = 0:0.8:800;

figure


ax(1)= plot(a,T.t_ao,'o','MarkerFaceColor',[150/255 185/255 220/255],'MarkerSize',8,'MarkerEdgeColor',[150/255 185/255 220/255]);
ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Aortic closing time';'(ms)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [236:4:248];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];


box off

hold off

%%% Half_sarc

load SBV_both_2018_no_fs

figure

t = 0:0.8:800;

ax(1) = plot(a,LM,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
% ax(2)= plot(a(4),LM(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Preload';'(µm)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off

hold off

%%% Fmax

load SBV_both_2018_no_fs

figure

t = 0:0.8:800;

ax(1) = plot(a,Fmax,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
%ax(2)= plot(a(4),Fmax(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
%hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Maximal normalized active force';'(mN/mm²)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

ax(1).YLim=[7.5 10];
ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off

hold off

%%% Pmax

load SBV_both_2018_no_fs

figure

t = 0:0.8:800;

ax(1) = plot(a,Pmax,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
% ax(2)= plot(a(4),Pmax(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Maximal left ventricular pressure';'(mmHg)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off

hold off

%%% Delta P

load SBV_both_2018_no_fs

figure

t = 0:0.8:800;

ax(1) = plot(a,Pmax-Afterload,'--o','Color',[150/255 185/255 220/255],'LineWidth',1.5,'MarkerSize',10);
hold on
% ax(2)= plot(a(4),Pmax(4),'ko','MarkerFaceColor',[192/255 0 0],'MarkerSize',10,'MarkerEdgeColor',[192/255 0 0]);
% hold on


ax(1) = gca; 
ax(1).XColor = 'k';
ax(1).YColor = 'k';
ax(1).YLabel.String = {'Pressure difference';'(mmHg)'};
ax(1).YLabel.FontSize = 18;
ax(1).YLabel.FontName='Candara';
ax(1).XLabel.String = 'Stressed blood volume (ml)';
ax(1).XLabel.FontSize = 18;
ax(1).XLabel.FontName='Candara';

ax(1).FontName = 'Candara';
ax(1).FontSize = 18;

%ax(1).YLim=[50 130];
ax(1).XLim=[250 1750];
%ax(1).YTick = [20:20:140];
ax(1).XTick=[250,500,750,1000,1250,1500,1750];

box off

hold off

%%% Export_fig

% FigHandle = figure(1);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(1);
% export_fig -painters -r400 -transparent VF_SV.png;

% FigHandle = figure(2);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(2);
% export_fig -painters -r400 -transparent VF_preload.png;

% FigHandle = figure(3);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(3);
% export_fig -painters -r400 -transparent VF_preload_SV.png;
% 
% FigHandle = figure(4);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(4);
% export_fig -painters -r400 -transparent VF_Pmax.png;
% 
% FigHandle = figure(5);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(5);
% export_fig -painters -r400 -transparent VF_Pao.png;
% 
% FigHandle = figure(6);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(6);
% export_fig -painters -r400 -transparent VF_tejc.png;
% 
% FigHandle = figure(7);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(7);
% export_fig -painters -r400 -transparent VF_tdur.png;
% 
% FigHandle = figure(11);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(11);
% export_fig -painters -r400 -transparent VF_Fmax.png;

% FigHandle = figure(13);set(FigHandle, 'Position', [100, 100, 600, 400]);hFigure = figure(13);
% export_fig -painters -r400 -transparent VF_diff.png;









