% close all; 38652
% clear all;
% clc; 
IdRes=[];
save IdRes IdRes
c=clock;
% format long
% copyfile('FinalConditions_NLD7ch_NewHeart_BCKP.mat','FinalConditions_NLD7ch_NewHeart.mat')
NbFactToIdentify=9;
FactToIdentify=ones(1,NbFactToIdentify);
% A=-eye(NbFactToIdentify);b=-0.25*ones(NbFactToIdentify,1);
% A=[ ];b=[ ];
% Aeq=[ ];beq=[ ];
% lb=0.1*ones(NbFactToIdentify,1);
% ub=Inf*ones(NbFactToIdentify,1);
% nonlcon=[ ];
% options_fmincon = optimset('Algorithm','interior-point','TolFun',4e-3,'TolX',1e-2,'Display','iter','UseParallel','always');
% options_fmincon = optimset('Algorithm','sqp','TolFun',1e-3,'TolX',1e-3,'Display','iter','UseParallel','always');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options_fminsearch = optimset('TolFun',1e-4,'TolX',1e-4,'Display','iter');
% options_fminsearch = optimset('TolFun',0.5e-3,'TolX',1e-3,'Display','iter');
% options_fminsearch = optimset('TolFun',8e-7,'TolX',8e-7,'Display','iter');
% options_fminsearch = optimset('TolFun',5e-6,'TolX',5e-6,'Display','iter');
% options_fminsearch = optimset('TolFun',5e-6,'TolX',5e-6,'Display','iter');
%options_fminsearch = optimset('TolFun',1e-6,'TolX',1e-6,'Display','iter');
% options_lsqnonlin = optimset('algorithm','trust-region-reflective');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('------------------------------------------------------'));
disp(sprintf('-----------------Start of minimization----------------'));
disp(sprintf('- %d -',c));
disp(sprintf('------------------------------------------------------'));

[FactIdentified,fval] = myfminsearch(@objectif,FactToIdentify,options_fminsearch)
% [FactIdentified,fval] = lsqnonlin(@objectif,FactToIdentify,...
%    0.1*FactToIdentify,10*FactToIdentify,options_lsqnonlin)

disp(sprintf('------------------------------------------------------'));
disp(sprintf('-----------------  End of minimization----------------'));
disp(sprintf('------------------------------------------------------'));
for i=1:NbFactToIdentify;
disp(sprintf('                   FINAL Values of the unknowns n°: %d = %0.16e',i,FactIdentified(i)));
end
