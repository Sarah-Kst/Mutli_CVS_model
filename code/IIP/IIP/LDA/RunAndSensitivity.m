function [] = RunAndSensitivity()


% global tVaoOpen;
% global tVpvOpen;

imax=1;  % to run the Multiscale Program once, with the parameters defined in NLD7ch_NewHeart_OPTIM.m
% imax=15; % to do a sensitivity analysis with respect to the 14 parameters
% imax=99;   % to study SV(SBV)

LoadVar=0; %(1/0) : (yes/no)

if imax==1 && LoadVar==0 ;typerun=1;end;  % unique calculation without Load variation 
if imax==1 && LoadVar==1 ;typerun=4;end;  % unique calculation with    Load variation 
if imax==15;typerun=2;end;  % sensitivity analysis
if imax==99;typerun=5;end;  % to study SV(SBV)
%%%typerun=3 : used in identification 

k = imax;
save k k 


% if typerun==2
res=[];
save res res
DataLoops=[];
save DataLoops DataLoops
DataFinalLoops=[];
save DataFinalLoops DataFinalLoops
% INLV=[];
% save INLV INLV;
% end

% fOptimRmt=FactOptim(1);
% fOptimRvaor=FactOptim(2);
% fOptimRsys=FactOptim(3);
% fOptimRtc=FactOptim(4);
% fOptimRpv=FactOptim(5);
% fOptimRpul=FactOptim(6);
% fOptimCaor=FactOptim(7);
% fOptimCvc=FactOptim(8);
% fOptimCpa=FactOptim(9);
% fOptimCpu=FactOptim(10);
% fOptimSBV=FactOptim(11);
% fOptimVw_lv=FactOptim(12);
% fOptimVw_rv=FactOptim(13);
% fOptimf_PassiveForce=FactOptim(14);

tic;
if typerun~=5 
    for i=1:imax
       if imax~=1;disp(sprintf('  i =   %d', i));end;
       FactOptim=ones(1,14); 
       if i~=1;FactOptim(i-1)=1.002*FactOptim(i-1);end     
       CalculatedVals=NLD7ch_NewHeart_OPTIM(FactOptim,typerun);
    end;
else
    for i=4:1:12
        k = i;
        save k k 
       if imax~=1;disp(sprintf('  i =   %d', i));end;
       FactOptim=ones(1,14); 
       if i~=1;FactOptim(11)=0.5+(i-1)*0.1;end     
       CalculatedVals=NLD7ch_NewHeart_OPTIM(FactOptim,typerun);
    end;
end

c=toc;
disp(sprintf('  duration for running "RunAndSensitivity"  =   %d', c))
end


