function output_NLD7ch_NewHeart_OPTIM = NLD7ch_NewHeart_OPTIM(FactOptim,typerun)

close all;

nRes=1;IdRes(nRes,1:40)=1;
iBeatLoadVar=1e69;


if typerun==1 || typerun==2 || typerun==4 || typerun==5; % always ecept for identification
%   load IdRes; % read identified factors in file IdRes and use line nRes
%   icolumn=size(IdRes,2)-7;
%   nRes=find(IdRes(:,icolumn)== min(IdRes(:,icolumn)))    % if nRes=-1, all factors are changed to 1
end

if typerun==1; % calcul unique without Load Variation
  IntermediatePlot=1;
  NbBeats = 10;
  deltaSVmax=0.0001;
%deltaSVmax=0.00025*10^9;disp(sprintf('  Attention in NLD7ch_NewHeart_OPTIM.m, line 29  '))
end;

if typerun==4; % calcul unique with Load Variation
  iBeatLoadVar=2;
  IntermediatePlot=1;
  NbBeats = 24;
  deltaSVmax=1e20;
end;

if typerun==2; % sensibility analysis
  IntermediatePlot=0;
  NbBeats = 40;
  deltaSVmax=0.00015;
end;

if typerun==3;  % identification
   IntermediatePlot=0;
  NbBeats = 15;
  deltaSVmax=0.0001;
end;

if typerun==5; % DataLoops
  IntermediatePlot=-1;
  NbBeats = 10;
  deltaSVmax=0.00015;
  
end;

j=1;
fOptimRmt=FactOptim(1)*0.387350554472977;%*IdRes(nRes,j);j=j+1;
fOptimRvao=FactOptim(2)*1.034575712449004;%*IdRes(nRes,j);j=j+1;
fOptimRsys=FactOptim(3)*1.668626835531069*1.00580716212824*1.00553719931612*IdRes(nRes,j);j=j+1;
fOptimRtc=FactOptim(4)*1.085878487374071;%*IdRes(nRes,j);j=j+1;
fOptimRpv=FactOptim(5)*0.411067820428896;%*IdRes(nRes,j);j=j+1;
fOptimRpul=FactOptim(6)*0.416440230111030*1.53385508620215*1.01890752200915*IdRes(nRes,j);j=j+1;
fOptimCao=FactOptim(7)*1.879407473472535*0.600392587239388*0.977576940890959*IdRes(nRes,j);j=j+1;
fOptimCvc=FactOptim(8)*3.409415892361760*1.01220961191878*0.974715728991473*IdRes(nRes,j);j=j+1;
fOptimCpa=FactOptim(9)*1.283176081892981*0.933929517415392*0.860463801218685*IdRes(nRes,j);j=j+1;
fOptimCpu=FactOptim(10)*2.422202057600711*0.922932805287515*1.03966274421367*IdRes(nRes,j);j=j+1;
fOptimSBV=FactOptim(11)*1.610763727230768*0.857826462983265*0.943098006553602*IdRes(nRes,j);j=j+1;
fOptimVw_lv=FactOptim(12)*0.111753863404668*1.66372932480168*0.989415047493162*IdRes(nRes,j);j=j+1;
fOptimVw_rv=FactOptim(13)*0.052161988343849*1.00335816001536*1.53728768007854*IdRes(nRes,j);j=j+1;

fOptimf_PassiveForce=FactOptim(14)*0.3;

HR = 75;   

f_PassiveForce=fOptimf_PassiveForce;
f_Rref_SARC=1;
TAU_Bioch_lv=1;
fGeom_lv=1.; fGeom_rv=1.; fGeom_la=1;
factGamav=1;


Inert_lv=0;
Inert_rv=Inert_lv/10;

fR=1.2;  
    fRpul=fR;
    fRsys=fR;
    fRvlv=fR;
fc=1;
    fcPUL=fc;
    fcSYS=fc;
fW=1;

VoT=5000; 

SBV=fOptimSBV*(1145.928421156985-423.93);%*0.5;  %SBV=485;

Rsystot=(600)*fRsys; 
Caotot=1.377*fcSYS; 
Fact_CaoCas=1;
Fact_RsysRas=1;

Rsys=fOptimRsys*Fact_RsysRas*Rsystot*0.989476210397302*1.1448789996563795e+000;

fOptimRas=1;Ras=fOptimRas*(1-Fact_RsysRas)*Rsystot;
Cao=fOptimCao*Fact_CaoCas*Caotot*6.2871437024222943e-001;
fOptimCas=1;Cas=fOptimCas*(1-Fact_CaoCas)*Caotot;

fOptimRprox=1;Rprox=fOptimRprox*(1e69)*fRpul;
Rpul=fOptimRpul*140.2926*fRpul;

Rmt=fOptimRmt*17.28243*fRvlv*2.75;
Rvao=fOptimRvao*38.6313*fRvlv;
Rtc=fOptimRtc*24.3987*fRvlv/2.75;
Rpv=fOptimRpv*7.11629*fRvlv;


Cvc=fOptimCvc*20.24*fcSYS*0.986158709869453*1.1780663249742223e+000; 
Cpa=fOptimCpa*2.356*fcPUL; 
Cpu=fOptimCpu*fcPUL*50.215644054204908/5.;

alfaIN=0.001; betaIN=4; 
alfaOUT=0.001; betaOUT=4; 

betaINleft=betaIN/12;
betaINright=betaIN/12;
betaOUTleft=betaOUT;
betaOUTright=betaOUT*12;

Gamav_lv=0.00005*factGamav; Vo_lv=80; 
Gamav_rv=0.00005*factGamav; Vo_rv=80; 
Gamav_la=0.00005*factGamav; Vo_la=25; 

fOptimVw_la=1;Vw_la=fOptimVw_la*1000*fW;
Lr_la=0.97;
WantedLm_la_m=0.98;
WantedLm_la_M=1.115;
WantedVla_m=40; 
RrefSARC_la=2.6*f_Rref_SARC;
Vw2refSARC_la=Vsph(RrefSARC_la)-WantedVla_m;
NSARC_la=2*pi*RrefSARC_la/(WantedLm_la_m*1e-4);

Vw_lv=fOptimVw_lv*1500*fW*1.03706734222564*1.1701588668467933e+000; 
Lr_lv=0.97;
WantedLm_lv_m=0.93;
WantedLm_lv_M=1.08;
WantedVlv_m=60; 
% RrefSARC_lv=3.485*f_Rref_SARC;
%RrefSARC_lv=3.1174*f_Rref_SARC;

RrefSARC_lv=(60*(2*pi)^3/(4/3*pi*(WantedLm_lv_M^3-WantedLm_lv_m^3)))^(1/3)*WantedLm_lv_m/2/pi;
Vw2refSARC_lv=Vsph(RrefSARC_lv)-WantedVlv_m;
NSARC_lv=2*pi*RrefSARC_lv/(WantedLm_lv_m*1e-4);

Vw_rv=fOptimVw_rv*600*fW*1.02853834556684*9.7291704365546816e-001; 
Lr_rv=0.97;
WantedVrv_m=60; 
WantedLm_rv_m=0.93;
WantedLm_rv_M=1.08;
% RrefSARC_rv=3.485*f_Rref_SARC;
%RrefSARC_rv=3.1174*f_Rref_SARC;
RrefSARC_rv=(60*(2*pi)^3/(4/3*pi*(WantedLm_rv_M^3-WantedLm_rv_m^3)))^(1/3)*WantedLm_rv_m/2/pi;
Vw2refSARC_rv=Vsph(RrefSARC_rv)-WantedVrv_m;
NSARC_rv=2*pi*RrefSARC_rv/(WantedLm_rv_m*1e-4);
 
paramsHemo=[
f_PassiveForce...
f_Rref_SARC...
TAU_Bioch_lv...
fGeom_lv...
fGeom_rv...
Inert_lv...
Inert_rv...
SBV...  
Rsys...
Ras...
Cao...
Cas...
Rpul... 
Rprox...
Rmt...
Rvao...
Rtc...
Rpv...
Cvc... 
Cpa...
Cpu...
alfaIN...
alfaIN...
alfaOUT...
alfaOUT...
betaINleft...
betaINright...
betaOUTleft...
betaOUTright...
Gamav_lv... 
Vo_lv... 
Gamav_rv... 
Vo_rv... 
Gamav_la... 
Vo_la... 
Vw_la...
Lr_la...
Vw2refSARC_la...
NSARC_la...
Vw_lv... 
Lr_lv...
Vw2refSARC_lv...
NSARC_lv...
Vw_rv... 
Lr_rv...
Vw2refSARC_rv...
NSARC_rv...  %47
WantedLm_lv_m...  %48
WantedLm_lv_M...  %49
WantedVlv_m... 
WantedLm_rv_m...
WantedLm_rv_M...
WantedVrv_m... 
WantedLm_la_m...
WantedLm_la_M...
WantedVla_m... 
iBeatLoadVar... %57
];

load FinalConditions_NLD7ch_NewHeart_BL     
%load Copy_of_FinalConditions_NLD7ch_NewHeart

y0=yfinal;
nbODEs=length(y0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% mo;%----------------(1)INa
% ho;%------------------(2)INa
% jo;%-------------------(3)INa 
% do;%------------------(4)ICaL 
% fo;%------------------(5)ICaL 
% fgo;%----------------(6)ICaL
% fcao;%---------------(7)ICaL
% ro;%------------------(8)Ito
% so;%---------------------(9)Ito
% xso;%---------------(10)IKs
% xrao;%--------------(11)IKr
% xrbo;%-----------------(12)IKr
% Reo;%----------------(13)Irel
% Casro;%---------------(14)[Ca]r
% Cajo;%------------------(15)[Ca]j
% Caio;%----------------(16)[Ca]i
% Kio;%---------------(17)[K]i
% Naio;%---------------(18)[Na]i
% Vmo;%------------------(19)Vm

% TSao; %--------------(20)Mech (Left ventricle)
% TSpo; %-------------(21)Mech
% TSwo; %-------------(22)Mech
% TSro; %-------------(23)Mech
% Xpo;%-----------------(24)Mech
% Xwo;%-----------------(25)Mech

% Vlao;%-----------------(26)Load
% Vlvo;%-------------------(27)Load
% Paoo;%-----------------(28)Load
% Faoo;%--------------(29)Load
% Parso;%-----------------(30)Load

% Oeo;%-----------------(31)Irel
% Ieo;%-----------------(32)Ir

% Vrvo;%-------------------(33)Load
% Ppao;%-------------------(34)Load
% Ppuo;%--------------------(35)Load
% Pvs;%-------------------(36)Load

% TSao_la; %--------------(37)Mech_la (Left atrium)
% TSpo_la; %-------------(38)Mech_la
% TSwo_la; %-------------(39)Mech_la
% TSro_la; %-------------(40)Mech_la
% Xpo_la;%-----------------(41)Mech_la
% Xwo_la;%-----------------(42)Mech_la

% TSao_rv; %--------------(43)Mech (Right ventricle)
% TSpo_rv; %-------------(44)Mech
% TSwo_rv; %-------------(45)Mech
% TSro_rv; %-------------(46)Mech
% Xpo_rv;%-----------------(47)Mech
% Xwo_rv;%-----------------(48)Mech


% WantedPrvmean=17.2964;
% WantedVrvmean=95.9417;
% factVw_rv=1;
% save WantedPrvmean WantedPrvmean; 
% save WantedVrvmean WantedVrvmean; 
% save factVw_rv factVw_rv; 

deltaSV=2*deltaSVmax;
flag=1;
flagmax=40;

%while abs(deltaSV)>deltaSVmax && flag~=flagmax
if flag~=1; NbBeats=10; end

NbBeats = 5;

for iBeat=1:NbBeats
    outputCalculateNextBeat=CalculateNextBeat(iBeat,NbBeats,HR,y0,paramsHemo,IntermediatePlot); 
    y0=outputCalculateNextBeat(1:nbODEs);
    
    disp(sprintf('Beat n° %d',iBeat));
end

yfinal=y0;
save FinalConditions_NLD7ch_NewHeart yfinal

j=1;
Plv_min =outputCalculateNextBeat(j+nbODEs);j=j+1;
Plv_max =outputCalculateNextBeat(j+nbODEs);j=j+1;
Plv_mean=outputCalculateNextBeat(j+nbODEs);j=j+1;
Plv_ampl=outputCalculateNextBeat(j+nbODEs);j=j+1;
Vlv_min =outputCalculateNextBeat(j+nbODEs);j=j+1;
Vlv_mean =outputCalculateNextBeat(j+nbODEs);j=j+1;
SVlv    =outputCalculateNextBeat(j+nbODEs);j=j+1;
Pao_min =outputCalculateNextBeat(j+nbODEs);j=j+1;
Pao_max =outputCalculateNextBeat(j+nbODEs);j=j+1;
Pao_mean=outputCalculateNextBeat(j+nbODEs);j=j+1;
Pao_ampl=outputCalculateNextBeat(j+nbODEs);j=j+1;
Vao_mean=outputCalculateNextBeat(j+nbODEs);j=j+1;
Pvc_min =outputCalculateNextBeat(j+nbODEs);j=j+1; 
Pvc_max =outputCalculateNextBeat(j+nbODEs);j=j+1;
Pvc_mean=outputCalculateNextBeat(j+nbODEs);j=j+1; 
Pvc_ampl=outputCalculateNextBeat(j+nbODEs);j=j+1;
Vvc_mean=outputCalculateNextBeat(j+nbODEs);j=j+1; 
Prv_min =outputCalculateNextBeat(j+nbODEs);j=j+1;
Prv_max =outputCalculateNextBeat(j+nbODEs);j=j+1;
Prv_mean=outputCalculateNextBeat(j+nbODEs);j=j+1;
Prv_ampl=outputCalculateNextBeat(j+nbODEs);j=j+1;
Vrv_min =outputCalculateNextBeat(j+nbODEs);j=j+1;
Vrv_mean =outputCalculateNextBeat(j+nbODEs);j=j+1;
SVrv    =outputCalculateNextBeat(j+nbODEs);j=j+1;
Ppa_min =outputCalculateNextBeat(j+nbODEs);j=j+1; 
Ppa_max =outputCalculateNextBeat(j+nbODEs);j=j+1; 
Ppa_mean=outputCalculateNextBeat(j+nbODEs);j=j+1; 
Ppa_ampl=outputCalculateNextBeat(j+nbODEs);j=j+1; 
Vpa_mean=outputCalculateNextBeat(j+nbODEs);j=j+1; 
Ppu_min =outputCalculateNextBeat(j+nbODEs);j=j+1;
Ppu_max =outputCalculateNextBeat(j+nbODEs);j=j+1; 
Ppu_mean=outputCalculateNextBeat(j+nbODEs);j=j+1;
Ppu_ampl=outputCalculateNextBeat(j+nbODEs);j=j+1; 
Vpu_mean=outputCalculateNextBeat(j+nbODEs);j=j+1;
Qmt_max=outputCalculateNextBeat(j+nbODEs);j=j+1;
Qvao_max=outputCalculateNextBeat(j+nbODEs);j=j+1;
Qtc_max=outputCalculateNextBeat(j+nbODEs);j=j+1;
Qpv_max=outputCalculateNextBeat(j+nbODEs);j=j+1;
dtLvEjection=outputCalculateNextBeat(j+nbODEs);j=j+1;

deltaSV=(SVrv-SVlv)/SVlv;
flag=flag+1;
%end

if flag==flagmax; disp('------abs(deltaSV)>deltaSVmax------ SVlv, SVrv = ');disp(SVlv);disp(SVrv); end

output_NLD7ch_NewHeart_OPTIM=[...
Plv_min Plv_max Plv_mean Plv_ampl...
Vlv_min SVlv ...
Pao_min Pao_max Pao_mean Pao_ampl...
Pvc_min Pvc_max Pvc_mean Pvc_ampl...
Prv_min Prv_max Prv_mean Prv_ampl...
Vrv_min SVrv...
Ppa_min Ppa_max Ppa_mean Ppa_ampl... 
Ppu_min Ppu_max Ppu_mean Ppu_ampl...
Qmt_max,Qvao_max,Qtc_max,Qpv_max,dtLvEjection...
                           ];
save NLD7ch_NewHeart  
load res;
QmtRelQvao=Qmt_max/Qvao_max;
QtcRelQvao=Qtc_max/Qvao_max;
QpvRelQvao=Qpv_max/Qvao_max;

if  size(res,1)==0
    res=[...
        (SVrv-SVlv)/SVlv...
        Plv_min Plv_max Plv_mean Plv_ampl...
        Vlv_min SVlv ...
        Pao_min Pao_max Pao_mean Pao_ampl...
        Pvc_min Pvc_max Pvc_mean Pvc_ampl...
        Prv_min Prv_max Prv_mean Prv_ampl...
        Vrv_min SVrv...
        Ppa_min Ppa_max Ppa_mean Ppa_ampl...
        Ppu_min Ppu_max Ppu_mean Ppu_ampl...
        Qmt_max Qvao_max Qtc_max Qpv_max...
        QmtRelQvao QtcRelQvao QpvRelQvao dtLvEjection...
        Rmt Rvao Rsys Rtc Rpv Rpul...
        Cao Cvc Cpa Cpu ...
        SBV  Vw_lv Vw_rv f_PassiveForce];
elseif  typerun==2; 
    j=2;
    Plv_minR=res(1,j);j=j+1;
    Plv_maxR=res(1,j);j=j+1;
    Plv_meanR=res(1,j);j=j+1;
    Plv_amplR=res(1,j);j=j+1;
    Vlv_minR=res(1,j);j=j+1;
    SVlvR=res(1,j);j=j+1;
    Pao_minR=res(1,j);j=j+1;
    Pao_maxR=res(1,j);j=j+1;
    Pao_meanR=res(1,j);j=j+1;
    Pao_amplR=res(1,j);j=j+1;
    Pvc_minR=res(1,j);j=j+1;
    Pvc_maxR=res(1,j);j=j+1;
    Pvc_meanR=res(1,j);j=j+1;
    Pvc_amplR=res(1,j);j=j+1;
    Prv_minR=res(1,j);j=j+1;
    Prv_maxR=res(1,j);j=j+1;
    Prv_meanR=res(1,j);j=j+1;
    Prv_amplR=res(1,j);j=j+1;
    Vrv_minR=res(1,j);j=j+1;
    SVrvR=res(1,j);j=j+1;
    Ppa_minR=res(1,j);j=j+1;
    Ppa_maxR=res(1,j);j=j+1;
    Ppa_meanR=res(1,j);j=j+1;
    Ppa_amplR=res(1,j);j=j+1;
    Ppu_minR=res(1,j);j=j+1;
    Ppu_maxR=res(1,j);j=j+1;
    Ppu_meanR=res(1,j);j=j+1;
    Ppu_amplR=res(1,j);j=j+1;
    Qmt_maxR=res(1,j);j=j+1;
    Qvao_maxR=res(1,j);j=j+1;
    Qtc_maxR=res(1,j);j=j+1;
    Qpv_maxR=res(1,j);j=j+1;
    QmtRelQvaoR=res(1,j);j=j+1;
    QtcRelQvaoR=res(1,j);j=j+1;
    QpvRelQvaoR=res(1,j);j=j+1;
    dtLvEjectionR=res(1,j);j=j+1;

    epsM1=(max(FactOptim)-1)^(-1);
    
    DelSVrel=(SVrv-SVlv)/SVlv;
    Plv_min=(Plv_min/Plv_minR-1)*epsM1;
    Plv_max=(Plv_max/Plv_maxR-1)*epsM1;
    Plv_mean=(Plv_mean/Plv_meanR-1)*epsM1;
    Plv_ampl=(Plv_ampl/Plv_amplR-1)*epsM1;
    Vlv_min=(Vlv_min/Vlv_minR-1)*epsM1;
    SVlv=(SVlv/SVlvR-1)*epsM1;
    Pao_min=(Pao_min/Pao_minR-1)*epsM1;
    Pao_max=(Pao_max/Pao_maxR-1)*epsM1;
    Pao_mean=(Pao_mean/Pao_meanR-1)*epsM1;
    Pao_ampl=(Pao_ampl/Pao_amplR-1)*epsM1;
    Pvc_min=(Pvc_min/Pvc_minR-1)*epsM1;
    Pvc_max=(Pvc_max/Pvc_maxR-1)*epsM1;
    Pvc_mean=(Pvc_mean/Pvc_meanR-1)*epsM1;
    Pvc_ampl=(Pvc_ampl/Pvc_amplR-1)*epsM1;
    Prv_min=(Prv_min/Prv_minR-1)*epsM1;
    Prv_max=(Prv_max/Prv_maxR-1)*epsM1;
    Prv_mean=(Prv_mean/Prv_meanR-1)*epsM1;
    Prv_ampl=(Prv_ampl/Prv_amplR-1)*epsM1;
    Vrv_min=(Vrv_min/Vrv_minR-1)*epsM1;
    SVrv=(SVrv/SVrvR-1)*epsM1;
    Ppa_min=(Ppa_min/Ppa_minR-1)*epsM1;
    Ppa_max=(Ppa_max/Ppa_maxR-1)*epsM1;
    Ppa_mean=(Ppa_mean/Ppa_meanR-1)*epsM1;
    Ppa_ampl=(Ppa_ampl/Ppa_amplR-1)*epsM1;
    Ppu_min=(Ppu_min/Ppu_minR-1)*epsM1;
    Ppu_max=(Ppu_max/Ppu_maxR-1)*epsM1;
    Ppu_mean=(Ppu_mean/Ppu_meanR-1)*epsM1;
    Ppu_ampl=(Ppu_ampl/Ppu_amplR-1)*epsM1;
    Qmt_max=(Qmt_max/Qmt_maxR-1)*epsM1;
    Qvao_max=(Qvao_max/Qvao_maxR-1)*epsM1;
    Qtc_max=(Qtc_max/Qtc_maxR-1)*epsM1;
    Qpv_max=(Qpv_max/Qpv_maxR-1)*epsM1;
    QmtRelQvao=(QmtRelQvao/QmtRelQvaoR-1)*epsM1;
    QtcRelQvao=(QtcRelQvao/QtcRelQvaoR-1)*epsM1;
    QpvRelQvao=(QpvRelQvao/QpvRelQvaoR-1)*epsM1;
    dtLvEjection=(dtLvEjection/dtLvEjectionR-1)*epsM1;
   
    res=[res;...
        DelSVrel...
        Plv_min Plv_max Plv_mean Plv_ampl...
        Vlv_min SVlv ...
        Pao_min Pao_max Pao_mean Pao_ampl...
        Pvc_min Pvc_max Pvc_mean Pvc_ampl...
        Prv_min Prv_max Prv_mean Prv_ampl...
        Vrv_min SVrv...
        Ppa_min Ppa_max Ppa_mean Ppa_ampl...
        Ppu_min Ppu_max Ppu_mean Ppu_ampl...
        Qmt_max Qvao_max Qtc_max Qpv_max...
        QmtRelQvao QtcRelQvao QpvRelQvao dtLvEjection...
        Rmt Rvao Rsys Rtc Rpv Rpul...
        Cao Cvc Cpa Cpu ...
        SBV  Vw_lv Vw_rv f_PassiveForce];
end

save res res

% load k
% 
% filename = ['Data' num2str(k) '.mat' ];
% save(filename);

%XXXXXXXXXXXXXXXXX End of Script_function XXXXXXXXXXXXXXXXXXX

%XXXXXXXXXXXXXXXXX Start of model function XXXXXXXXXXXXXXXXXX

function outputCalculateNextBeat= CalculateNextBeat(iBeat,NbBeats,HR,y0,paramsHemo,IntermediatePlot)

global tVaoOpen;
global tVpvOpen;

lengthbeat=60000/HR;
tstart=(iBeat-1)*lengthbeat;
tend = tstart+lengthbeat;

if IntermediatePlot~=0  ||  iBeat==NbBeats ; dtoutput=lengthbeat/10000;
else dtoutput=lengthbeat;
end;
tspan = tstart:dtoutput:tend; 
options = odeset('RelTol',1e-5,'MaxStep',1);
%options = odeset('RelTol',1e-13,'AbsTol',1e-10);%'MaxStep',1);
% options = odeset('RelTol',1e-6,'MaxStep',1);

%xxxxxxxxxxxxxxxxx Load variations xxxxxxxxxxxxxxx
iBeatLoadVar=paramsHemo(57);    
%if iBeat>=iBeatLoadVar; HR = 85; end; % Preload 
if iBeat>=iBeatLoadVar; Rtc=10*paramsHemo(17); paramsHemo(17)=Rtc; end; % Preload 
% if iBeat>=iBeatLoadVar; Cvc=2*paramsHemo(19); paramsHemo(19)=Cvc; end; % Preload 
% if iBeat>=iBeatLoadVar; Rpul=3.5*paramsHemo(13); paramsHemo(13)=Rpul; end; % Preload 
%if iBeat>=iBeatLoadVar; Rmt=4*paramsHemo(15); paramsHemo(15)=Rmt; end; % Preload 
% if iBeat>=iBeatLoadVar; Rsys=2*paramsHemo(9); paramsHemo(9)=Rsys; end; % Preload 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
tVaoOpen=0;
tVpvOpen=0;

[t,y] = ode15s(@rhs,tspan,y0,options,iBeat,NbBeats,HR,paramsHemo);
yfinal = y(end,:);

outputCalculateNextBeat=yfinal;

if IntermediatePlot==0  &&  iBeat~=NbBeats
    Plv_min=1e69; Plv_max=1e69; Plv_mean=1e69; Plv_ampl=1e69;
    Vlv_min=1e69; Vlv_mean=1e69; SVlv=1e69;
    Pao_min=1e69; Pao_max=1e69; Pao_mean=1e69; Pao_ampl=1e69;
    Vao_mean=1e69; 
    Pvc_min=1e69; Pvc_max=1e69; Pvc_mean=1e69; Pvc_ampl=1e69;
    Vvc_mean=1e69; 
    Prv_min=1e69; Prv_max=1e69; Prv_mean=1e69; Prv_ampl=1e69;
    Vrv_min=1e69; Vrv_mean=1e69; SVrv=1e69;
    Ppa_min=1e69; Ppa_max=1e69; Ppa_mean=1e69; Ppa_ampl=1e69;
    Vpa_mean=1e69; 
    Ppu_min=1e69; Ppu_max=1e69; Ppu_mean=1e69; Ppu_ampl=1e69;
    Vpu_mean=1e69; 
    Qmt_max=1e69;
    Qvao_max=1e69;
    Qtc_max=1e69;
    Qpv_max=1e69;
    dtLvEjection=1e69;
    SVlv=1e69;
    SVrv=1e69;

    outputCalculateNextBeat= [outputCalculateNextBeat...
        Plv_min Plv_max Plv_mean Plv_ampl...
        Vlv_min Vlv_mean SVlv ...
        Pao_min Pao_max Pao_mean Pao_ampl...
        Vao_mean...
        Pvc_min Pvc_max Pvc_mean Pvc_ampl...
        Vvc_mean...
        Prv_min Prv_max Prv_mean Prv_ampl...
        Vrv_min Vrv_mean SVrv...
        Ppa_min Ppa_max Ppa_mean Ppa_ampl...
        Vpa_mean...
        Ppu_min Ppu_max Ppu_mean Ppu_ampl...
        Vpu_mean...
        Qmt_max,Qvao_max,Qtc_max,Qpv_max,dtLvEjection...
        ];
    return;
end

for i=1:length(t)
    temp(i,:)=rhs(t(i),y(i,:),iBeat,NbBeats,HR,paramsHemo,'PostProcessing'); 
end
j=1;
Vla(:,iBeat)=temp(:,j); j=j+1;
Pla(:,iBeat)=temp(:,j); j=j+1;
Vlv(:,iBeat)=temp(:,j); j=j+1;
load Baseline

V_LV = [V_LV temp(:,j-1)];

Plv(:,iBeat)=temp(:,j); j=j+1;

P_LV = [P_LV temp(:,j-1)];
% Q_Inert_lv(:,iBeat)=temp(:,j); j=j+1;
j=j+1;
Vao(:,iBeat)=temp(:,j); j=j+1;
Pao(:,iBeat)=temp(:,j); j=j+1;

P_AO = [P_AO temp(:,j-1)];

Vas(:,iBeat)=temp(:,j); j=j+1;
Pas(:,iBeat)=temp(:,j); j=j+1;
Vvc(:,iBeat)=temp(:,j); j=j+1;
Pvc(:,iBeat)=temp(:,j); j=j+1;
Vrv(:,iBeat)=temp(:,j); j=j+1;

V_RV = [V_RV temp(:,j-1)];

Prv(:,iBeat)=temp(:,j); j=j+1;

P_RV = [P_RV temp(:,j-1)];
Vpa(:,iBeat)=temp(:,j); j=j+1;
Ppa(:,iBeat)=temp(:,j); j=j+1;
Vpu(:,iBeat)=temp(:,j); j=j+1;
Ppu(:,iBeat)=temp(:,j); j=j+1;

P_PU = [P_PU temp(:,j-1)];

if NbBeats == 5
    
save Baseline P_PU P_LV P_RV V_LV V_RV P_AO
    
end

Lm_la(:,iBeat)=temp(:,j); j=j+1;
L_la(:,iBeat)=temp(:,j); j=j+1;
Fm_la(:,iBeat)=temp(:,j); j=j+1;
Pla_parall(:,iBeat)=temp(:,j); j=j+1;
Lm_lv(:,iBeat)=temp(:,j); j=j+1;
L_lv(:,iBeat)=temp(:,j); j=j+1;
Fm_lv(:,iBeat)=temp(:,j); j=j+1;
Fparall_lv(:,iBeat)=temp(:,j); j=j+1;
Plv_parall(:,iBeat)=temp(:,j); j=j+1;
Lm_rv(:,iBeat)=temp(:,j); j=j+1;
L_rv(:,iBeat)=temp(:,j); j=j+1;
Fm_rv(:,iBeat)=temp(:,j); j=j+1;
Prv_parall(:,iBeat)=temp(:,j); j=j+1;
Qin_la(:,iBeat)=temp(:,j); j=j+1;
Qmt(:,iBeat)=temp(:,j); j=j+1;
Qvao(:,iBeat)=temp(:,j); j=j+1;
Qsys(:,iBeat)=temp(:,j); j=j+1;
Qtc(:,iBeat)=temp(:,j); j=j+1;
Qpv(:,iBeat)=temp(:,j); j=j+1;
Qpul(:,iBeat)=temp(:,j); j=j+1;
Irel(:,iBeat)=temp(:,j); j=j+1;
INa(:,iBeat)=temp(:,j); j=j+1;
INCx(:,iBeat)=temp(:,j); j=j+1;
IK1(:,iBeat)=temp(:,j); j=j+1;
ICaL(:,iBeat)=temp(:,j); j=j+1;
INaH(:,iBeat)=temp(:,j); j=j+1;
IKs(:,iBeat)=temp(:,j); j=j+1;
Ileak(:,iBeat)=temp(:,j); j=j+1;
Ilerk(:,iBeat)=temp(:,j); j=j+1;
Iup(:,iBeat)=temp(:,j); j=j+1;

% Hemodynamic parameters
f_PassiveForce=paramsHemo(1);
f_Rref_SARC=paramsHemo(2);
TAU_Bioch_lv=paramsHemo(3);
fGeom_lv=paramsHemo(4);
fGeom_rv=paramsHemo(5);
Inert_lv=paramsHemo(6);
Inert_rv=paramsHemo(7);
SBV=paramsHemo(8);  
Rsys=paramsHemo(9);
Ras=paramsHemo(10);
Cao=paramsHemo(11);
Cas=paramsHemo(12);
Rpul=paramsHemo(13); 
Rprox=paramsHemo(14);
Rmt=paramsHemo(15);
Rvao=paramsHemo(16);
Rtc=paramsHemo(17);
Rpv=paramsHemo(18);
Cvc=paramsHemo(19); 
Cpa=paramsHemo(20);
Cpu=paramsHemo(21);
alfaIN=paramsHemo(22);
alfaIN=paramsHemo(23);
alfaOUT=paramsHemo(24);
alfaOUT=paramsHemo(25);
betaINleft=paramsHemo(26);
betaINright=paramsHemo(27);
betaOUTleft=paramsHemo(28);
betaOUTright=paramsHemo(29);
Gamav_lv=paramsHemo(30); 
Vo_lv=paramsHemo(31); 
Gamav_rv=paramsHemo(32); 
Vo_rv=paramsHemo(33); 
Gamav_la=paramsHemo(34); 
Vo_la=paramsHemo(35); 
Vw_la=paramsHemo(36);
Lr_la=paramsHemo(37);
Vw2refSARC_la=paramsHemo(38);
NSARC_la=paramsHemo(39);
Vw_lv=paramsHemo(40); 
Lr_lv=paramsHemo(41);
Vw2refSARC_lv=paramsHemo(42);
NSARC_lv=paramsHemo(43);
Vw_rv=paramsHemo(44); 
Lr_rv=paramsHemo(45);
Vw2refSARC_rv=paramsHemo(46);
NSARC_rv=paramsHemo(47);    
WantedLm_lv_m=paramsHemo(48);  
WantedLm_lv_M=paramsHemo(49);

save Variables_BL

% load k
% 
% filename = ['Variables' num2str(k) '.mat' ];
% save(filename);

%xxxxxxxxx Mean, max, min values of the functions xxxxxxxxxxxxxx
Pao_mean=mean(Pao(1:end-1,iBeat));
Plv_mean=mean(Plv(1:end-1,iBeat));
Vlv_mean=mean(Vlv(1:end-1,iBeat));
Vao_mean=mean(Vao(1:end-1,iBeat));
Pvc_mean=mean(Pvc(1:end-1,iBeat));
Vvc_mean=mean(Vvc(1:end-1,iBeat));
Prv_mean=mean(Prv(1:end-1,iBeat));
Vrv_mean=mean(Vrv(1:end-1,iBeat));
Ppa_mean=mean(Ppa(1:end-1,iBeat));
Vpa_mean=mean(Vpa(1:end-1,iBeat));
Ppu_mean=mean(Ppu(1:end-1,iBeat));
Vpu_mean=mean(Vpu(1:end-1,iBeat));
Qmtmean=mean(Qmt(1:end-1,iBeat));
CO=60*Qmtmean;
Qtcmean=mean(Qtc(1:end-1,iBeat));

Qmt_max=max(Qmt(1:end-1,iBeat));
Qvao_max=max(Qvao(1:end-1,iBeat));
Qtc_max=max(Qtc(1:end-1,iBeat));
Qpv_max=max(Qpv(1:end-1,iBeat));

LvEjection=find(Qvao(:,iBeat)> 0);
iStartLvEjection=min(LvEjection);
iEndLvEjection=max(LvEjection);

tPositive=dtoutput*(iStartLvEjection-1);
%tStartLvEjection=tPositive-Qvao(iStartLvEjection,iBeat)*dtoutput/(Qvao(iStartLvEjection+1,iBeat)-Qvao(iStartLvEjection,iBeat));

tPositive=dtoutput*(iEndLvEjection-1);
tEndLvEjection=tPositive;%+Qvao(iEndLvEjection,iBeat)*dtoutput/(Qvao(iEndLvEjection-1,iBeat)-Qvao(iEndLvEjection,iBeat));

dtLvEjection=tEndLvEjection;%-tStartLvEjection;



L_la_max=max(L_la(:,iBeat));
L_la_min=min(L_la(:,iBeat));
L_lv_max=max(L_lv(:,iBeat));
L_lv_min=min(L_lv(:,iBeat));
L_rv_max=max(L_rv(:,iBeat));
L_rv_min=min(L_rv(:,iBeat));

Lm_la_max=max(Lm_la(:,iBeat));
Lm_la_min=min(Lm_la(:,iBeat));
Lm_lv_max=max(Lm_lv(:,iBeat));
Lm_lv_min=min(Lm_lv(:,iBeat));
Lm_rv_max=max(Lm_rv(:,iBeat));
Lm_rv_min=min(Lm_rv(:,iBeat));

Vla_max=max(Vla(:,iBeat));
Vlv_max=max(Vlv(:,iBeat));
Vrv_max=max(Vrv(:,iBeat));
Pla_max=max(Pla(:,iBeat));
Plv_max=max(Plv(:,iBeat));
Prv_max=max(Prv(:,iBeat));

Pao_max=max(Pao(:,iBeat));
Pao_min=min(Pao(:,iBeat));
Pao_ampl=Pao_max-Pao_min;

Pvc_max=max(Pvc(:,iBeat));
Pvc_min=min(Pvc(:,iBeat));
Pvc_ampl=Pvc_max-Pvc_min;

Ppa_max=max(Ppa(:,iBeat));
Ppa_min=min(Ppa(:,iBeat));
Ppa_ampl=Ppa_max-Ppa_min;

Ppu_max=max(Ppu(:,iBeat));
Ppu_min=min(Ppu(:,iBeat));
Ppu_ampl=Ppu_max-Ppu_min;

Vla_min=min(Vla(:,iBeat));
Vlv_min=min(Vlv(:,iBeat));
Vrv_min=min(Vrv(:,iBeat));
Pla_min=min(Pla(:,iBeat));
Plv_min=min(Plv(:,iBeat));
Prv_min=min(Prv(:,iBeat));

Plv_ampl=Plv_max-Plv_min;
Prv_ampl=Prv_max-Prv_min;

SVlv=Vlv_max-Vlv_min;
SVrv=Vrv_max-Vrv_min;

outputCalculateNextBeat= [outputCalculateNextBeat...
Plv_min Plv_max Plv_mean Plv_ampl...
Vlv_min Vlv_mean SVlv ...
Pao_min Pao_max Pao_mean Pao_ampl...
Vao_mean...
Pvc_min Pvc_max Pvc_mean Pvc_ampl...
Vvc_mean...
Prv_min Prv_max Prv_mean Prv_ampl...
Vrv_min Vrv_mean SVrv...
Ppa_min Ppa_max Ppa_mean Ppa_ampl...
Vpa_mean...
Ppu_min Ppu_max Ppu_mean Ppu_ampl...
Vpu_mean...
Qmt_max,Qvao_max,Qtc_max,Qpv_max,dtLvEjection...
                           ];

                                           
if IntermediatePlot~=1; return; end
save Variables

