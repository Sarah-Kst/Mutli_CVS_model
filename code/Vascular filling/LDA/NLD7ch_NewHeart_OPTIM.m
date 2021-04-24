function output_NLD7ch_NewHeart_OPTIM = NLD7ch_NewHeart_OPTIM(FactOptim,typerun)
%

%XXXXXXXXXXXXXXX Start of Script_function XXXXXXXXXXXXXXXXXX
% Reference: ten Tusscher and Panfilov. AJP 291:H1088-H1100, 2006.
% Changes: Irel(Shannon),Ileak(Shannon+Rstate),ICaL(Severi,Negroni) 
% Number of variables: 36+6. Units: ms (t), uM ([]), micron (L_lv), mN/mm2 (Fm)
% Last Change: 27-2-2014. 
close all;
% clear all;
% clc; 

nRes=1;IdRes(nRes,1:40)=1;
iBeatLoadVar=1e69;


if typerun==1 || typerun==2 || typerun==4 || typerun==5; % always ecept for identification
%   load IdRes; % read identified factors in file IdRes and use line nRes
%   icolumn=size(IdRes,2)-7;
%   nRes=find(IdRes(:,icolumn)== min(IdRes(:,icolumn)))    % if nRes=-1, all factors are changed to 1
end


if typerun==1; % calcul unique without Load Variation
  IntermediatePlot=1;
  NbBeats = 5;%10;
 % deltaSVmax=0.00025;
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
%   deltaSVmax=1e20;
%   deltaSVmax=0.00008;
  deltaSVmax=0.00015;
end;
if typerun==3;  % identification
   IntermediatePlot=0;
  NbBeats = 15;
  deltaSVmax=0.0001;
end;
if typerun==5; % DataLoops
  IntermediatePlot=-1;
  NbBeats = 50;
%   deltaSVmax=1e20;
%   deltaSVmax=0.00008;
  deltaSVmax=0.00015;
  
end;

% FactOptim=ones(17,1);
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
% disp('                                                          ');
% disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
% disp('Total number of beats:'); disp(NbBeats);
% disp('-------------------'); 


% Corrective factors
% facteur=max(0,1-((iBeat-1)*1.0)/(0.8*NbBeats)); % facteur décroissant à partir de 1 -> 0
f_PassiveForce=fOptimf_PassiveForce;
f_Rref_SARC=1;
TAU_Bioch_lv=1;
fGeom_lv=1.; fGeom_rv=1.; fGeom_la=1;
factGamav=1;


% Inert_lv=0.1;
Inert_lv=0;
Inert_rv=Inert_lv/10;

fR=1.2;   % *1e69
    fRpul=fR;
    fRsys=fR;
    fRvlv=fR;
fc=1;
    fcPUL=fc;
    fcSYS=fc;
fW=1;

VoT=5000; SBV=fOptimSBV*(1145.928421156985-423.93);%*0.5;  %SBV=485;

Rsystot=(600)*fRsys; 
Caotot=1.377*fcSYS; 
%     if (Inert_lv~=0)
%         Fact_CaoCas=0.01;
%         Fact_RsysRas=0.9;
%     else
        Fact_CaoCas=1;
        Fact_RsysRas=1;
%     end
    Rsys=fOptimRsys*Fact_RsysRas*Rsystot*0.989476210397302*1.1448789996563795e+000;
    %Rsys = Rsys*0.8;
    fOptimRas=1;Ras=fOptimRas*(1-Fact_RsysRas)*Rsystot;
    Cao=fOptimCao*Fact_CaoCas*Caotot*6.2871437024222943e-001;
    fOptimCas=1;Cas=fOptimCas*(1-Fact_CaoCas)*Caotot;
% Rprox=fOptimRprox*(5)*fRpul;
fOptimRprox=1;Rprox=fOptimRprox*(1e69)*fRpul;
Rpul=fOptimRpul*140.2926*fRpul; 


% Rvao=fOptimRvao*47.960405664397449*fRvlv;
% Rmt=fOptimRmt*0.460617125237590*47.960405664397449*fRvlv;
% Rtc=fOptimRtc*0.241053599408589*47.960405664397449*fRvlv;
% Rpv=fOptimRpv*0.073192320523130*47.960405664397449*fRvlv;

Rmt=fOptimRmt*17.28243*fRvlv*2.75;
Rvao=fOptimRvao*38.6313*fRvlv;
Rtc=fOptimRtc*24.3987*fRvlv/2.75;
Rpv=fOptimRpv*7.11629*fRvlv;

% Rsys=Rsys+0.8*Rvao;Rvao=0.2*Rvao;

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
% RrefSARC_la=3.1174*f_Rref_SARC;
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
WantedVrv_m=60; %;%61.86;
WantedLm_rv_m=0.93;
WantedLm_rv_M=1.08;
% RrefSARC_rv=3.485*f_Rref_SARC;
%RrefSARC_rv=3.1174*f_Rref_SARC;
RrefSARC_rv=(60*(2*pi)^3/(4/3*pi*(WantedLm_rv_M^3-WantedLm_rv_m^3)))^(1/3)*WantedLm_rv_m/2/pi;
Vw2refSARC_rv=Vsph(RrefSARC_rv)-WantedVrv_m;
NSARC_rv=2*pi*RrefSARC_rv/(WantedLm_rv_m*1e-4);


% % Corrective terms from JacMatrix

% facteur=0.8;
% Rmt=Rmt-6.672907647888428*facteur;
% Rvao=Rvao-0.419787277135136*facteur;
% Rsys=Rsys-12.489836849864201*facteur;
% Rtc=Rtc-0.852694951514154*facteur;
% Rpv=Rpv+0.054093896236126*facteur;
% Rpul=Rpul+0.337471241000097*facteur;
% 
% Cao=Cao+0.012731247891722*facteur;
% Cvc=Cvc+0.305980006686756*facteur;
% Cpa=Cpa+0.077923233664425*facteur;
% Cpu=Cpu+1.834516749706274*facteur;
% 
% SBV=SBV+18.727731428526447*facteur;
% Vw_lv=Vw_lv+29.178914057482995*facteur;
% Vw_rv=Vw_rv+2.790947729749009*facteur;
% 



  
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read initial conditions from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load FinalConditions_NLD7ch_NewHeart_BL    
load FinalConditions_NLD7ch_NewHeart   


% 
% %Loading the y0 vector
% yORIGINAL=[mo ho jo do fo fgo fcao ro so xso xrao xrbo Reo Casro Cajo Caio Kio Naio Vmo];
% poidsORIGINAL=1;


% load ElectrophysioCI    
% yElectrophysio=yfinal;

% if typerun==5; 
%     load FinalConditions_NLD7ch_NewHeart    
% else
%     %load FinalConditions_NLD7ch_NewHeart_BCKP
%     % load FinalConditions_NLD7ch_NewHeart_BaseLine
%     load FinalConditions_NLD7ch_NewHeart 
% end

% yfinal(1:19)=yElectrophysio(1:19);
% yfinal(31)=yElectrophysio(31);
% yfinal(32)=yElectrophysio(32);


% load Copy_of_FinalConditions_NLD7ch_NewHeart_REF4
% load Copy_of_FinalConditions_NLD7ch_NewHeart_Inert    
% load Copy_of_FinalConditions_NLD7ch_NewHeart_NoInert    
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
% Ally0=[ ];
% yrest=[];
% save yrest yrest;
while abs(deltaSV)>deltaSVmax && flag~=flagmax
if flag~=1; NbBeats=10; end
for iBeat=1:NbBeats
    outputCalculateNextBeat=CalculateNextBeat(iBeat,NbBeats,HR,y0,paramsHemo,IntermediatePlot); 
    y0=outputCalculateNextBeat(1:nbODEs);
    
    disp(sprintf('Beat n° %d',iBeat));
%     
%     LgthoutputCalculateNextBeat=size(outputCalculateNextBeat,2);
%     yrest2add=outputCalculateNextBeat(LgthoutputCalculateNextBeat-nbODEs+1:LgthoutputCalculateNextBeat);
%     load yrest;
%     yrest=[yrest;...
%     yrest2add];
%     save yrest yrest;
% 
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
end

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
% if typerun == 2
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
% end;


load k

filename = ['Data' num2str(k) '.mat' ];
save(filename);



 
%XXXXXXXXXXXXXXXXX End of Script_function XXXXXXXXXXXXXXXXXXX

%XXXXXXXXXXXXXXXXX Start of model function XXXXXXXXXXXXXXXXXX

function outputCalculateNextBeat= CalculateNextBeat(iBeat,NbBeats,HR,y0,paramsHemo,IntermediatePlot)

global tVaoOpen;
global tVpvOpen;

lengthbeat=60000/HR;
tstart=(iBeat-1)*lengthbeat;
tend = tstart+lengthbeat;

if IntermediatePlot~=0  ||  iBeat==NbBeats ; dtoutput=lengthbeat/1000;
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
% tic;
tVaoOpen=0;
tVpvOpen=0;
[t,y] = ode15s(@rhs,tspan,y0,options,iBeat,NbBeats,HR,paramsHemo);
% [t,y] = ode45(@rhs,tspan,y0,options,iBeat,NbBeats,HR,paramsHemo);
% toc
% save y y;
yfinal = y(end,:);
Vrest=min(y(:,19));
% jrest=find(y(:,19)==Vrest);
% yrest = y(jrest,:);
% disp(sprintf('      i, Vrest = : %d    %0.16e',iBeat,Vrest));
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
Plv(:,iBeat)=temp(:,j); j=j+1;
% Q_Inert_lv(:,iBeat)=temp(:,j); j=j+1;
j=j+1;
Vao(:,iBeat)=temp(:,j); j=j+1;
Pao(:,iBeat)=temp(:,j); j=j+1;
Vas(:,iBeat)=temp(:,j); j=j+1;
Pas(:,iBeat)=temp(:,j); j=j+1;
Vvc(:,iBeat)=temp(:,j); j=j+1;
Pvc(:,iBeat)=temp(:,j); j=j+1;
Vrv(:,iBeat)=temp(:,j); j=j+1;
Prv(:,iBeat)=temp(:,j); j=j+1;
Vpa(:,iBeat)=temp(:,j); j=j+1;
Ppa(:,iBeat)=temp(:,j); j=j+1;
Vpu(:,iBeat)=temp(:,j); j=j+1;
Ppu(:,iBeat)=temp(:,j); j=j+1;
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

Plv_save=Plv(:,iBeat);
save Plv Plv_save;

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

save Variables

load k

filename = ['Variables' num2str(k) '.mat' ];
save(filename);


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


% 
% 
% 
% dPlv=Plv(iAfterLoad+1,iBeat)-Plv(iAfterLoad,iBeat);
% dPao=Pao(iAfterLoad+1,iBeat)-Pao(iAfterLoad,iBeat);
% x=-(Plv(iAfterLoad+1,iBeat)-Pao(iAfterLoad+1,iBeat))/(dPlv-dPao);
% AfterLoad_lv=Plv(iAfterLoad,iBeat)+x*dPlv;
% PreLoad_lv=Vlv_max;




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

if  IntermediatePlot==-1 % generate DataLoops
load DataLoops;
    
Ke=105000;            %[mN/mm2/um5]
Lz=0.97;              %[um]
Le=10;                %[mN/mm2/um]
Le_lv=Le*f_PassiveForce;
Ke_lv=Ke*f_PassiveForce;
Le_rv=Le*f_PassiveForce;
Ke_rv=Ke*f_PassiveForce;

Fm_max_lv=max(Fm_lv(:,iBeat));
iFm_max=find(Fm_lv(:,iBeat)== Fm_max_lv);
PassiveForceAtFmMax_lv=Ke_lv*(L_lv(iFm_max,iBeat)-Lz)^5+Le_lv*(L_lv(iFm_max,iBeat)-Lz);
ActForceAtFmMax_lv=Fm_max_lv-PassiveForceAtFmMax_lv;

iAfterLoad=min(find(Plv(:,iBeat)> Pao(:,iBeat)))-1;
dPlv=Plv(iAfterLoad+1,iBeat)-Plv(iAfterLoad,iBeat);
dPao=Pao(iAfterLoad+1,iBeat)-Pao(iAfterLoad,iBeat);
x=-(Plv(iAfterLoad+1,iBeat)-Pao(iAfterLoad+1,iBeat))/(dPlv-dPao);
AfterLoad_lv=Plv(iAfterLoad,iBeat)+x*dPlv;
PreLoad_lv=Vlv_max;

Fm_max_rv=max(Fm_rv(:,iBeat));
iFm_max=find(Fm_rv(:,iBeat)== Fm_max_rv);
PassiveForceAtFmMax_rv=Ke_rv*(L_rv(iFm_max,iBeat)-Lz)^5+Le_rv*(L_rv(iFm_max,iBeat)-Lz);
ActForceAtFmMax_rv=Fm_max_rv-PassiveForceAtFmMax_rv;

iAfterLoad=min(find(Prv(:,iBeat)> Ppa(:,iBeat)))-1;
dPrv=Prv(iAfterLoad+1,iBeat)-Prv(iAfterLoad,iBeat);
dPpa=Ppa(iAfterLoad+1,iBeat)-Ppa(iAfterLoad,iBeat);
x=-(Prv(iAfterLoad+1,iBeat)-Ppa(iAfterLoad+1,iBeat))/(dPrv-dPpa);
AfterLoad_rv=Prv(iAfterLoad,iBeat)+x*dPrv;
PreLoad_rv=Vrv_max;

DataLoops=[DataLoops;...
        (SVrv-SVlv)/SVlv...
        Plv_min Plv_max Plv_mean Plv_ampl...
        Vlv_min SVlv ...
        Pao_min Pao_max Pao_mean Pao_ampl...
        Pvc_min Pvc_max Pvc_mean Pvc_ampl...
        Prv_min Prv_max Prv_mean Prv_ampl...
        Vrv_min SVrv...
        Ppa_min Ppa_max Ppa_mean Ppa_ampl...
        Ppu_min Ppu_max Ppu_mean Ppu_ampl...
        Qmt_max,Qvao_max,Qtc_max,Qpv_max,dtLvEjection,...
        Lm_lv_max  Lm_lv_max-Lm_lv_min...
        L_lv_max   L_lv_max-L_lv_min...
        Fm_max_lv  ActForceAtFmMax_lv  PassiveForceAtFmMax_lv... 
        Lm_lv(iFm_max,iBeat) L_lv(iFm_max,iBeat)... 
        PreLoad_lv,AfterLoad_lv...     
        Lm_rv_max  Lm_rv_max-Lm_rv_min...
        L_rv_max   L_rv_max-L_rv_min...
        Fm_max_rv  ActForceAtFmMax_rv  PassiveForceAtFmMax_rv... 
        Lm_rv(iFm_max,iBeat) L_rv(iFm_max,iBeat)... 
        PreLoad_rv,AfterLoad_rv...     
        Rmt Rvao Rsys Rtc Rpv Rpul...
        Cao Cvc Cpa Cpu ...
        SBV  Vw_lv Vw_rv f_PassiveForce];
%save DataLoops DataLoops;
if iBeat==NbBeats
load DataFinalLoops;
DataFinalLoops=[DataFinalLoops;...
        (SVrv-SVlv)/SVlv...
        Plv_min Plv_max Plv_mean Plv_ampl...
        Vlv_min SVlv ...
        Pao_min Pao_max Pao_mean Pao_ampl...
        Pvc_min Pvc_max Pvc_mean Pvc_ampl...
        Prv_min Prv_max Prv_mean Prv_ampl...
        Vrv_min SVrv...
        Ppa_min Ppa_max Ppa_mean Ppa_ampl...
        Ppu_min Ppu_max Ppu_mean Ppu_ampl...
        Qmt_max,Qvao_max,Qtc_max,Qpv_max,dtLvEjection,...
        Lm_lv_max  Lm_lv_max-Lm_lv_min...
        L_lv_max   L_lv_max-L_lv_min...
        Fm_max_lv  ActForceAtFmMax_lv  PassiveForceAtFmMax_lv... 
        Lm_lv(iFm_max,iBeat) L_lv(iFm_max,iBeat)... 
        PreLoad_lv,AfterLoad_lv...     
        Lm_rv_max  Lm_rv_max-Lm_rv_min...
        L_rv_max   L_rv_max-L_rv_min...
        Fm_max_rv  ActForceAtFmMax_rv  PassiveForceAtFmMax_rv... 
        Lm_rv(iFm_max,iBeat) L_rv(iFm_max,iBeat)... 
        PreLoad_rv,AfterLoad_rv...     
        Rmt Rvao Rsys Rtc Rpv Rpul...
        Cao Cvc Cpa Cpu ...
        SBV  Vw_lv Vw_rv f_PassiveForce];
%save DataFinalLoops DataFinalLoops;
end
end


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

                       
% outputCalculateNextBeat= [outputCalculateNextBeat...
% Plv_min Plv_max Plv_mean Plv_ampl...
% Vlv_min SVlv ...
% Pao_min Pao_max Pao_mean Pao_ampl...
% Pvc_min Pvc_max Pvc_mean Pvc_ampl...
% Prv_min Prv_max Prv_mean Prv_ampl...
% Vrv_min SVrv...
% Ppa_min Ppa_max Ppa_mean Ppa_ampl... 
% Ppu_min Ppu_max Ppu_mean Ppu_ampl...
% yrest...
%                            ];

if IntermediatePlot~=1; return; end

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% tBeat=[0:dtoutput:lengthbeat];
% 
% Vmax=max(y(:,19)); Vmin=min(y(:,19)); w=0.1*(Vmax-Vmin)+Vmin;   
% APD=0; ivmax=find(y(:,19)==Vmax); ist=ivmax;
% for i=ivmax:length(t) 
%     if y(i,19)<w && APD==0
%         APD=t(i)-t(ivmax); ist=i; 
%     end
% end
% 
% if iBeat>=max(0,iBeatLoadVar-1)
% figure(47)
% set(gcf,'Position',[1025 35 800 550])
% hold on
% plot(Vlv(:,iBeat),Plv(:,iBeat),'r',Vla(:,iBeat),Pla(:,iBeat),'b',Vrv(:,iBeat),Prv(:,iBeat),'g');
% axis([30, 150, -20, 140]); 
% xlabel ('volume (ml)')
% ylabel ('Plv(r),Pla(b),Prv(b)')
% hold on
% end
% 
% figure(1); set(gcf,'Position',[5 35 1000 550])
% subplot (1,2,1)
% plot(tBeat,y(:,16),'g',tBeat,0.1*Fm_lv(:,iBeat),'r',tBeat,Lm_lv(:,iBeat),'b',tBeat,y(:,19)/75+1.2,'k');
% hold on 
% plot (tBeat,y(:,14)/3000,'-.b',tBeat,y(:,15)/1000,'r',tBeat,L_lv(:,iBeat),'b');
% hold on 
% plot (tBeat,Vlv(:,iBeat)/70,'c',tBeat,Pao(:,iBeat)/70,'c',tBeat,Plv(:,iBeat)/70,'c',tBeat,Pla(:,iBeat)/70,'c',tBeat,Qmt(:,iBeat)/70,'c');
% plot (tBeat,Plv_parall(:,iBeat)/10,'--c')
% axis([0,tBeat(end),-.1,2.5]);
% xlabel ('time (ms)')
% ylabel ('Fm(r), [Ca]i(g), V(b)/75+1.2, [Ca]sr(b), Lm_lv(b), L_lv(b), hp(m)')
% ytext=2.35; desc=0.12; xtext=tBeat(nearest(0.5*end));
% 
% 
% j=-1;
% j=j+1;text(xtext,(ytext-j*desc),['BeatNo= ',num2str(iBeat)])
% j=j+1;text(xtext,(ytext-j*desc),['Fm_{la,M/m}= ',num2str(max(Fm_la(:,iBeat))),' / ',num2str(min(Fm_la(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['Fm_{lv,M/m}= ',num2str(max(Fm_lv(:,iBeat))),' / ',num2str(min(Fm_lv(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['Fm_{rv,M/m}= ',num2str(max(Fm_rv(:,iBeat))),' / ',num2str(min(Fm_rv(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['Lm_{lv,M/m}= ',num2str(Lm_lv_max),' / ',num2str(Lm_lv_min)])
% j=j+1;text(xtext,(ytext-j*desc),['Lm_{rv,M/m}= ',num2str(Lm_rv_max),' / ',num2str(Lm_rv_min)])
% j=j+1;text(xtext,(ytext-j*desc),['Lm_{la,M/m}= ',num2str(max(Lm_la(:,iBeat))),' / ',num2str(min(Lm_la(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['L_{lv,M/m}= ',num2str(L_lv_max),' / ',num2str(L_lv_min)])
% j=j+1;text(xtext,(ytext-j*desc),['L_{rv,M/m}= ',num2str(L_rv_max),' / ',num2str(L_rv_min)])
% j=j+1;text(xtext,(ytext-j*desc),['L_{la,M/m}= ',num2str(max(L_la(:,iBeat))),' / ',num2str(min(L_la(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['[Ca]srmax= ',num2str(max(y(:,14)))])
% j=j+1;text(xtext,(ytext-j*desc),['[Ca]imax= ',num2str(max(y(:,16)))])
% j=j+1;text(xtext,(ytext-j*desc),['[Ca]imin= ',num2str(min(y(:,16)))])
% j=j+1;text(xtext,(ytext-j*desc),['APD90= ',num2str(APD)])
% j=j+1;text(xtext,(ytext-j*desc),['Vrest= ',num2str(min(y(:,19)))])
% j=j+1;text(xtext,(ytext-j*desc),['[Na]i= ',num2str(max(y(:,18)))])
% j=j+1;text(xtext,(ytext-j*desc),['[K]i= ',num2str(max(y(:,17)))])
% j=j+1;text(xtext,(ytext-j*desc),['[Ca]cmax= ',num2str(max(y(:,15)))])
% j=j+1;text(xtext,(ytext-j*desc),['[Ca]cmin= ',num2str(min(y(:,15)))])
% hold off
% 
% subplot(1,2,2)
% 
% v=0:0.5:120; 
% PlaG=Gamav_la*(v-Vo_la).^3; plot(v,PlaG,':b');
% hold on
% PlvG=Gamav_lv*(v-Vo_lv).^3; plot(v,0,'k',v,PlvG,':r');
% hold on
% PrvG=Gamav_rv*(v-Vo_rv).^3; plot(v,PrvG,':g');
% hold on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% plot(v,0,'k',Vlv(:,iBeat),Plv(:,iBeat),'r',Vla(:,iBeat),Pla(:,iBeat),'b',Vrv(:,iBeat),Prv(:,iBeat),'g');
% plot(Vla(:,iBeat),Pla_parall(:,iBeat),'--b');
% plot(Vlv(:,iBeat),Plv_parall(:,iBeat),'--r');
% plot(Vrv(:,iBeat),Prv_parall(:,iBeat),'--g');
% axis([30, 150, -20, 140]); 
% xlabel ('volume (ml)')
% ylabel ('Plv(r),Pla(b),Prv(b)')
% desc=0.05*150; xtext=0.55*120; ytext=150-2*desc;
% j=0;
% j=j+1;text(xtext,(ytext-j*desc),['P_{lv,M/m}= ',num2str(max(Plv(:,iBeat))),' / ',num2str(min(Plv(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['P_{rv,M/m}= ',num2str(max(Prv(:,iBeat))),' / ',num2str(min(Prv(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['V_{lv,m}/V_{rv,m}= ',num2str(Vlv_min),' / ',num2str(Vrv_min)])
% j=j+1;text(xtext,(ytext-j*desc),['Pao_{ampl}= ',num2str(Pao_ampl)])
% j=j+1;text(xtext,(ytext-j*desc),['Pvc_{ampl}= ',num2str(Pvc_ampl)])
% j=j+1;text(xtext,(ytext-j*desc),['Ppa_{ampl}= ',num2str(Ppa_ampl)])
% j=j+1;text(xtext,(ytext-j*desc),['Ppu_{ampl}= ',num2str(Ppu_ampl)])
% j=j+1;text(xtext,(ytext-j*desc),['SV= ',num2str(max(Vlv(:,iBeat))-min(Vlv(:,iBeat)))])
% j=j+1;text(xtext,(ytext-j*desc),['del%(r-l) SV= ',num2str(100*(SVrv-SVlv)/SVlv)])
% j=j+1;text(xtext,(ytext-j*desc),['CO(L/min)= ',num2str(CO)])
% w=max(y(:,27)); EF=(w-min(y(:,27)))/w;
% j=j+1;text(xtext,(ytext-j*desc),['EF(%)= ',num2str(100*EF)])
% 
% hold on
% 
% hold off
% pause (0.01);  
% 
% if iBeat<NbBeats; return; end
% % Continuar=input('Next?(Press Enter): ');
% % differents figures corresponding to the last calculated beat
% 
% figure(12); 
% set(gcf,'Position',[20 560 600 400])
% hold on
% plot(tBeat,Lm_la(:,iBeat),'b')
% plot(tBeat,L_la(:,iBeat),':b')
% plot(tBeat,Lm_lv(:,iBeat),'r')
% plot(tBeat,L_lv(:,iBeat),':r')
% plot(tBeat,Lm_rv(:,iBeat),'g')
% plot(tBeat,L_rv(:,iBeat),':g')
% plot([0 tBeat(end)],[WantedLm_lv_M WantedLm_lv_M],'k')
% plot([0 tBeat(end)],[WantedLm_lv_m WantedLm_lv_m],'k')
% legend('Lm_{la}','L_{la}','Lm_{lv}','L_{lv}','Lm_{rv}','L_{rv}','WtdLm_{lv,M}','WtdLm_{lv,m}')
% hold off
% 
% % figure(93); 
% % hold on
% % plot(tBeat,Lm_lv(:,iBeat),'b')
% % plot(tBeat,L_lv(:,iBeat),'--b')
% % plot(tBeat,Fm_lv(:,iBeat)/2.5,'k')
% % plot(tBeat,Fparall_lv(:,iBeat)/2.5,'--k')
% % plot(tBeat,Plv(:,iBeat)/100,'r')
% % plot(tBeat,Plv_parall(:,iBeat)/100,'--r')
% % plot(tBeat,Vlv(:,iBeat)/120,'m')
% % XBp=y(:,21)+y(:,23);
% % maxXBp=max(XBp)
% % hp=L_lv(:,iBeat)-y(:,24);
% % maxhp=max(hp)
% % plot(tBeat,(XBp)/maxXBp,'g',tBeat,(hp)/maxhp,'--g');
% % plot([0 tBeat(end)],[0 0],'k')
% % legend('Lm_{lv}','L_{lv}','Fm_{lv}/2.5','Fparall_{lv}/2.5','Plv/100','Plv_{parall}/100','Vlv/120','XBp/max','hp/max')
% % hold off
% % 
% % 
% 
% % figure(92); % L-V loops
% % hold on
% % plot(Vlv(:,iBeat),L_lv(:,iBeat),'r')
% % plot(Vlv(:,iBeat),Lm_lv(:,iBeat),'b')
% % legend('L_{lv}-V_{lv} (r)','Lm_{lv}-V_{lv} (b)')
% % hold off
% % 
% 
% % max(Qvao(:,iBeat))
% figure(13); 
% set(gcf,'Position',[650 560 600 400])
% hold on
% plot(tBeat,Pla(:,iBeat),'k')
% plot(tBeat,Plv(:,iBeat),'b')
% plot(tBeat,Pao(:,iBeat),'r')
% plot(tBeat,Pas(:,iBeat),':r')
% plot(tBeat,Pvc(:,iBeat),'g')
% plot(tBeat,Prv(:,iBeat),'c')
% plot(tBeat,Ppa(:,iBeat),'m')
% plot(tBeat,Ppu(:,iBeat),'y')
% plot(tBeat,100*Qvao(:,iBeat)/max(Qvao(:,iBeat)),'-k')
% % plot(tBeat,100*Q_Inert_lv(:,iBeat)/max(Q_Inert_lv(:,iBeat)),'--k')
% legend('Pla','Plv','Pao','Pas','Pvc','Prv','Ppa','Ppu','100 Qvao/max')
% hold off
% 
% figure(83); 
% set(gcf,'Position',[660 570 600 400])
% Qvao_M=max(Qvao(:,iBeat));
% hold on
% plot(tBeat,Qmt(:,iBeat)/Qvao_M,'--k')
% plot(tBeat,Qvao(:,iBeat)/Qvao_M,'k')
% % plot(tBeat,Q_Inert_lv(:,iBeat)/Qvao_M,'b')
% plot(tBeat,Qsys(:,iBeat)/Qvao_M,'--b')
% plot(tBeat,Qtc(:,iBeat)/Qvao_M,'--r')
% plot(tBeat,Qpv(:,iBeat)/Qvao_M,'r')
% plot(tBeat,Qpul(:,iBeat)/Qvao_M,'m')
% plot(tBeat,Qin_la(:,iBeat)/Qvao_M,':k')
% legend('Qmt','Qvao','Qsys','Qtc','Qpv','Qpul','Qin_{la}')
% hold off
% 
% % figure(14); % plot all volumes
% % hold on
% % xx=Vla(:,iBeat);
% % plot(tBeat,xx,'k')
% % xx=xx+Vlv(:,iBeat);
% % plot(tBeat,xx,'b')
% % xx=xx+Vao(:,iBeat);
% % plot(tBeat,xx,'r')
% % xx=xx+Vvc(:,iBeat);
% % plot(tBeat,xx,'g')
% % xx=xx+Vrv(:,iBeat);
% % plot(tBeat,xx,'c')
% % xx=xx+Vpa(:,iBeat);
% % plot(tBeat,xx,'m')
% % xx=xx+Vpu(:,iBeat);
% % plot(tBeat,xx,'y')
% % legend('Vla','Vlv','Vao','Vvc','Vrv','Vpa','Vpu')
% % hold off
% 
% TSt=70;
% hpr=0.006;            %[um]
% hwr=0.0001;           %[um]
% 
% % figure(5); % plot all concentrations etc for LA
% % 
% % set(gcf,'Position',[850 85 800 400])
% % hold on
% % axis([0 tBeat(end) -0.5 1.5]);
% % title('LA')
% % 
% % plot([0 tBeat(end)],[1 1],'-.k')
% % 
% % TSCap=y(:,38);
% % maxTSCap=max(TSCap);
% % plot(tBeat,TSCap/maxTSCap,'r')
% % TSp=y(:,40);
% % maxTSp=max(TSp);
% % plot(tBeat,TSp/maxTSp,'--r')
% % TSCaw=y(:,39);
% % maxTSCaw=max(TSCaw);
% % plot(tBeat,TSCaw/maxTSCaw,'b')
% % TSCa=y(:,37);
% % maxTSCa=max(TSCa);
% % plot(tBeat,TSCa/maxTSCa,'g')
% % hp=L_la(:,iBeat)-y(:,41);
% % maxhp=max(hp);
% % plot(tBeat,hp/hpr,'m')
% % hw=L_la(:,iBeat)-y(:,42);
% % maxhw=max(hw);
% % plot(tBeat,hw/hwr,'--m')
% % maxFm=max(Fm_la(:,iBeat));
% % plot(tBeat,Fm_la(:,iBeat)/maxFm,'k')
% % 
% % maxLm=max(Lm_la(:,iBeat));
% % plot(tBeat,Lm_la(:,iBeat)/maxLm,'--k')
% % 
% % legend('','TSCa*/max','TS*/max','TSCA~/max','TSCa/max','hp/hpr','hw/hwr','Fm/max','Lm/max')
% % ytext=0.2; desc=0.085; xtext=tBeat(nearest(0.6*end));
% % j=-1;
% % j=j+1;text(xtext,(ytext-j*desc),['maxTSCa~= ',num2str(maxTSCaw)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxTSCa*= ',num2str(maxTSCap)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxTS*= ',num2str(maxTSp)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxTSCa= ',num2str(maxTSCa)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxhp= ',num2str(maxhp)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxhw= ',num2str(maxhw)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxFm= ',num2str(maxFm)])
% % j=j+1;text(xtext,(ytext-j*desc),['maxLm= ',num2str(maxLm)])
% % 
% % hold off
% 
% 
% figure(6); % plot all concentrations etc for LV
% 
% set(gcf,'Position',[850 55 800 400])
% hold on
% axis([0 tBeat(end) -0.5 1.5]);
% title('LV')
% Ap=950*3;
% Aw=Ap/5;
% plot([0 tBeat(end)],[1 1],'-.k')
% 
% Cai=y(:,16);
% maxCai=max(Cai);
% plot(tBeat,Cai/maxCai,'--g')
% TSCap=y(:,21);
% maxTSCap=max(TSCap);
% plot(tBeat,TSCap/maxTSCap,'r')
% TSp=y(:,23);
% maxTSp=max(TSp);
% plot(tBeat,TSp/maxTSp,'--r')
% TSCaw=y(:,22);
% maxTSCaw=max(TSCaw);
% plot(tBeat,TSCaw/maxTSCaw,'b')
% TSCa=y(:,20);
% maxTSCa=max(TSCa);
% plot(tBeat,TSCa/maxTSCa,'g')
% hp=L_lv(:,iBeat)-y(:,24);
% maxhp=max(hp);
% plot(tBeat,hp/hpr,'m')
% hw=L_lv(:,iBeat)-y(:,25);
% maxhw=max(hw);
% plot(tBeat,hw/hwr,'--m')
% maxFm=max(Fm_lv(:,iBeat));
% plot(tBeat,Fm_lv(:,iBeat)/maxFm,'k')
% 
% maxLm=max(Lm_lv(:,iBeat));
% plot(tBeat,Lm_lv(:,iBeat)/maxLm,'--k')
% 
% F_Act=Ap.*(TSCap+TSp).*hp+Aw*TSCaw.*hw;
% plot(tBeat,F_Act/maxFm,'--b')
% 
% 
% legend('','Cai/max','TSCa*/max','TS*/max','TSCA~/max','TSCa/max','hp/hpr','hw/hwr','Fm/max','Lm/max','F_{Act}/maxFm')
% ytext=0.2; desc=0.085; xtext=tBeat(nearest(0.6*end));
% j=-1;
% j=j+1;text(xtext,(ytext-j*desc),['maxTSCa~= ',num2str(maxTSCaw)])
% j=j+1;text(xtext,(ytext-j*desc),['maxTSCa*= ',num2str(maxTSCap)])
% j=j+1;text(xtext,(ytext-j*desc),['maxTS*= ',num2str(maxTSp)])
% j=j+1;text(xtext,(ytext-j*desc),['maxTSCa= ',num2str(maxTSCa)])
% j=j+1;text(xtext,(ytext-j*desc),['maxhp= ',num2str(maxhp)])
% j=j+1;text(xtext,(ytext-j*desc),['maxhw= ',num2str(maxhw)])
% j=j+1;text(xtext,(ytext-j*desc),['maxFm= ',num2str(maxFm)])
% j=j+1;text(xtext,(ytext-j*desc),['maxLm= ',num2str(maxLm)])
% 
% hold off
% 
% figure(7); % plot all concentrations etc for RV
% 
% set(gcf,'Position',[20 55 800 400])
% hold on
% axis([0 tBeat(end) -0.5 1.5]);
% title('RV')
% 
% plot([0 tBeat(end)],[1 1],'-.k')
% 
% TSCap=y(:,44);
% maxTSCap=max(TSCap);
% plot(tBeat,TSCap/maxTSCap,'r')
% TSp=y(:,46);
% maxTSp=max(TSp);
% plot(tBeat,TSp/maxTSp,'--r')
% TSCaw=y(:,45);
% maxTSCaw=max(TSCaw);
% plot(tBeat,TSCaw/maxTSCaw,'b')
% TSCa=y(:,43);
% maxTSCa=max(TSCa);
% plot(tBeat,TSCa/maxTSCa,'g')
% hp=L_rv(:,iBeat)-y(:,47);
% maxhp=max(hp);
% plot(tBeat,hp/hpr,'m')
% hw=L_rv(:,iBeat)-y(:,48);
% maxhw=max(hw);
% plot(tBeat,hw/hwr,'--m')
% maxFm=max(Fm_rv(:,iBeat));
% plot(tBeat,Fm_rv(:,iBeat)/maxFm,'k')
% 
% maxLm=max(Lm_rv(:,iBeat));
% plot(tBeat,Lm_rv(:,iBeat)/maxLm,'--k')
% 
% legend('','TSCa*/max','TS*/max','TSCA~/max','TSCa/max','hp/hpr','hw/hwr','Fm/max','Lm/max')
% ytext=0.2; desc=0.085; xtext=tBeat(nearest(0.6*end));
% j=-1;j=j+1;text(xtext,(ytext-j*desc),['maxTSCa~= ',num2str(maxTSCaw)])
% j=j+1;text(xtext,(ytext-j*desc),['maxTSCa*= ',num2str(maxTSCap)])
% j=j+1;text(xtext,(ytext-j*desc),['maxTS*= ',num2str(maxTSp)])
% j=j+1;text(xtext,(ytext-j*desc),['maxTSCa= ',num2str(maxTSCa)])
% j=j+1;text(xtext,(ytext-j*desc),['maxhp= ',num2str(maxhp)])
% j=j+1;text(xtext,(ytext-j*desc),['maxhw= ',num2str(maxhw)])
% j=j+1;text(xtext,(ytext-j*desc),['maxFm= ',num2str(maxFm)])
% j=j+1;text(xtext,(ytext-j*desc),['maxLm= ',num2str(maxLm)])
% hold off

save Variables