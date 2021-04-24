function output = rhs(t,y,iBeat,NbBeats,HR,paramsHemo,runType)

ydot = zeros(size(y));
t = mod(t,60000/HR);

MagicFactor=1.;

%% Model Parameters
% Physical Constants
R = 8314.472;%8314;       % [J/kmol*K] 
Frdy =96485.3415;% 96500;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp; %(EX in previous)
Cmem = 1.85e-4;   % [microF] membrane capacitance
Volcel = 16404e-15;% Myoplasmic volume (Liters)
MV=1; NV=(1094e-15)/Volcel; CV=(54.68e-15)/Volcel;% Fractions of Volcel
ICON=Cmem/Frdy/Volcel/1000;

% External ion concentrations     
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
CabL = 1.15;  % Extracellular blood Ca  [mM]
Cao = 1.74*CabL;  % Extracellular Ca  [mM] 

% Hemodynamic parameters
f_PassiveForce=paramsHemo(1);
f_Rref_SARC=paramsHemo(2);
TAU_Bioch_lv=paramsHemo(3);
fGeom_lv=paramsHemo(4);
fGeom_rv=paramsHemo(5);
Inert_lv=paramsHemo(6);
Inert_rv=paramsHemo(7);
Inert_ao=0;
SBV=paramsHemo(8) ;
Rsys=paramsHemo(9);
Ras=paramsHemo(10);
Caor=paramsHemo(11);
Cas=paramsHemo(12);
Rpul=paramsHemo(13);
Rprox=paramsHemo(14);
Rmt=paramsHemo(15);
Rvaor=paramsHemo(16);
Rtc=paramsHemo(17);
Rpv=paramsHemo(18);
Cvc=paramsHemo(19) ;
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
Vw_lv=paramsHemo(40)/MagicFactor;
Lr_lv=paramsHemo(41);
Vw2refSARC_lv=paramsHemo(42);
NSARC_lv=paramsHemo(43);
Vw_rv=paramsHemo(44); 
Lr_rv=paramsHemo(45);
Vw2refSARC_rv=paramsHemo(46);
NSARC_rv=paramsHemo(47);    

% Mechanical parameters
fYv_la=1; ffa_la=1; fYv_lv=1; ffa_lv=1; fYv_rv=1; ffa_rv=1;fgama_lv=1;fgama_rv=1;fgama_la=1;

nc=3;                 %Addim.(=1,2,3...)
Ap=950*nc;            %[mN/mm2/um/uM]; 1400*nc;
Aw=Ap/5;              %[mN/mm2/um/uM]
alfa=0.5*f_PassiveForce;             %[mN/mm2]
bet=80*1;               %[1/um]
Bp=0.5*3.5;               %[1/ms]
Bw=0.35*3.5;              %[1/ms]
Bw_la=Bw;            %[1/ms]
Bw_lv=Bw;          %[1/ms]
Bw_rv=Bw;            %[1/ms]
f=0.0023;             %[1/ms]
gama=28000; %/1.2;           %[1/um2]
hpr=0.006;            %[um]
hwr=0.0001;           %[um]
Ke=105000;            %[mN/mm2/um5]
Lz=0.97;              %[um]
La=1.15;%1.15; %1.125;%1.15;              %[um]
Lc=1.05;              %[um]
Le=10;                %[mN/mm2/um]
RLa=15;%/1.2;               %[1/um2]
TSt=70/nc;            %[mM]
Yb=0.1816;            %[1/uM^nc/ms]
Yc=4;                 %[1/um]
Yd=0.028;             %[1/ms]
Yp=0.1397;            %[1/ms]  
Yr=0.1397;            %[1/ms]
Yv=0.9;               %[1/ms]
Za=0.0023;            %[1/ms]
Zb=0.1397;            %[1/ms]
Zp=0.2095;            %[1/ms]
Zr=7.2626;            %[1/uM^nc/ms]
Fh=200;              %[ms]

Le_lv=Le*f_PassiveForce*MagicFactor;
% Le_lv=Le*MagicFactor;
Ke_lv=Ke*f_PassiveForce*MagicFactor;
% alfa_lv=alfa*f_PassiveForce;
alfa_lv=alfa;
bet_lv=bet;

Le_rv=Le*f_PassiveForce;
% Le_rv=Le;
Ke_rv=Ke*f_PassiveForce;
% alfa_rv=alfa*f_PassiveForce;
alfa_rv=alfa;
bet_rv=bet;

Le_la=Le*f_PassiveForce;
% Le_la=Le;
Ke_la=Ke*f_PassiveForce;
% alfa_la=alfa*f_PassiveForce;
alfa_la=alfa;
bet_la=bet;

GNa = 14.838;    % INa  [nS/pF]
PNaK = 2.724;    % INaK [pA/pF]
KmK = 1;         % INaK [uM]
KmN = 40000;     % INaK [uM]

gam = 0.35;      % INaCa [addim]
KmNa = 87.5;     % INaCa [mM]
KmCa = 1.38;     % INaCa [mM] 
ksat = 0.1;      % INaCa [addim]
kNC = 1000;      % INaCa [pA/pF]

Vup = 6.375;%*0.9;     % Iup [uM/ms]
Kup = 0.25;      % Iup [uM]
Vleak = 0.00036; % Ileak [1/ms] 

Gto = 0.294;     % Ito [nS/pF]
GK1 = 5.405;     % IK1 [nS/pF]
GKr = 0.153;     % IKr [nS/pF]
GKs = 0.098;     % IKs  [nS/pF]
pKN = 0.03;      % IKs [addim]

GCaL = 0.0000398; % ICaL [cm/ms/uF]

GKp = 0.0146;     % IKp  [nS/pF]

GCab = 0.000592;  % ICab [nS/pF]
GCap = 0.1238;    % ICap [nS/pF]
KCap = 0.5;       % ICap [uM]
GNab = 0.00029;   % INab [nS/pF]

Vrel = 0.102;%*0.9;     % Irel [1/ms]
k1p = 0.15e-6;    % Irel [1/uM^2/ms]
k2p = 0.000045;   % Irel [1/uM/ms] 
k3 = 0.06;        % Irel [1/ms]
k4 = 0.005;       % Irel [1/ms]
EC = 1500;        % Irel [uM]
maxsr = 2.5;      % Irel [addim] 
minsr = 1;        % Irel [addim]
Vlc = Vleak;  % Irel [1/ms]

Vxfer = 0.0038;   % Ixfer [1/ms]

Bufi = 130;%200;%175;  % Buffer(i) [uM]; 105(VHHSS2010)
Kbufi = 1;        % Buffer(i) [uM]
Bufsr = 10000;    % Buffer(sr) [uM]
Kbufsr = 300;     % Buffer(sr) [uM] 
Bufj = 400;       % Buffer(j) [uM]
Kbufj = 0.25;     % Buffer(j) [uM]

Vnh = 4.5;        % INaH [uM] %(Crampin: 3.1) 
Phm = 0.0121;     % INaH [1/ms]
Phs = 0.00329;    % INaH [1/ms] 
Pnm = 0.000733;   % INaH [1/ms]
Pns = 0.00269;    % INaH [1/ms]
Kh = 0.1648;      % INaH [uM]
Kna = 33600;      % INaH [uM]
pHo = 7.4;
ALo = 10^(-pHo)*10^6/Kh;    % INaH [addim]
BTo = 1000*Nao/Kna;         % INaH [addim]
pHr=7.15; 
App=9.5;


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%---------Myocyte adjust---------
Vlc=0.05*Vleak; Vleak=0.95*Vleak; k6=1; Khf=1; pH=pHr;
% Kup=0.8*Kup; Gto=1.4*Gto; GNab=0.4*GNab;  
% Bufi=1.25*Bufi; 
%---------Multiscale adjust---------
% Ap=1.*Ap; Aw=Ap/5; bet=1*bet;
% Ap=2.*Ap; Aw=Ap/5; bet=0.8*bet;
f_HeartFailure=1;  % 1 for No HF 
%xxxxxxxxxxxxxxxxx Heart faillure xxxxxxxxxxxxxxx
% GK1=0.51*GK1; Gto=0.64*Gto; kNC=2*kNC; Vup=0.76*Vup; f_HeartFailure=0.7;
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% f_HeartFailure=0.7;
% u=1.8;
% GK1=(1+u*(0.51-1))*GK1; Gto=(1+u*(0.64-1))*Gto; 
% kNC=(1+u*(2-1))*kNC; Vup=(1+u*(0.76-1))*Vup;

% xxxxxxxxxxxxxxx Computation xxxxxxxxxxxxxxxxxx

% Simulation type
% App=0;
dtStim=6;
protocol = 'pace1';
switch lower(protocol)
    case {'none',''},
        I_Fh = 0;
    case 'pace1',        % Stabilization 
		if t>Fh && t<Fh+dtStim
            I_app = App; %pH=pHr;
%  		elseif t>Fh-dtStim && t<Fh
%             I_app = -2.7*App*(y(19)<0); %pH=pHr;   
        else
            I_app = 0.0; %pH=pHr;
        end
		 
     case 'pace3',        % Cafeine pulse
		Vup = 0.072*Vup; pH=pHr; Stac=80;
        if iBeat<Stac
            if t>10 && t<16
                I_app = App;
            else
                I_app = 0.0; 
            end
        elseif iBeat==Stac
            I_app = 0.0;  
            if t>50 && t<900
                block = 0.9; Vup = 0.1*Vup; 
            end
        else    
        end
end

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Nernst Potentials
ENa = (1/FoRT)*log(1000*Nao/y(18));    % [mV]
EK = (1/FoRT)*log(1000*Ko/y(17));	   % [mV]
ECa = (1/FoRT/2)*log(1000*Cao/y(16));  % [mV]

% INa: Fast Na Current
minf=1/(1+exp((-56.86-y(19))/9.03))^2;
AL1=1/(1+exp((-60-y(19))/5));BT1=0.1/(1+exp((35+y(19))/5))+0.1/(1+exp((-50+y(19))/200));
tam=AL1*BT1;

hinf=1/(1+exp((71.55+y(19))/7.43))^2;
AL2=0.057*exp((-80-y(19))/6.8); BT2=2.7*exp(0.079*y(19))+310000*exp(0.3485*y(19));

jinf=1/(1+exp((71.55+y(19))/7.43))^2;
w=-25428*exp(.2444*y(19))-6.948e-06*exp(-.04391*y(19));
AL3=w*(y(19)+37.78)/(1+exp(.311*(y(19)+79.23))); 
BT3=0.02424*exp(-0.01052*y(19))/(1+exp(-0.1378*(y(19)+40.14))); % 0.02424 au lieu de 0.0242
if y(19)>=-40
    AL2=0;AL3=0;
    BT2=0.77/(0.13*(1+exp((y(19)+10.66)/-11.1)));
    BT3=0.6*exp(0.057*y(19))/(1+exp(-.1*(y(19)+32)));
end
tah=1/(AL2+BT2);
taj=1/(AL3+BT3);
ydot(1) = (minf-y(1))/tam;
ydot(2) = (hinf-y(2))/tah;
ydot(3) = (jinf-y(3))/taj;
INa = GNa*y(1)^3*y(2)*y(3)*(y(19)-ENa);

% ICaL: L_lv-type Calcium Current
d1=1/(1+exp((-8-y(19))/7.5)); AL1=1.4/(1+exp((-35-y(19))/13))+.25; 
BT1=1.4/(1+exp((5+y(19))/5)); GA1=1/(1+exp((50-y(19))/20));
tad=AL1*BT1+GA1;
ydot(4)=(d1-y(4))/tad;
              
f1=1/(1+exp((20+y(19))/7)); AL2=1102.5*exp(-((y(19)+27)/15)^2); 
BT2=200/(1+exp((13-y(19))/10)); taf=AL2+BT2+20+180/(1+exp((30+y(19))/10));
ydot(5)=(f1-y(5))/taf;
       
% g1=0.3+0.7/(1+exp((50+y(19))/7)); 
g1=0.33+0.67/(1+exp((35+y(19))/7)); 
AL3=562*exp(-((y(19)+27)^2)/240); 
BT3=31/(1+exp((25-y(19))/10)); tag=AL3+BT3+80/(1+exp((30+y(19))/10)); 
ydot(6)=(g1-y(6))/tag;

% fc1=.1+.9/(1+exp((y(15)/1000-1.5)/.15)); tac=1+80/(1+(y(15)/50)^2);
% u=1; if (fc1>y(7)) && y(19)>-60; u=0; end; %(60)
% ydot(7)=u*(fc1-y(7))/tac; 

fc1=0.6/(1.0+(y(15)/50)^2.0)+0.4; tac=2+80/(1+(y(15)/50)^2);
ydot(7)=(fc1-y(7))/tac; 


w=exp(2*FoRT*(y(19)-15));
IM1=GCaL*4*FoRT*Frdy*(y(19)-15)*(.25*y(15)*w/1000-Cao)/(w-1);
ICaL=y(4)*y(5)*y(6)*y(7)*IM1;

% Ito: Transient Outward K Current (epicardial and M cells)
B1=1/(1+exp((20-y(19))/6)); BT1=9.5*exp(-(y(19)+40)^2/1800)+0.8;
G1=1/(1+exp((20+y(19))/5)); 
BT2=85*exp(-(y(19)+45)^2/320)+5/(1+exp((y(19)-20)/5))+3;
ydot(8)=(B1-y(8))/BT1;
ydot(9)=(G1-y(9))/BT2;       
Ito=Gto*y(8)*y(9)*(y(19)-EK);

% IKs: Slow delayed rectifier current
D1=1/(1+exp((-5-y(19))/14)); u=sqrt(1+exp((5-y(19))/6));
AL1=1400/u; BT1=1/(1+exp((-35+y(19))/15));
TA1=AL1*BT1+80; 
ydot(10)=(D1-y(10))/TA1;       
EKs=log(1000*(Ko+pKN*Nao)/(y(17)+pKN*y(18)))/FoRT;
IKs=GKs*y(10)^2*(y(19)-EKs); 

% IKr: Rapid delayed rectifier current
Q1=1/(1+exp((-26-y(19))/7)); S1=1/(1+exp((88+y(19))/24));
AL1=450/(1+exp((-45-y(19))/10)); BT1=6/(1+exp((30+y(19))/11.5));
AL2=3/(1+exp((-60-y(19))/20)); BT2=1.12/(1+exp((-60+y(19))/20));
TA1=AL1*BT1; TA2=AL2*BT2;
ydot(11)=(Q1-y(11))/TA1;
ydot(12)=(S1-y(12))/TA2; w=(Ko/5.4)^.5;      
IKr=GKr*w*y(11)*y(12)*(y(19)-EK);

% IK1: Inward rectifier K+ current
w=y(19)-EK;
AL1=0.1/(1+exp(0.06*(w-200)));
w1=3*exp(.0002*(w+100))+exp(.1*(w-10)); BT1=w1/(1+exp(-.5*w));
IK1=GK1*AL1/(AL1+BT1)*(y(19)-EK);

% INCx: Na/Ca Exchanger flux
LRK=kNC/(KmNa^3+Nao^3)/(KmCa+Cao);
ALFA=(gam-1)*FoRT*y(19); BETA=gam*FoRT*y(19); 
PFA=(y(18)/1000)^3*Cao*exp(BETA)-Nao^3*y(16)*2.5/1000*exp(ALFA);
INCx=PFA*LRK/(1+ksat*exp(ALFA)); 

% INaK: Na/K pump current
w=1+.1245*exp(-.1*FoRT*y(19))+.0353*exp(-FoRT*y(19));
INaK=PNaK*Ko*y(18)/(KmK+Ko)/(KmN+y(18))/w;        

% ICap: Sarcolemmal Ca pump current
ICap=GCap*y(16)/(KCap+y(16));

% IKp: Sarcolemmal K pump current
IKp=GKp*(y(19)-EK)/(1+exp((25-y(19))/5.98));

% Background currents
ICab=GCab*(y(19)-ECa);
INab=GNab*(y(19)-ENa);

% Iup, Ileak, Ixfer 
Ileak=Vleak*(y(14)-y(16)); Iup=Vup/(1+(Kup/y(16))^2);  
Ixfer=Vxfer*(y(15)-y(16)); 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% IREL: CICR current
u=maxsr-(maxsr-minsr)/(1+(EC/y(14))^2);  
k1=k1p/u; k2=k2p*u; O=y(31); Rr=y(13);      
RI=1-Rr-O-y(32);  
ydot(13)=(k4*RI-k2*Rr*y(15))-k6*(k1*Rr*y(15)^2-k3*O);    %   R
ydot(31)=k6*(k1*Rr*y(15)^2-k3*O)-(k2*O*y(15)-k4*y(32));  %   O
ydot(32)=(k2*O*y(15)-k4*y(32))-(k3*y(32)-k1*RI*y(15)^2); %   I
% ---------------------------------------------
Irel=(Vrel*O+Vlc*Rr)*(y(14)-y(15));
Ilerk=Vlc*Rr*(y(14)-y(15));

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% INaH
% ALi=(10^(-pH))*10^6/Kh; BTi=y(18)/Kna;  
% w=(Phs*ALo+Pnm*BTo)*(1+ALi+BTi)+(Phm*ALi+Pns*BTi)*(1+ALo+BTo);
% INaH=Vnh*(Phm*Pnm*ALi*BTo-Phs*Pns*ALo*BTi)/w; 
INaH=0.; 

% Ion concentrations in the intracellular compartments
% Caid=-(ICab+ICap-2*INCx)*ICON/MV/2+(Ileak-Iup)*NV/MV+Ixfer-nc*(ydot(20)+ydot(21)+ydot(22));
% ydot(16)=Caid/(1+Bufi*Kbufi/(y(16)+Kbufi)^2);
% 
% 
% ydot(14)=(Iup-Ileak-Irel)/(1+Bufsr*Kbufsr/(y(14)+Kbufsr)^2);
% 
% Cajd=-ICaL*ICON/CV/2+Irel*NV/CV-Ixfer*MV/CV;
% ydot(15)=Cajd/(1+Bufj*Kbufj/(y(15)+Kbufj)^2);
% 
% ydot(17)=-(-2*INaK+IKs+IKr+IK1+Ito+IKp-1*I_app)*ICON/MV;
% ydot(18)=-(3*INaK+3*INCx+INa+INab)*ICON/MV + INaH;


%xxxxxxxxxxxxxxxx Force development Left Ventricule xxxxxxxxxxxxxxxxx
% Lm_lv=1.05;
rint_lv=Rsph(y(27));
rout_lv=Rsph(y(27)+Vw_lv);
CurrentRrefSARC_lv=Rsph(y(27)+Vw2refSARC_lv); % radial position of the part of the muscle which was at rmwr for the reference distension
Lm_lv=1.e4*2*pi*CurrentRrefSARC_lv/NSARC_lv;  % NewHeartHeart  

% EpsilonL=1e-5;
% 
% con=0; corL=10; L_lv=0.9884*Lm_lv;
% while (abs(corL)>EpsilonL)
%    Fm_lv=MagicFactor*alfa_lv*(exp(bet_lv*(Lm_lv-L_lv))-1); 
%    FB=MagicFactor*(Ap*(y(21)+y(23))*(L_lv-y(24))+Aw*y(22)*(L_lv-y(25)));
%    w=FB+Ke_lv*(L_lv-Lz)^5+Le_lv*(L_lv-Lz)-Fm_lv;
%    w1=Ap*(y(21)+y(23))+Aw*y(22)+5*Ke_lv*(L_lv-Lz)^4+Le_lv+bet_lv*(Fm_lv+alfa_lv);
%    corL=-w/w1; L_lv=L_lv+0.1*corL; 
%    con=con+1; if con>200; break; end;
% end; %while corL;   

L_lv0=0.9884*Lm_lv;
% optionsfzero = optimset('TolX',1.e-14);
% L_lv = fzero(@(L_lv)Ap*(y(21)+y(23))*(L_lv-y(24))+Aw*y(22)*(L_lv-y(25))+Ke_lv*(L_lv-Lz)^5+Le_lv*(L_lv-Lz)-alfa_lv*(exp(bet_lv*(Lm_lv-L_lv))-1),L_lv0,optionsfzero);
L_lv = fzero(@(L_lv)Ap*(y(21)+y(23))*(L_lv-y(24))+Aw*y(22)*(L_lv-y(25))+Ke_lv*(L_lv-Lz)^5+Le_lv*(L_lv-Lz)-alfa_lv*(exp(bet_lv*(Lm_lv-L_lv))-1),L_lv0);

Fm_lv=alfa_lv*(exp(bet_lv*(Lm_lv-L_lv))-1); 
Fparall_lv=Ke_lv*(L_lv-Lz)^5+Le_lv*(L_lv-Lz);
% F_act=(Ap*(y(21)+y(23))*(L_lv-y(24))+Aw*y(22)*(L_lv-y(25)));
% if F_act > Fm_lv;
%     disp(sprintf('F_act > Fm_lv   %d %d %d',F_act,Fm_lv,L_lv)); 
% end

%L_inter_lv(t)

TS=TSt-y(20)-y(21)-y(22)-y(23); ER=exp(-RLa*(L_lv-La)^2); 
Yh=fYv_lv*Yv*(1-exp(-gama*fgama_lv*(L_lv-y(25)-hwr)^2)); if (L_lv-y(25))>hwr; Yh=0.1*Yh; end;
fa=f*ER*ffa_lv; ga=Za+Yh; gd=Yd*exp(-Yc*(L_lv-Lc)); 
ydot(23)=TAU_Bioch_lv*(-gd*y(23)+Yr*y(21)-Zr*y(23)*y(16)^nc);%---------------TS*
ydot(21)=TAU_Bioch_lv*(Zr*y(23)*y(16)^nc-Yr*y(21)+Yp*y(22)-Zp*y(21));%-------TSCa*
ydot(22)=TAU_Bioch_lv*(Zp*y(21)-Yp*y(22)+fa*y(20)-ga*y(22));%----------------TSCa#
ydot(20)=TAU_Bioch_lv*(ga*y(22)-fa*y(20)+Yb*TS*y(16)^nc-Zb*y(20));%----------TSCa
ydot(24)=TAU_Bioch_lv*(Bp*(L_lv-y(24)-hpr));
ydot(25)=TAU_Bioch_lv*(Bw_lv*(L_lv-y(25)-hwr));


%xxxxxxxxxxxxxxxx Force development Right ventricule xxxxxxxxxxxxxxxxx
rint_rv=Rsph(y(33));
rout_rv=Rsph(y(33)+Vw_rv);
CurrentRrefSARC_rv=Rsph(y(33)+Vw2refSARC_rv); % radial position of the part of the muscle which was at rmwr for the reference distension
Lm_rv=1.e4*2*pi*CurrentRrefSARC_rv/NSARC_rv;  % NewHeartHeart  

% Lm_rv=((y(33)+.4*OLD_Vw_rv)/Kv_rv)^(1/3);
% disp(sprintf('Lm_rv1, Lm_rv =   %d %d',Lm_rv1,Lm_rv)); 
% 
% con=0; corL=10; L_rv=0.9881*Lm_rv;
% while (abs(corL)>EpsilonL)
%    Fm_rv=alfa_rv*(exp(bet_rv*(Lm_rv-L_rv))-1); 
%    FB=Ap*(y(44)+y(46))*(L_rv-y(47))+Aw*y(45)*(L_rv-y(48));
%    w=FB+Ke_rv*(L_rv-Lz)^5+Le_rv*(L_rv-Lz)-Fm_rv;  
%    w1=Ap*(y(44)+y(46))+Aw*y(45)+5*Ke_rv*(L_rv-Lz)^4+Le_rv+bet*(Fm_rv+alfa_rv);  
%    corL=-w/w1; L_rv=L_rv+0.1*corL; 
%    con=con+1; if con>200; break; end;
% end; %while corL;   

% % L_rv0=L_rv;
L_rv0=0.98*Lm_rv;
% L_rv = fzero(@(L_rv)Ap*(y(44)+y(46))*(L_rv-y(47))+Aw*y(45)*(L_rv-y(48))+Ke_rv*(L_rv-Lz)^5+Le_rv*(L_rv-Lz)-alfa_rv*(exp(bet_rv*(Lm_rv-L_rv))-1),L_rv0,optionsfzero);
L_rv = fzero(@(L_rv)Ap*(y(44)+y(46))*(L_rv-y(47))+Aw*y(45)*(L_rv-y(48))+Ke_rv*(L_rv-Lz)^5+Le_rv*(L_rv-Lz)-alfa_rv*(exp(bet_rv*(Lm_rv-L_rv))-1),L_rv0);

Fm_rv=alfa_rv*(exp(bet_rv*(Lm_rv-L_rv))-1); 
Fparall_rv=Ke_rv*(L_rv-Lz)^5+Le_rv*(L_rv-Lz);

%L_inter_rv(t)

TS_rv=TSt-y(43)-y(44)-y(45)-y(46); ER=exp(-RLa*(L_rv-La)^2); 
Yh=fYv_rv*Yv*(1-exp(-gama*fgama_rv*(L_rv-y(48)-hwr)^2)); if (L_rv-y(48))>hwr; Yh=0.1*Yh; end;
fa=f*ER*ffa_rv; ga=Za+Yh; gd=Yd*exp(-Yc*(L_rv-Lc)); 
ydot(46)=-gd*y(46)+Yr*y(44)-Zr*y(46)*y(16)^nc;%---------------TS*_rv
ydot(44)=Zr*y(46)*y(16)^nc-Yr*y(44)+Yp*y(45)-Zp*y(44);%-------TSCa*_rv
ydot(45)=Zp*y(44)-Yp*y(45)+fa*y(43)-ga*y(45);%----------------TSCa#_rv
ydot(43)=ga*y(45)-fa*y(43)+Yb*TS_rv*y(16)^nc-Zb*y(43);%----------TSCa_rv
ydot(47)=Bp*(L_rv-y(47)-hpr);%----------Xp_rv
ydot(48)=Bw_rv*(L_rv-y(48)-hwr);%----------Xw_rv


% 
% % Ion concentrations in the intracellular compartments
Caid=-(ICab+ICap-2*INCx)*ICON/MV/2+(Ileak-Iup)*NV/MV+Ixfer...
    -0.5*nc*(ydot(20)+ydot(21)+ydot(22)+ydot(43)+ydot(44)+ydot(45)); % mean contribution of the lv and rv
ydot(16)=Caid/(1+Bufi*Kbufi/(y(16)+Kbufi)^2);

ydot(14)=(Iup-Ileak-Irel)/(1+Bufsr*Kbufsr/(y(14)+Kbufsr)^2);

Cajd=-ICaL*ICON/CV/2+Irel*NV/CV-Ixfer*MV/CV;
ydot(15)=Cajd/(1+Bufj*Kbufj/(y(15)+Kbufj)^2);

ydot(17)=-(-2*INaK+IKs+IKr+IK1+Ito+IKp-I_app)*ICON/MV;
ydot(18)=-(3*INaK+3*INCx+INa+INab)*ICON/MV + INaH;


%xxxxxxxxxxxxxxxx Force development Left Atrium  xxxxxxxxxxxxxxxxx
rint_la=Rsph(y(26));
rout_la=Rsph(y(26)+Vw_la);
CurrentRrefSARC_la=Rsph(y(26)+Vw2refSARC_la); % radial position of the part of the muscle which was at rmwr for the reference distension
% Lm_la=1.e4*2*pi*CurrentRrefSARC_la/NSARC_la;  % NewHeartHeart  
% 
% con=0; corL=10; L_la=1.03;
% while (abs(corL)>.00001)
%    Fm_la=alfa_la*(exp(bet_la*(Lm_la-L_la))-1); 
%    FB=Ap*(y(38)+y(40))*(L_la-y(41))+Aw*y(39)*(L_la-y(42));
%    w=FB+Ke_la*(L_la-Lz)^5+Le_la*(L_la-Lz)-Fm_la;
%    w1=Ap*(y(38)+y(40))+Aw*y(39)+5*Ke_la*(L_la-Lz)^4+Le_la+bet_la*(Fm_la+alfa_la);
%    corL=-w/w1; L_la=L_la+0.1*corL; 
%    con=con+1; if con>200; break; end;
% end; %while corL;   
 
L_la=1;
Lm_la=1;
Fm_la=0;
Fparall_la=0;
% Fparall_la=Ke_la*(L_la-Lz)^5+Le_la*(L_la-Lz);

% L_la = fzero(@(L_la)Ap*(y(38)+y(40))*(L_la-y(41))+Aw*y(39)*(L_la-y(42))+Ke_la*(L_la-Lz)^5+Le_la*(L_la-Lz)-alfa_la*(exp(bet_la*(Lm_la-L_la))-1),1.03);
% Fm_la=alfa_la*(exp(bet_la*(Lm_la-L_la))-1); 

TS_la=TSt-y(37)-y(38)-y(39)-y(40); ER=exp(-RLa*(L_la-La)^2); 
Yh=fYv_la*Yv*(1-exp(-gama*fgama_la*(L_la-y(42)-hwr)^2)); if (L_la-y(42))>hwr; Yh=0.1*Yh; end;
fa=f*ER*ffa_la; ga=Za+Yh; gd=Yd*exp(-Yc*(L_la-Lc)); 
Cai_la=CaDrLA2(t,60000/HR,f_HeartFailure);
ydot(40)=-gd*y(40)+Yr*y(38)-Zr*y(40)*Cai_la^nc;%---------------TS*_la
ydot(38)=Zr*y(40)*Cai_la^nc-Yr*y(38)+Yp*y(39)-Zp*y(38);%-------TSCa*_la
ydot(39)=Zp*y(38)-Yp*y(39)+fa*y(37)-ga*y(39);%----------------TSCa#_la
ydot(37)=ga*y(39)-fa*y(37)+Yb*TS_la*Cai_la^nc-Zb*y(37);%----------TSCa_la
ydot(41)=Bp*(L_la-y(41)-hpr);%----------Xp_la
ydot(42)=Bw_la*(L_la-y(42)-hwr);%----------Xw_la

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%--------Central venous pressure
Pvc=(SBV-y(28)*Caor-y(30)*Cas-y(34)*Cpa-y(35)*Cpu-0*y(26)-y(27)-y(33))/Cvc;

%-------Left atrial pressure 
Pla=0;
Pla_parall=0;

%--------Left ventricular pressure   new HEART 
Plv=7.5*(Fm_lv*Lm_lv/Lr_lv)*((rout_lv/rint_lv)^2-1)+Gamav_lv*(y(27)-Vo_lv)^3;
Plv_parall=7.5*(Fparall_lv*Lm_lv/Lr_lv)*((rout_lv/rint_lv)^2-1)+Gamav_lv*(y(27)-Vo_lv)^3;

%-------Rigth ventricular pressure 
Prv=7.5*(Fm_rv*Lm_rv/Lr_rv)*((rout_rv/rint_rv)^2-1)+Gamav_rv*(y(33)-Vo_rv)^3;
Prv_parall=7.5*(Fparall_rv*Lm_rv/Lr_rv)*((rout_rv/rint_rv)^2-1)+Gamav_rv*(y(33)-Vo_rv)^3;

Paor=y(28);

load k

if (0<t) && (t<200) && iBeat==1
    %y(35) = y(35) + 9.5; %10
    %Qmt=(y(35)+9 - Plv)/Rmt ; Qmt=Qmt*(Qmt > 0); %7
    %Qmt=(y(35) - Plv)/Rmt - 0.3; Qmt=Qmt*(Qmt > 0); % 7d +0.46
    %Qmt = Qmt+0.2; % 7b
    Qmt = (k-1)*0.032;
    %Qmt = 0.9;
    
else
    Qmt=(y(35) - Plv)/Rmt ; Qmt=Qmt*(Qmt > 0);
end

%Qmt=(y(35) - Plv)/Rmt ; Qmt=Qmt*(Qmt > 0);

dp_lv=Plv-Paor;
Rvaor_eff=Rvaor*1;

if (Inert_lv==0);
Qvaor=dp_lv/Rvaor_eff; Qvaor=Qvaor*(Qvaor > 0);
end
if (Inert_lv~=0);
Qvaor=y(36);
Inert_lv_eff=Inert_lv;
ydot(36)=((((Qvaor>0)+(dp_lv>0))>0)*(dp_lv-Qvaor*Rvaor_eff))/Inert_lv_eff;
Qvaor=Qvaor*(Qvaor>0);
end
ydot(27)=Qmt-Qvaor; %-----------------Left ventricular volume

Ppa=y(34);
Qtc=(Pvc-Prv)/Rtc; Qtc= Qtc*(Qtc > 0);

dp_rv=Prv-Ppa;
Rpv_eff=Rpv*1;

if (Inert_rv==0);
Qpv=dp_rv/Rpv_eff; Qpv=Qpv*(Qpv > 0);
end
if (Inert_rv~=0);
Qpv=y(29);
ydot(29)=((((Qpv>0)+(dp_rv>0))>0)*(dp_rv-Qpv*Rpv_eff))/Inert_rv;
Qpv=Qpv*(Qpv>0);
end
ydot(33)=Qtc-Qpv; % ----------------Rigth ventricular volume

Pla=0.;
ydot(26)=0*(y(35)-Pla)/Rprox-Qmt; % ----Left atrial volume

%-------Pulmonary pressures
Qpul=(y(34)-y(35))/Rpul;
ydot(34)=(Qpv-Qpul)/Cpa;%Pulmonary Pap
ydot(35)=(Qpul-Qmt)/Cpu;%-Pulmonary Ppu
%-------Systemic pressures
if (Inert_ao==0);
Qsys=(y(28)-Pvc)/Rsys;
ydot(30)=0;%------------Pas
% ydot(36)=0; %((y(30)-y(36))/Rsys-(y(36)-Pvc)/Rvs)/Cvs;%-Pvs
ydot(28)=(Qvaor-Qsys)/Caor; %------------------------Paor 
end

Qin_la=(y(35)-Pla)/Rprox;

%% Vm:  Intracellular voltage

ydot(19)=I_app-INa-IK1-Ito-IKr-IKs-ICaL-INCx-INaK-ICap-IKp-ICab-INab;  


%% ----- END EC COUPLING MODEL ---------------
% adjust output depending on the function call

if (nargin == 6)
    output = ydot;
elseif (nargin-1 == 6) && strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 7) && strcmp(runType,'PostProcessing')
       output=[...
       y(26)  ...   % Vla
       Pla  ...   %
       y(27)  ...   % Vlv
       Plv  ...   %
       y(29)  ...   %  Q_Inert_ao
       y(28)*Caor  ...   % Vaor
       y(28)  ...   % Paor
       y(30)*Cas  ...   % Vas
       y(30)  ...   % Pas
       Pvc*Cvc  ...   %  Vvc
       Pvc  ...   %
       y(33)  ...   % Vrv
       Prv  ...   %
       y(34)*Cpa  ...   %  Vpa
       y(34)  ...   % Ppa
       y(35)*Cpu  ...   %  Vpu
       y(35)  ...   % Ppu
       Lm_la  ...   %  Lm_la
       L_la  ...   %  L_la
       Fm_la  ...   %  Fm_la
       Pla_parall   ...   %
       Lm_lv  ...   %  Lm_lv
       L_lv  ...   %  L_lv
       Fm_lv  ...   %  Fm_lv
       Fparall_lv  ...   %  Fparall_lv
       Plv_parall   ...   %
       Lm_rv  ...   %  Lm_rv
       L_rv  ...   %  Lm_rv
       Fm_rv  ...   %  Fm_rv
       Prv_parall   ...   %
       Qin_la  ...   %
       Qmt ...   %
       Qvaor ...   %
       Qsys ...   %
       Qtc ...   %
       Qpv ...   %
       Qpul ...   %
       Irel ...   %
       INa ...   %
       INCx ...   %
       IK1 ...   %
       ICaL ...   %
       INaH ...   %
       IKs ...   %
       Ileak ...   %
       Ilerk ...   %
       Yh ...   %     
        ];
end