Preload = [];

Afterload = [];

for m =1:1:15
    
dossier = ['Variables_NO_FS',num2str(m),'.mat'];
load(dossier)


scolor = jet(15);
dossier = ['Data',num2str(m),'.mat'];
load(dossier)
t = 0:0.08:800;

if size(Plv,2) ~= 10

P_LV = Plv(:,50);
P_AO = Pao(:,50);
V_LV = Vlv(:,50);
P_PU = Ppu(:,50);
L_m = Lm_lv(:,50);
F_m = Fm_lv(:,50);

else
    
P_LV = Plv(:,10);
P_AO = Pao(:,10);
V_LV = Vlv(:,10);
P_PU = Ppu(:,10);
L_m = Lm_lv(:,10);
F_m = Fm_lv(:,10);

end

Nres = find((P_PU<P_LV));
n1 = min(Nres);

Nres = find((P_AO<P_LV));
n2 = min(Nres);

plot(V_LV,P_LV,'g')

hold on

plot(V_LV(n1),P_LV(n1),'ko')

hold on

plot(V_LV(n2),P_LV(n2),'ro')

Preload(m) = V_LV(n1);

Afterload(m) = P_AO(n2);

LM(m) = L_m(n1);

Pmax(m) = max(P_LV);

EDP(m) = P_LV(n1);
EDV(m) = V_LV(n1);

deltL(m) = max(L_m) - min(L_m);

S_V(m) = max(V_LV) - min(V_LV);

a(m) = SBV;

time(m) = t(n1);

duree(m) = t(n2) - t(n1);

Fmax(m) = max(F_m);


end


save SBV_both_2018 Preload Afterload LM Pmax EDP EDV deltL S_V a time duree Fmax

clear
close all

load SBV_both_2018

plot(a,S_V)

hold on

plot(a(4),S_V(4),'bo')

figure

plot(EDV,S_V)





