for m =1:1:15
    
dossier = ['Variables',num2str(m),'.mat'];
load(dossier)

scolor = jet(15);
dossier = ['Data',num2str(m),'.mat'];
load(dossier)
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

t_init(m,1) = t(n1)-min(t);

Nres = find((P_AO<P_LV));
n2 = min(Nres);
n3 = max(Nres);

t_end(m,1) = t(n2)-min(t);

t_dur(m,1) = t(n2) - t(n1);

t_ao(m,1) = t(n3)-min(t);

t_ejc(m,1) = t(n3) - t(n2);

Preload(m,1) = L_m(n1);

Afterload(m,1) = P_AO(n2);

Max_pressure(m,1) = max(P_LV);
Max_force(m,1) = max(F_m);

Force_mc(m,1) = F_m(n1);
Force_ao(m,1) = F_m(n2);

Cal_mc(m,1) = y(n1,16);
Cal_ao(m,1) = y(n2,16);

Max_Ca(m,1) = max(y(:,16));

S_V(m,1) = max(V_LV)-min(V_LV);

end

T = table(t_init,t_end,t_dur,t_ao,t_ejc,Preload,Afterload,Max_pressure,Max_force,Force_mc,Force_ao,Cal_mc,Cal_ao,Max_Ca,S_V);%,'RowNames',{'Baseline','Increased Ppv'});

save SBV_2018_bas T

clear

load SBV_2018_bas








