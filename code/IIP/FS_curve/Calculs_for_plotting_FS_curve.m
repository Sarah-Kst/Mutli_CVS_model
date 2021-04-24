
for l = 1:1:11
    
    dossier = ['Variables',num2str(l),'.mat'];
    load(dossier)
    
    P_LV = Plv;
    P_AO = Pao;
    V_LV = Vlv;
    P_PU = Ppu;
    L_m = Lm_lv;
    F_m = Fm_lv-Fparall_lv;
    
    m = l;
    
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
    
    Preload(m,1) = max(L_m);
    
    Afterload(m,1) = P_AO(n2);
    
    Max_pressure(m,1) = max(P_LV);
    Max_force(m,1) = max(F_m);
    
    Force_mc(m,1) = F_m(n1);
    Force_ao(m,1) = F_m(n2);
    
%     Cal_mc(m,1) = y(n1,16);
%     Cal_ao(m,1) = y(n2,16);
%     
%     Max_Ca(m,1) = max(y(:,16));
    
    S_V(m,1) = max(V_LV)-min(V_LV);
    S_W(m,1) = -trapz(V_LV,P_LV);

 EDV(m,1) = max(V_LV);
 EDP(m,1) = P_LV(n1);


end

T = table(t_init,t_end,t_dur,t_ao,t_ejc,Preload,Afterload,Max_pressure,Max_force,Force_mc,Force_ao,S_V,EDV,S_W,EDP);%,'RowNames',{'Baseline','Increased Ppv'});

T_FS_curve = T;

save T_FS_curve T_FS_curve

clear

load T_FS_curve
