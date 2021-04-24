load T_BL
load T_FS
load T_NO_FS

T_test = [T_BL;T_FS;T_NO_FS];


save T_tot T_test

clear

load T_tot

T_test.t_init = (T_test.t_init/T_test.t_init(1)-1)*100;
T_test.t_end = (T_test.t_end/T_test.t_end(1)-1)*100;
T_test.t_dur = (T_test.t_dur/T_test.t_dur(1)-1)*100;
T_test.t_ao = (T_test.t_ao/T_test.t_ao(1)-1)*100;
T_test.t_ejc = (T_test.t_ejc/T_test.t_ejc(1)-1)*100;

T_test.Preload = (T_test.Preload/T_test.Preload(1)-1)*100;
T_test.Afterload = (T_test.Afterload/T_test.Afterload(1)-1)*100;

T_test.Max_pressure = (T_test.Max_pressure/T_test.Max_pressure(1)-1)*100;
T_test.Max_force = (T_test.Max_force/T_test.Max_force(1)-1)*100;

T_test.S_V = (T_test.S_V/T_test.S_V(1)-1)*100;
T_test.EDV = (T_test.EDV/T_test.EDV(1)-1)*100;

T_test.Force_mc = (T_test.Force_mc/T_test.Force_mc(1)-1)*100;
T_test.Force_ao = (T_test.Force_ao/T_test.Force_ao(1)-1)*100;


save T_tot_rel T_test

clear

load T_tot_rel

filename = 'Tableau_rel_BL.xlsx';
writetable(T_test,filename)

clear

%%

load T_FS
load T_NO_FS

T_test = [T_FS;T_NO_FS];


save T_tot T_test

clear

load T_tot

T_test.t_init = (T_test.t_init/T_test.t_init(1)-1)*100;
T_test.t_end = (T_test.t_end/T_test.t_end(1)-1)*100;
T_test.t_dur = (T_test.t_dur/T_test.t_dur(1)-1)*100;
T_test.t_ao = (T_test.t_ao/T_test.t_ao(1)-1)*100;
T_test.t_ejc = (T_test.t_ejc/T_test.t_ejc(1)-1)*100;

T_test.Preload = (T_test.Preload/T_test.Preload(1)-1)*100;
T_test.Afterload = (T_test.Afterload/T_test.Afterload(1)-1)*100;

T_test.Max_pressure = (T_test.Max_pressure/T_test.Max_pressure(1)-1)*100;
T_test.Max_force = (T_test.Max_force/T_test.Max_force(1)-1)*100;

T_test.S_V = (T_test.S_V/T_test.S_V(1)-1)*100;
T_test.EDV = (T_test.EDV/T_test.EDV(1)-1)*100;

T_test.Force_mc = (T_test.Force_mc/T_test.Force_mc(1)-1)*100;
T_test.Force_ao = (T_test.Force_ao/T_test.Force_ao(1)-1)*100;


save T_tot_rel T_test

clear

load T_tot_rel

filename = 'Tableau_rel_IIP.xlsx';
writetable(T_test,filename)

%%

load T_BL
load T_FS
load T_NO_FS

T_test = [T_BL;T_FS;T_NO_FS];

save T_tot T_test

clear

load T_tot

filename = 'Tableau_abs_BL.xlsx';
writetable(T_test,filename)