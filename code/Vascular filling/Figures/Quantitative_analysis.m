load T

T.temp = T.Max_pressure - T.Afterload;

T.t_init = (T.t_init-T.t_init(7))/T.t_init(7)*100;
T.t_end = (T.t_end-T.t_end(7))/T.t_end(7)*100;
T.t_dur = (T.t_dur-T.t_dur(7))/T.t_dur(7)*100;
T.t_ao = (T.t_ao-T.t_ao(7))/T.t_ao(7)*100;
T.t_ejc = (T.t_ejc-T.t_ejc(7))/T.t_ejc(7)*100;

T.Preload = (T.Preload-T.Preload(7))/T.Preload(7)*100;
T.Afterload = (T.Afterload-T.Afterload(7))/T.Afterload(7)*100;

T.Max_pressure = (T.Max_pressure-T.Max_pressure(7))/T.Max_pressure(7)*100;
T.Max_force = (T.Max_force-T.Max_force(7))/T.Max_force(7)*100;

T.S_V = (T.S_V-T.S_V(7))/T.S_V(7)*100;

T.Force_mc = (T.Force_mc-T.Force_mc(7))/T.Force_mc(7)*100;
T.Force_ao = (T.Force_ao-T.Force_ao(7))/T.Force_ao(7)*100;

T.Max_Ca = (T.Max_Ca-T.Max_Ca(7))/T.Max_Ca(7)*100;

T.Delta = (T.temp-T.temp(7))/T.temp(7)*100;


T_test = T;

save T_test T_test

clear

load T_test

filename = 'quant_bas.xls';

writetable(T_test,filename)
