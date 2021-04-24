function output = L_inter_lv(x)

load Data_input

output = interp1(time,L_sarc_lv,x);
