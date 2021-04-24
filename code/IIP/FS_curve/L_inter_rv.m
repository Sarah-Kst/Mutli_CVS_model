function output = L_inter_rv(x)

load Data_input

output = interp1(time,L_sarc_rv,x);