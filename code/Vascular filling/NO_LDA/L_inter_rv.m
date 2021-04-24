function output = L_inter_rv(x)

%load Data_Fzero_Length
load Data_length_rv

output = interp1(time,L_sarc,x);