function output = L_inter(x)

%load Data_Fzero_Length
load Data_length

output = interp1(time,L_sarc,x);