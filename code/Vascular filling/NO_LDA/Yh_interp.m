function output = Yh_interp(x)

%load Data_inter
%load Yh_inter
load Yh_data

output = interp1(time,Yhh,x);