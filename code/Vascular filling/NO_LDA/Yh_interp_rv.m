function output = Yh_interp_rv(x)

%load Data_inter
%load Yh_inter
load Yh_data_rv

output = interp1(time,Yh_rv,x);