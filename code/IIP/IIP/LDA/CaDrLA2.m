function Ca = CaDrLA2(t,period,f_HeartFailure)


Ca_i_max = 1*f_HeartFailure ;    % uM
Ca_i_rest = 0.1 ; % uM
tp = 27 ;          % time to peak

t = mod(t,period) ;
Ca = Ca_i_rest + Ca_i_max*(t/tp).^2.*exp(2*(1-t/tp)) ;

end


