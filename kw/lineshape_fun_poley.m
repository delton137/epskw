function y = lineshape_fun_poley(params,x)
% This is where we compute the sum of the square of the errors.

AP   = params(1);
tauP = params(2);
AD   = params(3);
tauD = params(4);

y =  AP*(3.141^(.5))*tauP^2*x.*exp(-.25*tauP^2*x.^2);
%    + AD*tauD*x.*(1 + x.^2*tauD^2).^(-1);

end