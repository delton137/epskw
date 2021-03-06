function y = lineshape_fun(params,x)
% This is where we compute the sum of the square of the errors.

A1    = params(1);
tau1  = params(2);
w01   = params(3);
A2    = params(4);
tau2  = params(5);
A3    = params(6);
tau3  = params(7);
w03   =  params(8);

% y =  .5*A1*tau1*x.*( (1 + tau1^2*(x + w01).^2 ).^(-1) + (1 + tau1^2*(x - w01).^2 ).^(-1) ) + A2*tau2*x.*(1 + x.^2*tau2^2).^(-1) + A3*tau3*x.*(1 + x.^2*tau3^2).^(-1) ;
 y =  .5*A1*tau1*x.*( (1 + tau1^2*(x + w01).^2 ).^(-1) + (1 + tau1^2*(x - w01).^2 ).^(-1) )  ...
       + A2*tau2*x.*(1 + x.^2*tau2^2).^(-1)  ...
       + .5*A3*tau3*x.*( (1 + tau3^2*(x + w03).^2 ).^(-1) + (1 + tau3^2*(x - w03).^2 ).^(-1) ) ;

end