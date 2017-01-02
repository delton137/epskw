function y = lineshape_fun2(params,x)
% This is where we compute the sum of the square of the errors.

A1    = abs(params(1));
tau1  = params(2);
w01   = params(3);
A2    = abs(params(4));
tau2  = params(5);

 y =  .5*A1*tau1*x.*( (1 + tau1^2*(x + w01).^2 ).^(-1) + (1 + tau1^2*(x - w01).^2 ).^(-1) )  ...
       + A2*tau2*x.*(1 + x.^2*tau2^2).^(-1) ;
end