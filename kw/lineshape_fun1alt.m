function y = lineshape_fun1(params,x)
% This is where we compute the sum of the square of the errors.

A1    = params(1);
tau1  = params(2);
w01   = params(3);

  y =  .5*A1*tau1*x.*( (1 + tau1^2*(x + w01).^2 ).^(-1) + (1 + tau1^2*(x - w01).^2 ).^(-1) );
%   y = A1*(1/tau1)*w01^2*x.*(   (w01.^2  - x.^2).^2    + ((1/tau1).*x).^2  ).^(-1) ;

 
end