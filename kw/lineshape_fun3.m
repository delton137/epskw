function y = lineshape_fun3(params,x)
% This is where we compute the sum of the square of the errors.

A1    = abs(params(1));
tau1  = params(2);
w01   = params(3);

A2    = abs(params(4));
tau2  = params(5);
w02   = params(6);

D1    = abs(params(7));
Dtau  = params(8);

%  A1*(1/tau1)*w01^2*x.*(   (w01.^2  - x).^2    + (1/tau1).*x  ).^(-1)   ...

 y =   .5*abs(A1)*tau1*x.*( (1 + tau1^2*(x + w01).^2 ).^(-1) + (1 + tau1^2*(x - w01).^2 ).^(-1) )  ...
     + .5*abs(A2)*tau1*x.*( (1 + tau2^2*(x + w02).^2 ).^(-1) + (1 + tau2^2*(x - w02).^2 ).^(-1) ) ; ...
     +  abs(D1)*tau2*x.*(1 + x.^2*Dtau^2).^(-1) ;
end