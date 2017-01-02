function sum_square_errors = sumsquares(params,x,y)
% This is where we compute the sum of the square of the errors.

A1    = params(1);
tau1  = params(2);
w01   = params(3);
A2    = params(4);
tau2  = params(5);


y_trial =  .5*A1*tau1*x.*( (1 + tau1^2*(x + w01).^2 ).^(-1) + (1 + tau1^2*(x - w01).^2 ).^(-1) ) + A2*tau2*x.*(1 + x.^2*tau2^2).^(-1) ;

diff = y - y_trial;   %errors

sum_square_errors = sum(diff.^2);  % The sum of the squared errors

end