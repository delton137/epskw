function sum_square_errors = sumsquares2(params,x,y)
% This is where we compute the sum of the square of the errors.
A = params(1);
tau1 = params(2);
B = params(3);
tau2 = params(4);
%C = params(5);

y_trial = A*exp(-x/tau1) + B*exp(-x/tau2); 
diff = y - y_trial;   %errors
sum_square_errors = sum(diff.^2);  % The sum of the squared errors
endfunction 
