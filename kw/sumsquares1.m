function sum_square_errors = sumsquares1(params1,x,y)
% This is where we compute the sum of the square of the errors.

tau1 = params1(2);
A = params1(1);

y_trial = A*exp(-x/tau1); 
diff = y - y_trial;   %errors
sum_square_errors = sum(diff.^2);  % The sum of the squared errors
endfunction 
