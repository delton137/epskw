function sum_square_errors = sumsquares_stretch(params_str,x,y)
% This is where we compute the sum of the square of the errors.
A = params_str(1);
tau1 = params_str(2);
beta1 = params_str(3);

y_trial = A*exp(-(x/tau1).^beta1); 
diff = y - y_trial;   %errors
sum_square_errors = sum(diff.^2);  % The sum of the squared errors
endfunction 
