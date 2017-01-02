times = [0:.01:2];
x = sin(times);


n = pow2(nextpow2(length(times)));  % Transform length
y = fft(corr_funs(:,5),n);           % DFT
y2 = ifft(corr_funs(:,5),n);           % DFT

f = (0:n-1)*(1/(n*timestep));      % Frequency range


plot(f,y)
figure(2); 
plot(f,y2) 
    