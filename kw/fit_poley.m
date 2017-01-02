% Poley absorption fit
 disp('fitting Poley');
AP = zeros(Nk,1);
tauP = zeros(Nk,1);
AD = zeros(Nk,1);
tauD = zeros(Nk,1);
for k = 100
    fit_start =.001;         % Start of fit in inverse cm
    fit_end  = 1000;         % End of fit in inverse cm
   
    freq_increment = freqs(num_points)*33.4/num_points;
   
    fit_start = ceil( fit_start/freq_increment );
    fit_end   = floor( fit_end/freq_increment );

    fit_window = fit_start:fit_end;
    freqs_fit = freqs(fit_window);
    chikw_fit = chikw(k,fit_window);

    mid = ceil((fit_end - fit_start)/2);
    
    if (k == 100) 
      paramsP = [.01,  2, .2, 1]; % Initial guesses for parameters (L) TTM3F
    else
      paramsP = [AP(k-1), tauP(k-1), AD(k-1), tauD(k-1) ]; % use values from previous fit  
    end
   
   options = optimoptions(@lsqcurvefit,'MaxFunEvals',100000,'TolFun',1e-13);
    
   paramsP = lsqcurvefit(@lineshape_fun_poley, paramsP, freqs_fit',chikw_fit',-1,1,options);

   AP(k)   = paramsP(1);
   tauP(k) = paramsP(2);
   AD(k)   = paramsP(3);
   tauD(k) = paramsP(4);
   
   part1 = AP(k)*(3.141^(.5))*tauP(k)^2*freqs.*exp(-.25*tauP(k)^2*freqs.^2);
   part2 = 0;%AD(k)*tauD(k)*freqs.*(1 + freqs.^2*tauD(k)^2).^(-1);
   fitcurve(k,:) = part1  + part2;
    
%  option to plot to see fits   
    figure(50)
    clf;
    plot(33.44*freqs_fit, chikw_fit,  33.44*freqs_fit, fitcurve(k,fit_window),'g',33.44*freqs_fit,part1(fit_window),'b',33.44*freqs_fit,part2(fit_window),'y')
    string = sprintf('%2i  k = %7.1e',k,k_values(k));
    title(string);%  
end