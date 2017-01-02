%% ----Fitting the lineshapes ---------------------------------------------
disp('fitting One Resonance');

fitcurve = zeros(Nk,num_points); 

A1 = zeros(Nk,1);
w01 = zeros(Nk,1);
tau1 = zeros(Nk,1);
gamma = zeros(Nk,1);
options = optimoptions(@lsqcurvefit,'MaxFunEvals',50000,'TolFun',1e-10,'MaxIter', 1000);


for k = 1
    fit_start = 1300;         % Start of fit in inverse cm
    fit_end   = 2200;         % End of fit in inverse cm
%     fit_start = 2500;         % Start of fit in inverse cm
%     fit_end   = 4000;         % End of fit in inverse cm
    mid_cm = (fit_end+fit_start)/2;
  
    freq_increment = freqs(num_points)*33.4/num_points;
   
    fit_start = ceil( fit_start/freq_increment );
    fit_end   = ceil( fit_end/freq_increment );

    fit_window = fit_start:fit_end;
    freqs_fit = freqs(fit_window);
    chikw_fit = chikw(k,fit_window);

    mid = ceil((fit_end - fit_start)/2);
    
    if (k == 1) 
        params = [.001, .005, mid_cm/33.4/2]; % BE 
    else
        params = [A1(k-1), tau1(k-1), w01(k-1) ]; % use values from previous fit  
    end
       %     fit_start = 2500;         % Start of fit in inverse cm
%     fit_end   = 4000;         % End of fit in inverse cm

%     params  = lsqcurvefit(@lineshape_fun1DHO, params, freqs_fit',chikw_fit',-1,1,options);
    params  = lsqcurvefit(@lineshape_fun1, params, freqs_fit',chikw_fit',-1,1,options);

    A1(k)    = params(1);
    tau1(k)  = params(2);
    w01(k)   = params(3);
    gamma(k) = 1./tau1(k); 
    
  fitcurve(k,:) = .5*A1(k)*tau1(k)*freqs.*( (1 + tau1(k)^2*(freqs + w01(k)).^2 ).^(-1) + (1 + tau1(k)^2*(freqs - w01(k)).^2 ).^(-1) );    

%    option to plot to see fits   
    figure(50)
    clf;
    plot(33.44*freqs_fit, chikw_fit,'b', 33.44*freqs_fit, fitcurve(k,fit_window),'g',33.44*freqs_fit, chikw_fit-fitcurve(k,fit_window),'r')
    string = sprintf('%2i  k = %7.1e, w0 = %7.4f, gamma = %7.4f',k,k_values(k), w01(1)*33.4, gamma(1)*33.4);
    title(string);%  
%   M(n)=getframe;
end
 


%% Plot w vs k  
% tau1_cm = tau1/33.44;
% w01_cm = w01*33.44;
% gamma_cm = gamma*33.44;
% w01_max = ((1 + (w01_cm.*tau1_cm).^2).^(1/2))./tau1_cm;
% %find propagation speed
% fit_lin_start = 1;
% fit_lin_end = 3;
% range_lin_fit = [fit_lin_start:fit_lin_end];
% lin_fit = polyfit(k_values(range_lin_fit),w01_cm(range_lin_fit),1);
% fit_values = lin_fit(1)*k_values(range_lin_fit) + lin_fit(2);
% speed = (lin_fit(2)/33.4)*100;
% 
% w01_cm = abs(w01_cm);
% gamma_cm = abs(gamma_cm);
% w01_max = abs(w01_max);
% 
% figure(3);clf;
% plot(k_values,w01_cm,'r+',k_values,gamma_cm,'g*',k_values,w01_max,k_values(range_lin_fit),fit_values)
% xlabel('k (ang^-1)');
% legend('w0','gamma','w0_central');
% string = sprintf('propataion speed = %7.1f m/s',speed);
% title(string);%     M(n)=getframe;
% xlim([0,8]);

 