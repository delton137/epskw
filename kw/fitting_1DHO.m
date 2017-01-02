%% ----Fitting the lineshapes ---------------------------------------------
disp('fitting Lorentzians');
 options = optimoptions(@lsqcurvefit,'MaxFunEvals',500000,'TolFun',1e-15,'MaxIter', 5000);


fitcurve = zeros(Nk,num_points); 

A1 = zeros(Nk,1);
w01 = zeros(Nk,1);
tau1 = zeros(Nk,1);
gamma = zeros(Nk,1);
A2 = zeros(Nk,1);
tau2 = zeros(Nk,1);
A3 = zeros(Nk,1);
tau3 = zeros(Nk,1);
w03 = zeros(Nk,1);
nframes = 1;

% writerObj = VideoWriter('fitting.avi');
% writerObj.FrameRate = 1;
% open(writerObj);

for k = 1
%    if (k_values(k) < 3) 
    fit_start = 200;         % Start of fit in inverse cm
    fit_end   = 1300;         % End of fit in inverse cm
%    end
%   if (3 <= k_values(k) < 7)
%     fit_start = 120;        
%     fit_end  = 800;     
%   end
%   if (7 <= k_values(k))
%     fit_start = 50;        
%     fit_end  = 800;     
%   end
   
    freq_increment = freqs(num_points)*33.4/num_points;% writerObj = VideoWriter('fitting.avi');
% writerObj.FrameRate = 1;
% open(writerObj);

   
    fit_start = ceil( fit_start/freq_increment );
    fit_end   = ceil( fit_end/freq_increment );

    fit_window = fit_start:fit_end;
    freqs_fit = freqs(fit_window);
    chikw_fit = chikw(k,fit_window);
   

    mid = ceil((fit_end - fit_start)/2);
    
    if (k == 1) 
%         params = [.004, .2, 410/33.4 ,  .14, 10,  0,0]; % Initial guesses for parameters (T) 300 K 
%           params = [.004, .2, 500/33.4 ,  .14, 10,  0, 0]; % Initial guesses for parameters (T) 300 K 
%         params = [.00005, .2, 550/33.4 ,  .0002, 10, .001, .00001]; % Initial guesses for parameters (T) 250 K 
%       params = [.5, .1, 750/33.4 , .2,  2,   .5, 1/(200/33.4),(200/33.4)]; % Initial guesses for parameters (L)
%       params = [.0012, .35,    730/33.4 ,     .0016,  1,    .000001, .5,180/33.4]; % Initial guesses for parameters (L)
      params = [.0014, .35,    630/33.4 ,     .0016,  1,    0 , 0 ]; % Initial guesses for parameters (L) TTM3F
% 
%     elseif (k == 6)
%         params = [.0013, .24, 680/33.4 ,     .0016,  1,    .00015, .5,180/33.4]; % L 350
%     elseif (k == 8)
%         params = [.0012, .19, 640/33.4 ,     .0015,  1,    .0003, .5,180/33.4]; % L 350
% %     elseif (k == 13)
%         params = [.06, .2, 480/33.4 ,         .8,  2,       .006, .5,200/33.4]; % L 350
    else
       params = [A1(k-1), tau1(k-1), w01(k-1), A2(k-1), tau2(k-1),A3(k-1),tau3(k-1) ]; % use values from previous fit  
    end
    
    
    
% 
%     if (1 <= k_values(k)) && (k_values(k) <= 2.1) 
%         
%         params = lsqcurvefit(@lineshape_fun3, params, freqs_fit',chikw_fit');
%     elseif (3 <= k_values(k)) && (k_values(k) <= 4.5) 
%         
%         params = lsqcurvefit(@lineshape_fun3, params, freqs_fit',chikw_fit',-1,1,options);
%     else
        params(6:8) = 0 ;
        params(1:5) = lsqcurvefit(@lineshape_fun2, params(1:5), freqs_fit',chikw_fit',-1,1,options);
        params(6:8) = 0 ;
% %    end
  
    A1(k)    = abs(params(1));
    tau1(k)  = params(2);
    w01(k)   = params(3);
    gamma(k) = 1./tau1(k); 
    A2(k)    = abs(params(4));
    tau2(k)  = params(5);
    A3(k)    = abs(params(6)); 
    tau3(k)  = params(7); 
    w03(k)   =  190/33.4;
   
    part1 = .5*A1(k)*tau1(k)*freqs.*( (1 + tau1(k)^2*(freqs + w01(k)).^2 ).^(-1) + (1 + tau1(k)^2*(freqs - w01(k)).^2 ).^(-1) );
    part2 = A2(k)*tau2(k)*freqs.*(1 + freqs.^2*tau2(k)^2).^(-1); 
%      part3 = A3(k)*tau3(k)*freqs.*(1 + freqs.^2*tau3(k)^2).^(-1); 
    part3 = .5*A3(k)*tau3(k)*freqs.*( (1 + tau3(k)^2*(freqs + w03(k)).^2 ).^(-1) + (1 + tau3(k)^2*(freqs - w03(k)).^2 ).^(-1) ) ;
    % A3(k)*tau3(k)*freqs.*(1 + freqs.^2*tau3(k)^2).^(-1);
    fitcurve(k,:) = part1 + part2 + part3;
%      
% %    option to plot to see fits   
figure1 = figure(100);
set(figure1,'Position',[10 10 800 700]);
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.164110429447853 0.796454445664105 0.760889570552148],...
    'FontSize',24);
box(axes1,'on');
hold(axes1,'all');
ffcm = 33.44*freqs_fit;
 plot1= plot(33.44*freqs_fit, chikw_fit,'k',33.44*freqs_fit, fitcurve(k,fit_window),'g',33.44*freqs_fit, chikw_fit-fitcurve(k,fit_window),'r',33.44*freqs_fit, part1(fit_window),'b', 33.44*freqs_fit, part2(fit_window)+  part3(fit_window),'y');
 
set(plot1(1),'DisplayName','data');
set(plot1(2),'DisplayName','fit');
set(plot1(3),'DisplayName','residual');
set(plot1(4),'DisplayName','DHO1');
set(plot1(5),'DisplayName','Debye');
 
xlabel({'\omega (cm^{-1})'},'FontSize',25);
string = sprintf('\\chi_T(k= %7.1e,\\omega)',k_values(k));
ylabel({string},'FontSize',25);
legend1 = legend(gca,'show');
set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.748701163910984 0.635253391304684 0.142151481888035 0.249233128834356],...
    'FontSize',18); 
 
% % frame = getframe;
% % writeVideo(writerObj,frame);

end


 
% close(writerObj);


%% Plot w2 vs k  
% figure(4);clf;
% plot(k_values,(1./tau2)*33.44,'r+')
% xlabel('k (ang^-1)');
% legend('w2');
% 
% %Plot w vs k  
tau1_cm = tau1/33.44;
w01_cm = w01*33.44;
gamma_cm = gamma*33.44;
w01_max = ((1 + (w01_cm.*tau1_cm).^2).^(1/2))./tau1_cm;
% find propagation speed
fit_lin_start = 1  ;
fit_lin_end   = 5 ;
range_lin_fit = [fit_lin_start:fit_lin_end];
lin_fit = polyfit(k_values(range_lin_fit),w01_cm(range_lin_fit),1);
fit_values = lin_fit(1)*k_values(range_lin_fit) + lin_fit(2);
speed = (lin_fit(2)/33.4)*100;

w01_cm = abs(w01_cm);
gamma_cm = abs(gamma_cm);
w01_max = abs(w01_max);

figure(3);clf;
plot(k_values,w01_cm,'r+',k_values,gamma_cm,'g*',k_values,w01_max,k_values(range_lin_fit),fit_values)
xlabel('k (ang^-1)');
legend('w0','gamma','w0_central');
string = sprintf('propataion speed = %7.1f m/s',speed);
title(string);%     M(n)=getframe;
xlim([0,8]);


%% %------ Print data ------------------------------------------- 
matrix2save = [k_values,w01_cm,gamma_cm,w01_max];

save -ascii TIP4P2005f_300_T_w_vs_k.dat matrix2save 	
% save -ascii TTM3F_350_L_w_vs_k.dat matrix2save 	



