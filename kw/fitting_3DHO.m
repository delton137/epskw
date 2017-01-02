%% ----Fitting the lineshapes (supplement to kw.m) ------------------------
disp('fitting Lorentzians');
fitcurve = zeros(Nk,num_points); 

A1 = zeros(Nk,1);
w01 = zeros(Nk,1);
tau1 = zeros(Nk,1);
A2 = zeros(Nk,1);
w02 = zeros(Nk,1);
tau2 = zeros(Nk,1);
A3 = zeros(Nk,1);
w03 = zeros(Nk,1);
tau3 = zeros(Nk,1);
D1 = zeros(Nk,1);
Dtau = zeros(Nk,1);

options = optimoptions(@lsqcurvefit,'MaxFunEvals',50000,'TolFun',1e-15,'MaxIter', 10000);

for k = 1
 
    fit_start = 10;         % Start of fit in inverse cm
    fit_end   = 1300;         % End of fit in inverse cm
 
   
    freq_increment = freqs(num_points)*33.4/num_points;
   
    fit_start = ceil( fit_start/freq_increment );
    fit_end   = ceil( fit_end/freq_increment );

    fit_window = fit_start:fit_end;
    freqs_fit = freqs(fit_window);
    chikw_fit = chikw(k,fit_window);
   

    mid = ceil((fit_end - fit_start)/2);
    
    if (k == 1) 
%       params = [.00201,  .35, 920/33.4 ,  ...
%                 .00009,   .9, 785/33.4,   ...
%                 .00007,   .9, 950/33.4,   ...
%                 .003, 2         ]; % Initial guesses for parameters (L)
      params = [0.004,.3   , 400/33.4 ,  ...  
                .00009,   .9, 785/33.4,   ...
                .00007,   .9, 950/33.4,   ...
                  0.1404 ,6.1407        ]; % Initial guesses for parameters (T) 18.4859
    else
       params = [A1(k-1), tau1(k-1), w01(k-1), A2(k-1), tau2(k-1), w03(k-1), A3(k-1), tau3(k-1), w03(k-1), D1(k-1), Dtau(k-1) ]; % use values from previous fit  
    end
    
%      params = lsqcurvefit(@lineshape_fun4, params, freqs_fit',chikw_fit',-1,1,options);
  
  
    A1(k)    = abs(params(1));
    tau1(k)  = params(2);  
    w01(k)   = params(3);

    A2(k)    = abs(params(4));
    tau2(k)  = params(5);
    w02(k)   = params(6);

    A3(k)    = abs(params(7));
    tau3(k)  = params(8);
    w03(k)   = params(9);

    D1(k)    = abs(params(10));
    Dtau(k)  = params(11);
       
    Debye = D1(k)*Dtau(k)*freqs.*(1 + freqs.^2*Dtau(k)^2).^(-1);
    DHO1 =  .5*A1(k)*tau1(k)*freqs.*( (1 + tau1(k)^2*(freqs + w01(k)).^2 ).^(-1) + (1 + tau1(k)^2*(freqs - w01(k)).^2 ).^(-1) ); 
%     DHO1 =  A1(k)*(1/tau1(k))*w01(k)^2*freqs.*(   (w01(k).^2  - freqs.^2).^2    + ((1/tau1(k)).*freqs).^2  ).^(-1);
    DHO2 =  .5*A2(k)*tau2(k)*freqs.*( (1 + tau2(k)^2*(freqs + w02(k)).^2 ).^(-1) + (1 + tau2(k)^2*(freqs - w02(k)).^2 ).^(-1) );
    DHO3 =  .5*A3(k)*tau3(k)*freqs.*( (1 + tau3(k)^2*(freqs + w03(k)).^2 ).^(-1) + (1 + tau3(k)^2*(freqs - w03(k)).^2 ).^(-1) );

    fitcurve(k,:) = Debye + DHO1 + DHO2 + DHO3;
     
%option to plot to see fits   
figure1 = figure(100+k);clf;    
set(figure1,'Position',[10 10 800 700]);
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.164110429447853 0.796454445664105 0.760889570552148],...
    'FontSize',24);
box(axes1,'on');
hold(axes1,'all');
ffcm = 33.44*freqs_fit;
plot1= plot(ffcm, chikw_fit,'k', ffcm,fitcurve(k,fit_window),'g',ffcm, chikw_fit-fitcurve(k,fit_window),'r',ffcm, Debye(fit_window),'b', ffcm, DHO1(fit_window),'o',ffcm, DHO2(fit_window),'y', ffcm,DHO3(fit_window));
set(plot1(1),'DisplayName','data');
set(plot1(2),'DisplayName','fit');
set(plot1(3),'DisplayName','residual');
set(plot1(4),'DisplayName','Debye');
set(plot1(5),'DisplayName','DHO1');
set(plot1(6),'DisplayName','DHO2');
set(plot1(7),'DisplayName','DHO3');


xlabel({'\omega (cm^{-1})'},'FontSize',25);
string = sprintf('\\chi_L(k= %7.1e,\\omega)',k_values(k));
ylabel({string},'FontSize',25);
legend1 = legend(gca,'show');
set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.748701163910984 0.635253391304684 0.142151481888035 0.249233128834356],...
    'FontSize',18);   
     
%   M(n)=getframe;
end


%% Plot w2 vs k  
% figure(4);clf;
% plot(k_values,(1./tau2)*33.44,'r+')
% xlabel('k (ang^-1)');
% legend('w2');
% 

% Plot w vs k  
tau1_cm = tau1/33.44;
w01_cm = w01*33.44;
gamma_cm = gamma*33.44;
w01_max = ((1 + (w01_cm.*tau1_cm).^2).^(1/2))./tau1_cm;
%find propagation speed
fit_lin_start = 1;
fit_lin_end = 5;
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