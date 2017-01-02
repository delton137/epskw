%Comparison of two different Lorentz oscillator 

A1 = 00201;
tau1 =  .35;
tau1p =  .7;
w01 = 920/33.4;  


DHO1 =  A1*(1/tau1)*w01^2*freqs.*(   (w01.^2  - freqs.^2).^2    + ((1/tau1).*freqs).^2  ).^(-1);
DHO2 =  .5*A1*tau1p*freqs.*( (1 + tau1p^2*(freqs + w01).^2 ).^(-1) + (1 + tau1p^2*(freqs - w01).^2 ).^(-1) ); 


fit_window = 1:2000;

figure1 = figure(200);clf;
set(figure1,'Position',[10 10 800 700]);
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.164110429447853 0.796454445664105 0.760889570552148],...
    'FontSize',24);
box(axes1,'on');
hold(axes1,'all');
ffcm  = 33.44*freqs;
plot1 = plot(ffcm, DHO1,'k', ffcm,DHO2);
set(plot1(1),'DisplayName','traditional DHO ');
set(plot1(2),'DisplayName','Cosine DHO ');
 
xlabel({'\omega (cm^{-1})'},'FontSize',25);
string = sprintf('\\chi_L(k= %7.1e,\\omega)',k_values(k));
ylabel({string},'FontSize',25);
legend1 = legend(gca,'show');
set(legend1,'YColor',[1 1 1],'XColor',[1 1 1],...
    'Position',[0.748701163910984 0.635253391304684 0.142151481888035 0.249233128834356],...
    'FontSize',18);   
     