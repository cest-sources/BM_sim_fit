% set initial values of CEST and MT pools (required for numeric solution)

% do the following to initialize, or call 

  Sim=init_Sim(struct());


%%

Sim.n_cest_pool=1;
Sim.dwA=0;
Sim.R1A=1/2;
Sim.R2A=1/1;
% Sim.R1A=1/1.820;  % GM3T
% Sim.R2A=1/0.099;  % GM3T

Sim.fB=0.0001406;
Sim.kBA=7000;
Sim.dwB=55; % ppm  relative to dwA
Sim.R1B=1; 
Sim.R2B=50;

Sim.MT        = 0; 
%MT
% Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
% Sim.MT_lineshape      = 'SuperLorentzian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
% Sim.R1C=1;
% Sim.fC=0.05;
% Sim.R2C=109890;  % 1/9.1µs
% Sim.kCA=40; Sim.kAC=Sim.kCA*Sim.fC;

Sim.FREQ=300;
Sim.Zi=0.9;
Sim.B1=25;		 % the saturation B1 in µT
Sim.tp=5;		 % [s]
Sim.n=1;	Sim.DC=1;
Sim.shape='block';		 
Sim.pulsed=0;		

Sim.xZspec = -80:0.1:80;


Sim.analytic          = 1;                % calculate analtical solution? 1=yes, 0=no
Sim.numeric           = 1;                % calculate numerical solution? 1=yes, 0=no
Sim.MT                = 0;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.Rex_sol           = 'Hyper';          % solution for Rex - cases: 'Hyper', 
Sim.MT_lineshape      = 'Gaussian';       % MT lineshape -SuperLorentzian, Gaussian, Lorentzian
Sim.MT_sol_type       = 'Rex_MT';         % Rex_MT solution type - cases: 'Rex_MT'
Sim.B1cwpe_quad   = -1;                     %XX
Sim.c     = 1;                              %XX
Sim.dummies=0; 
Sim.flipangle=5; 
Sim.readout='FID'; 
Sim.spoilf=0; 



num = NUMERIC_SIM(Sim);
figure(1), plot(num.x,num.zspec,'.'); hold on;
ana = ANALYTIC_SIM(Sim)
figure(1), plot(ana.x,ana.zspec,ana.x,ana.zspec(end:-1:1)-ana.zspec ); hold on;


expname='PARACEST'

if ~exist('expn')
    expn=1;
end;

T = table(num.x,num.zspec,num.zspec(end:-1:1)-num.zspec);
T.Properties.VariableNames = {'offsets' 'Z' 'Asym'}
filename = 'simulation_data.xlsx';
writetable(T,filename,'Sheet',sprintf('exp%d_%s',expn,expname));
expn=expn+1;
