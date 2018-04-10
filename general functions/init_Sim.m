function Sim=init_Sim(Sim)
%% init_pools(Sim)
% sets all pools to very small f values to avid zeroes in the BM matrix.
Sim.fA=1;
Sim.fB=1e-32;Sim.R1B=1;Sim.R2B=1;Sim.kBA=0;Sim.kAB=0;Sim.dwB=0;
Sim.fC=1e-32;Sim.R1C=1;Sim.R2C=1;Sim.kCA=0;Sim.kAC=0;Sim.dwC=0;
Sim.fD=1e-32;Sim.R1D=1;Sim.R2D=1;Sim.kDA=0;Sim.kAD=0;Sim.dwD=0;
Sim.fE=1e-32;Sim.R1E=1;Sim.R2E=1;Sim.kEA=0;Sim.kAE=0;Sim.dwE=0;
Sim.fF=1e-32;Sim.R1F=1;Sim.R2F=1;Sim.kFA=0;Sim.kAF=0;Sim.dwF=0;
Sim.fG=1e-32;Sim.R1G=1;Sim.R2G=1;Sim.kGA=0;Sim.kAG=0;Sim.dwG=0;


% set initial values for intramolecular exchange
Sim.kBD=0; Sim.kBE=0; Sim.kBF=0; Sim.kBG=0;
Sim.kDE=0; Sim.kDF=0; Sim.kDG=0;
Sim.kEF=0; Sim.kEG=0;
Sim.kFG=0;



% simulation initializations
Sim.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
Sim.MT_sol_type       = 'Rex_MT';         % Rex_MT solution type - cases: 'Rex_MT'
Sim.MT_lineshape= 'Lorentzian';
Sim.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
%XXXXXXXXXXXXXXX DONT CHANGE ! XXXXXXXXXXXXXXXX
Sim.B1cwpe_quad   = -1;                     % this is the B1 power equivalent flag: -1 means that P.B1 is the average B1 over a single pulse 
Sim.c             = 1;                      % this is a free parameter to play around in the solutions  
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% sequence/scanner parameters
% readout
Sim.dummies=0;  %% 
 Sim.TR  = 3/1000;
 Sim.linestomeasure=1;
 Sim.flipangle= 90;
 Sim.readout='FID';
 Sim.TTM_rep = 0;    
 Sim.DC=0.5;