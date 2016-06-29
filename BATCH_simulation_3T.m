%%%%%%%%%%% BATCH_simulation %%%%%%%%%%%
% - run it section by section
% - results are NOT saved automatically
% - numeric (BM) solution is calculated in NUMERIC_SIM
%   it supports 1 water (A), 1 MT (C) and up to 5 CEST/NOE pools (B,D,E,F,G)
% - analytic solution is calculated in ANALYTIC_SIM
%   it supports 1 water (A), 1 MT (C) and up to 5 CEST/NOE pools (B,D,E,F,G)
% - always BOTH solutions are calculated and plotted
% - always several analytic solutions are calculated, switch them by
%   changing 'Sim.modelfield'
%   cases:  zspec_santyr (n pools, most simple solution)
%           zspec_CERT (2 pools) 
%           zspec_pulsedCESL (2 pools)
%           zspec_pulsedCESL_approx (DEFAULT!, most accurate n-pool sol.)
% last change: 2014/07/03 by PS

%% clear all
clear all
clc


%% BASIC SIMULATION PARAMETERS
Sim.analytic          = 0;                % calculate analtical solution? 1=yes, 0=no
Sim.numeric           = 1;                % calculate numerical solution? 1=yes, 0=no
Sim.all_offsets       = 1;                % 1 = z-spectrum, 0 = simulate only +- offset, 2 = on-resonant
Sim.offset            = 3;                % offset range in ppm
Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
Sim.n_cest_pool       = 5;                % number of CEST/NOE pools (CEST pools: B,D,E,F,G...)
Sim.Rex_sol           = 'Hyper';          % solution for Rex - cases: 'Hyper', 'Lorentz' , 'minilorentz'
Sim.MT_lineshape      = 'Gaussian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.MT_sol_type       = 'Rex_MT';         % Rex_MT solution type - cases: 'Rex_MT'


%% SEQUENCE AND SCANNER PARAMETERS

Sim.TR  = 3/1000;
Sim.linestomeasure=1;
Sim.flipangle= 14;
Sim.readout='bssfp';
Sim.dummies=1;  %% or better shots?
    
    
Sim.FREQ          =     3*gamma_;           % frequency (=B0[T] * gamma)
Sim.B1            =     1.5;                 % standard B1 value in µT
Sim.Trec          =     2.5;                  % standard recover time in s
Sim.spoilf        =     0;                  % spoilingfactor (0 = full spoiling, 1= no spoiling)
Sim.Zi            =     0.5;                  % initial magnetisation (should be between -1 and +1)
Sim.shape         =     'SPINLOCK';            % cases: SPINLOCK, seq_gauss, block, AdiaSL, AdiaSinCos, AdiaInversion,
                                            % block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
                                            
Sim.pulsed        =     1;                  % 0 = cw saturation, 1 = pulsed saturation

% settings for pulsed saturation
if Sim.pulsed       
    Sim.n     = 7;              % number of saturation pulses
    Sim.tp    = 0.1;            % saturation time per pulse in s
    Sim.DC    = 0.70;             % duty cycle
    
% settings for continuous wave (cw) saturatation
else
    Sim.tp    = 20;           % saturation time in s
    Sim.n     = 1;              % choose n=1 for cw saturation
    Sim.shape = 'block';        % choose 'block' for cw saturation
    Sim.DC    = 1.0;            % choose DC=1 for cw saturation
end;

%XXXXXXXXXXXXXXX DONT CHANGE ! XXXXXXXXXXXXXXXX
Sim.B1cwpe_quad   = -1;                     %XX
Sim.td    = calc_td(Sim.tp,Sim.DC);         %XX
Sim.c     = 1;                              %XX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% CHOOSE TISSUE-TYPE AND CEST-AGENT AND LOAD PARAMETERS WITH getSim
% !!! dont change cases in getSim -  PLEASE create your own !!!

Sim.tissue        = 'GM';
Sim.CESTagent     = 'glucose';
Sim.pH            = 7;                        % pH only relevant for CESTagent 'creatine_pH_T'
Sim.temperature   = 37;                       % temperature only relevant for CESTagent 'creatine_pH_T'

Sim               = getSim(Sim);
SimStart          = Sim;


%% CHOOSE PARAMETER'S VALUES TO BE SIMULATED
% !!! dont delete parameters - just comment them out !!!

% settings for all_offsets = 1 (complete Z-spectrum)
if Sim.all_offsets == 1
    
     Space.n     = [2 3 4 5 6 7 10 15 30];
     Space.tp    = [0.04 0.06 0.08 0.1];
     Space.B1    = [1 1.5 2.2 2.4 3 3.2 3.4 3.6 4 4.2];
     Space.Trec  = [0:0.5:10 20];
     Space.DC  = [ 0.4 0.5 0.7 0.8 0.99];
      
%     Space.dummies =[1 2 3];
%     Space.TR  = [2 2.5 3 4]/1000;
%     Space.flipangle= [5 10 15 20];



else % settings for all_offsets = 0 / 2

     Space.n     = linspace(1,100,10);
%     Space.DC    = linspace(0.01,1,100);
%     Space.tp    = linspace(0.05,0.055,1000);
%      Space.B1    = logspace(-1,2,10);
     Space.B1    = [ 1 2 4 8];
     
     Space.Trec  = [0 0.5 2 10];

end;

names = fieldnames(Space);

for jj = 1:numel(names)
    Sim       = SimStart;   % set back to standard Parameter
    field     = names{jj}; 
    
    for ii = 1:numel(Space.(field))
        
        % do NOT change
        Sim.(field)   = Space.(field)(ii);
        Sim.kAB = Sim.kBA*Sim.fB;
        Sim.kAC = Sim.kCA*Sim.fC;
        Sim.kAD = Sim.kDA*Sim.fD;
        Sim.kAE = Sim.kEA*Sim.fE;
        Sim.kAF = Sim.kFA*Sim.fF;
        Sim.kAG = Sim.kGA*Sim.fG;
        Sim.td  = calc_td(Sim.tp,Sim.DC);
        
        % adjust your spectral resolution here
        if Sim.all_offsets == 1
            Sim.xZspec  = [Sim.dwA-Sim.offset:0.1:Sim.dwA+Sim.offset];
            Sim.xxZspec = [Sim.dwA-Sim.offset:0.1:Sim.dwA+Sim.offset];

        elseif Sim.all_offsets == 0 % !!! CHANGED FROM +- dwB to +- Offset !!!
            Sim.xZspec  = [Sim.offset -Sim.offset];
            Sim.xxZspec = [Sim.offset -Sim.offset];
            
        elseif Sim.all_offsets == 2
            Sim.xZspec  = 0;
            Sim.xxZspec = 0;
        end;
        
        % do NOT change
        if Sim.numeric && Sim.analytic
            NUMERIC_SPACE.(field){ii} = NUMERIC_SIM(Sim);
            ANALYTIC_SPACE.(field){ii} = ANALYTIC_SIM(Sim);
        elseif Sim.numeric && ~Sim.analytic
            NUMERIC_SPACE.(field){ii} = NUMERIC_SIM(Sim);
            ANALYTIC_SPACE = NUMERIC_SPACE;
        elseif ~Sim.numeric && Sim.analytic
            ANALYTIC_SPACE.(field){ii} = ANALYTIC_SIM(Sim);
            NUMERIC_SPACE = ANALYTIC_SPACE;
        end
    end
end

clear ii jj names field


%% PLOT SIMULATED ASYM- AND Z-SPECTRA
% - numerical solutions are displayed as diamonds
% - analytical solutions are displayed as solid lines
% - always several analytic solutions are calculated, switch them by
%   changing "Sim.modelfield"
%   cases:  zspec_santyr (n pools, most simple solution)
%           zspec_CERT (2 pools) 
%           zspec_pulsedCESL (2 pools)
%           zspec_pulsedCESL_approx (n pools, most accurate n-pool sol.)
Sim.modelfield = 'zspec'; %standard value = 'zspec'
PLOT_SPACE(Sim,Space,NUMERIC_SPACE,ANALYTIC_SPACE);

