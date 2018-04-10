%%%%%%%%%%% BATCH_simulation %%%%%%%%%%%
% - run it section by section
% - results are NOT saved automatically
% - numeric (BM) solution is calculated in NUMERIC_SIM
%   it supports 1 water (A), 1 MT (C) and up to 5 CEST/NOE pools (B,D,E,F,G)
% - analytic solution is calculated in ANALYTIC_SIM
%   it supports 1 water (A), 1 MT (C) and up to 5 CEST/NOE pools (B,D,E,F,G)
% - always BOTH solutions are calculated and plotted

% last change: 2018/04/10

%% clear all
clear all
clc



%% BASIC SIMULATION PARAMETERS
Sim=init_Sim(struct());
Sim.analytic          = 1;                % calculate analtical solution? 1=yes, 0=no
Sim.numeric           = 0;                % calculate numerical solution? 1=yes, 0=no
Sim.all_offsets       = 1;                % 1 = z-spectrum, 0 = simulate only +- offset, 2 = on-resonant
Sim.offset            = 70;                % offset range in ppm
Sim.n_cest_pool       = 5;                % number of CEST/NOE pools (CEST pools: B,D,E,F,G...)


%% SEQUENCE AND SCANNER PARAMETERS

   
Sim.FREQ          =     7*gamma_;           % frequency (=B0[T] * gamma)
Sim.B1            =     5;                 % standard B1 value in µT
Sim.Trec          =     2.5;                  % standard recover time in s
Sim.Zi            =     0.1;                  % initial magnetisation (should be between -1 and +1)
Sim.shape         =     'block';            % cases: SPINLOCK, seq_gauss, block, AdiaSL, AdiaSinCos, AdiaInversion,
                                            % block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
Sim.pulsed        =     0;                  % 0 = cw saturation, 1 = pulsed saturation

% settings for pulsed saturation
if Sim.pulsed       
    Sim.n     = 7;              % number of saturation pulses
    Sim.tp    = 0.1;            % saturation time per pulse in s
    Sim.DC    = 0.70;             % duty cycle
    
% settings for continuous wave (cw) saturatation
else
    Sim.tp    = 3;           % saturation time in s
    Sim.n     = 1;              % choose n=1 for cw saturation
    Sim.shape = 'block';        % choose 'block' for cw saturation
    Sim.DC    = 1.0;            % choose DC=1 for cw saturation
end;


%% CHOOSE TISSUE-TYPE AND CEST-AGENT AND LOAD PARAMETERS WITH getSim

% Pool system parameters  
% water pool A
Sim.dwA=0;
Sim.R1A=1/3;  % PBS
Sim.R2A=1/2;    % PBS
% Sim.R1A=1/1.820;  % GM3T
% Sim.R2A=1/0.099;  % GM3T
       

% first CEST pool B (paraCEST pool)
Sim.fB=0.00009;  % rel. conc 10mM/111M
Sim.kBA=1000;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwB=55;     % ppm  relative to dwA
Sim.R1B=1;      % R1B relaxation rate [Hz]
Sim.R2B=50;     % R2B relaxation rate [Hz]

% second CEST pool D (amide)
Sim.fD=0.00009;  % rel. conc 10mM/111M
Sim.kDA=6000;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwD=-5;     % ppm  relative to dwA
Sim.R1D=1;      % R1B relaxation rate [Hz]
Sim.R2D=50;     % R2B relaxation rate [Hz]

Sim.n_cest_pool=2;
Sim.MT        = 0; 


%MT
% Sim.MT                = 1;                % 1 = with MT pool (pool C), 0 = no MT pool
% Sim.MT_lineshape      = 'SuperLorentzian';       % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
% Sim.R1C=1;
% Sim.fC=0.05;
% Sim.R2C=109890;  % 1/9.1µs
% Sim.kCA=40; Sim.kAC=Sim.kCA*Sim.fC;

SimStart          = Sim;

%% CHOOSE PARAMETER'S VALUES TO BE SIMULATED
% !!! dont delete parameters - just comment them out !!!
clear Space
% settings for all_offsets = 1 (complete Z-spectrum)
if Sim.all_offsets == 1
    
   Space.tp    = [0.5 1 2 3 4];
   Space.B1    = [5 10 15 25 50 100];
   Space.Trec  = [0 1 2 5 10 20];


else % settings for all_offsets = 0 / 2

   Space.tp    = linspace(0.05,20,50);
   Space.B1    = logspace(-1,2,50);
   Space.Trec  = linspace(0,20,50);

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
        Sim.Zi= 1- (1-SimStart.Zi)*exp(-Sim.R1A*Sim.Trec);
        
        % adjust your spectral resolution here
        if Sim.all_offsets == 1
            Sim.xZspec  = [Sim.dwA-Sim.offset:0.5:Sim.dwA+Sim.offset];
            Sim.xxZspec = [Sim.dwA-Sim.offset:0.5:Sim.dwA+Sim.offset];

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
Sim.modelfield = 'zspec'; %standard value = 'zspec'
PLOT_SPACE(Sim,Space,NUMERIC_SPACE,ANALYTIC_SPACE);

