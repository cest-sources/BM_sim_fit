%% this file allows to create fit Z-spectrum data with BM simulation
%%it also allows to fit single offset QUESP and QUEST formulaes

%% 1.1A   load XLS_qCEST file 
% this creates a stack of Z-spectra Z_x and teh parameters struct P
 [w_x, Z_x, P, rowname, Ztab]=LOAD_xls_qCEST_2_Ztab({'15degC'});
 expname=Ztab{rowname,'exp'}{1};
P

%% 1.1B1 load QUEST and QUESP data from CESTtab structure
%first load a Ztab from folder
% Ztab=ZtabEu;
% or load it from multiple excel sheets
[w_x, Z_x, P, rowname, Ztab]=LOAD_xls_qCEST_2_Ztab();

%% 1.1B2 load dataset from Ztab by rowname
rowname='Dota37degC'  
%rowname= Ztab.Properties.RowNames{1};

expname=Ztab{rowname,'exp'}{1};  P=Ztab{rowname,'P'}{1};  
Ztab(rowname,:) = b0correct_run(Ztab(rowname,:),'B1_run') ;
[ Ztab(rowname,:), Mnorm] = norm_run(Ztab(rowname,:),'B1_run',P.normalized) ;
% Ztab(rowname,:) = exclude_run(Ztab(rowname,:),'B1_run',-80) 
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));

figure(2001), [w_x, Z_x, ~, ~,  varyval, P.vary, P]= plot_tab(Ztab,rowname,'B1_run');
P.varyval=varyval; clear varyval; % this overwrites varyval if only some B1 values are picked.

%% 1.2: create parameter struct "Sim" for simulation/fit model

% simulation parameters
Sim=init_Sim(struct()); % initializes all pools with zeroes and sets some standard values
Sim.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
Sim.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
Sim.MT_lineshape  = 'SuperLorentzian';    % ssMT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
Sim.n_cest_pool   = 2;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)

% MR and sequence parameters
Sim.FREQ          = 7*gamma_;           % frequency (=B0[T] * gamma)
Sim.B1            = 10;                   % B1 value in µT
Sim.Trec          = 3;                    % recover time in s

Sim.Zi            = 0    ;        % initial magnetisation (should be between -1 and +1)
Sim.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation

if Sim.pulsed
    Sim.shape = 'seq_gauss';              % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
    Sim.n     = 20;                      % number of saturation pulses
    Sim.tp    = 0.1;                      % saturation time per pulse in s
    Sim.DC    = 0.5;                   % duty cycle
else
    Sim.tp    = 10;                       % saturation time in s
    Sim.n     = 1;                        % choose n=1 for cw saturation
    Sim.shape = 'block';                  % choose 'block' for cw saturation
    Sim.DC    = 1.0;                      % choose DC=1 for cw saturation
end;

% Pool system parameters  
% water pool A
Sim.dwA=0;
Sim.R1A=0.5;
Sim.R2A=2;

% first CEST pool B (paraCEST pool)
Sim.n_cest_pool=1;
Sim.fB=0.00009;  % rel. conc 10mM/111M
Sim.kBA=1000;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwB=55;     % ppm  relative to dwA
Sim.R1B=1;      % R1B relaxation rate [Hz]
Sim.R2B=50;     % R2B relaxation rate [Hz]

% second CEST pool D (amide)
Sim.n_cest_pool=2;
Sim.fD=0.00009;  % rel. conc 10mM/111M
Sim.kDA=7000;   % exchange rate in Hz ( the fast one, kBA is calculated by this and fB)
Sim.dwD=-5;     % ppm  relative to dwA
Sim.R1D=1;      % R1B relaxation rate [Hz]
Sim.R2D=50;     % R2B relaxation rate [Hz]

%% 1.3 set and optimize start values and boundaries
Sim.n_cest_pool=2;
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL A XXXXXXXXXX
T.varyA       = [   1           0                1             ];
T.dep_varsA   = {   'dwA',      'R1A',          'R2A'          };                  
T.startA      = [   Sim.dwA       Sim.R1A       Sim.R2A        ];                  
T.lowerA      = [   Sim.dwA-1     1/4*Sim.R1A   1/10*Sim.R2A    ];
T.upperA      = [   Sim.dwA+1     2*Sim.R1A     20000*Sim.R2A  ];

[T.dep_varsA, T.startA, T.lowerA, T.upperA] = selectVars( T.varyA, T.dep_varsA, T.startA, T.lowerA, T.upperA );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL B XXXXXXXXXX
T.varyB       = [   1               1           1             1      ];
T.dep_varsB   = {   'dwB',        'fB',       'kBA',        'R2B'    };     
T.startB      = [   Sim.dwB       Sim.fB      Sim.kBA       Sim.R2B  ];
T.lowerB      = [   Sim.dwB-10    Sim.fB*0.1  Sim.kBA/1000     0      ];
T.upperB      = [   Sim.dwB+10    Sim.fB*10    Sim.kBA*1000    50     ];

[T.dep_varsB, T.startB, T.lowerB, T.upperB] = selectVars( T.varyB, T.dep_varsB, T.startB, T.lowerB, T.upperB );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL C XXXXXXXXXX  always reserved for MT in numeric BM, thus commented here
% T.varyC       = [ 0         0           0           0           0           ];
% T.dep_varsC   = {'dwC',     'fC',       'kCA',      'R1C',      'R2C'       };     
% T.startC      = [Sim.dwC    Sim.fC      Sim.kCA    Sim.R1C     Sim.R2C       ];
% T.lowerC      = [0          0.0001      10          1/20        1/20        ];
% T.upperC      = [9.9        0.1         10000       1/0.1       1/0.1       ];
% [T.dep_varsC, T.startC, T.lowerC, T.upperC] = selectVars( T.varyC, T.dep_varsC, T.startC, T.lowerC, T.upperC );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL D XXXXXXXXXX
T.varyD       = [ 1            1           1                0            ];
T.dep_varsD   = {'dwD',       'fD',       'kDA',          'R2D',         };     
T.startD      = [Sim.dwD      Sim.fD        Sim.kDA       Sim.R2D        ];
T.lowerD      = [Sim.dwD-10   Sim.fD*0.01    Sim.kDA/1000     0            ];
T.upperD      = [Sim.dwD+3   Sim.fD*100     Sim.kDA*1000    50            ];
 
[T.dep_varsD, T.startD, T.lowerD, T.upperD] = selectVars( T.varyD, T.dep_varsD, T.startD, T.lowerD, T.upperD );

FIT.T=T;   clear T;  % store start values and boundaries in FIT struct
FIT.Sim=Sim; % store Sim values in FIT struct

figure(2002), multiZplot(P,Sim,FIT.T,w_x,Z_x); % plot guess  with data

% Zknown=0.95; % this can be added to guess the initial magnetization Zi
% P.Zi= 1 - (1-Zknown)*exp(-P.R1A*P.Trec);

%% 1.4  multi-Z-BMfitting of the Z_x data (wherever it comes from)
% RUN full BM OPTIMIZATION 
% you need a startvalue, run 1.3 first!
Sim.analytic=1;  % set1 this to 1 if analytic fit should be used, numeric =0 can take forever
Sim.n_cest_pool=2;

% fit-options
[FIT] =multiZfit(P,Sim,FIT.T,w_x,Z_x);

figure,
FIT.Z_fit=multiZplot(P,Sim,FIT.T,w_x,Z_x,FIT.popt,FIT.pci);
title([expname ' : ' rowname]);
% savefig(['BM_FIT_3p' rowname '.fig']);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

FIT.P=P;
FIT.w_x=w_x;
FIT.Z_x=Z_x;

if exist('Ztab') % save fitresult in Ztable
Ztab(rowname,'FIT3p')={{FIT}}; % name the fit for saving in Ztab
end;

%% 1.5 single offset preparation ( QUESP and OMEGA plot)

single_offset=FIT.dwB(1);
figure(), subplot(2,1,1);
[Zlab Zref]=prep_Zlab_Zref(P,single_offset,w_x,Z_x);

% boundaries for   fb        kb
opts.Lower      = [0.0000135 0     ];
opts.StartPoint = [0.000135  4000  ];
opts.Upper      = [0.0135    1500000];

%some evaluations need a value for R2B, the 
P.R2B=0; P.B1cwpe_quad=-1;  P.dwB=single_offset;

[que, que_inv] =QUESP(Zlab,Zref,single_offset,P,opts);

subplot(2,3,6);
ome=OmegaPlot(Zlab,Zref,single_offset,P,opts);

Ztab(rowname,'QUESP')={{que}};
Ztab(rowname,'QUESP_inv')={{que_inv}};
Ztab(rowname,'Omega')={{ome}};

clear que ome;

%% 1.6 save Ztab

uisave('Ztab','Ztab_DOTA')

%% 1.7 reaload Ztab  (after this run 1.4 or 1.5 again)
% open Ztab and check for rowname and also column name of fit(FIT3p here)

if ~exist('Ztab') uiload; end;
rowname='Dota37degC';
Sim=Ztab.FIT3p{rowname}.Sim
FIT=Ztab.FIT3p{rowname}
expname=Ztab{rowname,'exp'}{1};  P=Ztab{rowname,'P'}{1};  
figure(2001), [w_x, Z_x, ~, ~,  varyval, P.vary, P]= plot_tab(Ztab,rowname,'B1_run');

figure(2002), multiZplot(P,Sim,FIT.T,w_x,Z_x,FIT.popt,FIT.pci);

%% PLOT field over Ztabfits
clear kBA fitkBA quespkBA omegakBA
% give rows you want to plot
row_ind=[1:3];   
% give x-axis for this plot, e.g. Temperature or 
% row_x= row_ind; % at elast a number
row_x=[15  25  37 ];


ind=1:numel(row_ind);
for ii=ind
   FITres= Ztab{row_ind(ii),'FIT3p'}{1};
   quesp= Ztab{row_ind(ii),'QUESP_inv'}{1};
   omega= Ztab{row_ind(ii),'Omega'}{1};
   fitkBA(ii,:)=FITres.kBA;
%    fitkDA(ii)=FITres.kDA;
   quespkBA(ii,:)=quesp.kBA;
   omegakBA(ii,:)=omega.kBA;
%    quesp_origkBA(ii)=quesp_orig.kBA;
   
% kBA(ii)=FITres.omega_Rex.kBA;
% kBA(ii)=FITres.omega_Dixon.kBA;
% kBA(ii)=FITres.omega_Rex.fB;
      
end;

figure();
errorbar(row_x,fitkBA(:,1),fitkBA(:,2),'DisplayName','BM fit, kBA'); hold on;
errorbar(row_x,quespkBA(:,1),quespkBA(:,2),'DisplayName','QUESP, kBA'); hold on;
errorbar(row_x,omegakBA(:,1),omegakBA(:,2),'DisplayName','\Omega-plot, kBA'); hold on;

ylabel('exchange rate'); 
xlabel('Temperature [°C]');
legend('Location','northwest');

