function [Label, Unit]=getSimLabelUnit()
% comments here
% last change: 2014/04/02 by PS

% XXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXX LABELS XXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXX

% scanner and sequence 
Label.n='n';
Label.theta='\theta';
Label.DC='DC';
Label.tp='t_{p}';
Label.B1='B_1';
Label.Trec='T_{rec}';
Label.Zi='Z_i';
Label.FREQ='B_0';
Label.shape='pulse shape';
Unit.shape='';

Label.TR  = 'TR';
Label.linestomeasure='linestomeasure';
Label.flipangle='FA';
Label.readout='redout type';
Label.dummies='dummies';

%adiabatic settings
Label.AdiaMu = 'AdiaMu';
Label.AdiaBw = 'AdiaBw';
Label.AdiaDur = 'AdiaDur';
Label.AdiaTSL = 'AdiaTSL';
Label.AdiaB1 = 'AdiaB1';
Label.AdiaFreqFactor = 'AdiaFreqFactor';

% TTM pulse block
Label.TTM_rep = 'TTM repetition';
Label.TTM_occurance   = 'TTM occurance';
Label.TTM_pulseshape = 'TTM shape';
Label.TTM_bandwidth = 'TTM bandwidth';
Label.TTM_B1 = 'TTM B1';
Label.TTM_type = 'TTM type';

% pool A (water pool)
Label.R1A='R_{1A}';
Label.R2A='R_{2A}';
Label.dwA='\delta_A';

% pool B (1st CEST pool)
Label.dwB='\delta_B';
Label.fB='f_B';
Label.kBA='k_{BA}';
Label.R1B='R_{1B}';
Label.R2B='R_{2B}';

% pool C (MT Pool !!!)
Label.dwC='\delta_C';
Label.fC='f_C';
Label.kCA='k_{CA}';
Label.R1C='R_{1C}';
Label.R2C='R_{2C}';

% pool D (2nd CEST pool)
Label.dwD='\delta_D';
Label.fD='f_D';
Label.kDA='k_{DA}';
Label.R1D='R_{1D}';
Label.R2D='R_{2D}';

% pool E (3rd CEST pool)
Label.dwE='\delta_E';
Label.fE='f_E';
Label.kEA='k_{EA}';
Label.R1E='R_{1E}';
Label.R2E='R_{2E}';

% pool F (4th CEST pool)
Label.dwF='\delta_F';
Label.fF='f_F';
Label.kFA='k_{FA}';
Label.R1F='R_{1F}';
Label.R2F='R_{2F}';

% pool G (5th CEST pool)
Label.dwG='\delta_G';
Label.fG='f_G';
Label.kGA='k_{GA}';
Label.R1G='R_{1G}';
Label.R2G='R_{2G}';

% intramolecular exchange
Label.kBD = 'k_{BD}';
Label.kBE = 'k_{BE}';
Label.kBF = 'k_{BF}';
Label.kBG = 'k_{BG}';

Label.kDE = 'k_{DE}';
Label.kDF = 'k_{DF}';
Label.kDG = 'k_{DG}';

Label.kEF = 'k_{EF}';
Label.kEG = 'k_{EG}';

Label.kFG = 'k_{FG}';

Label.anomer = 'anomer';

% pH exchaneg constants
Label.kB0 = 'k_{B0}';
Label.kB1 = 'k_{B1}';
Label.kD0 = 'k_{D0}';
Label.kD1 = 'k_{D1}';
Label.kE0 = 'k_{E0}';
Label.kE1 = 'k_{E1}';
Label.kF0 = 'k_{F0}';
Label.kF1 = 'k_{F1}';
Label.anomer0 = 'anomer_{0}';
Label.anomer1 = 'anomer_{1}';
Label.R2A0 = 'R_{2A0}';
Label.R2A1 = 'R_{2A1}';

% pH exchaneg constants
Unit.kB0 = '';
Unit.kB1 = '';
Unit.kBT0 = '';
Unit.kBT1 = '';
Unit.kD0 = '';
Unit.kD1 = '';
Unit.kDT0 = '';
Unit.kDT1 = '';
Unit.kE0 = '';
Unit.kE1 = '';
Unit.kET0 = '';
Unit.kET1 = '';
Unit.kF0 = '';
Unit.kF1 = '';
Unit.kFT0 = '';
Unit.kFT1 = '';
Unit.anomer0 = '';
Unit.anomer1 = '';
Unit.anomerT0 = '';
Unit.R2A0 = '';
Unit.R2A1 = '';
Unit.R2AT0 = '';
Unit.R2AT1 = '';

% XXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXX UNITS XXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXX

% scanner and sequence
Unit.DC='';
Unit.n='';
Unit.theta='degree';
Unit.tp='s';
Unit.B1='µT';
Unit.Trec='s';
Unit.Zi='';
Unit.FREQ='MHz';
Unit.AdiaTSL='s';

Unit.TR  = 's';
Unit.linestomeasure='';
Unit.flipangle='°';
Unit.readout='';
Unit.dummies='';

%adiabatic settings
Unit.AdiaMu = '';
Unit.AdiaBw = 'Hz';
Unit.AdiaDur = 's';
Unit.AdiaTSL = 's';
Unit.AdiaB1 = 'µT';
Unit.AdiaFreqFactor = '';

% TTM pulse block
Unit.TTM_rep='';
Unit.TTM_occurance='';
Unit.TTM_pulseshape='';
Unit.TTM_bandwidth='ppm';
Unit.TTM_B1='µT';
Unit.TTM_type='';

% pool A (water pool)
Unit.R1A='Hz';
Unit.R2A='Hz';
Unit.dwA='ppm';

% pool B (1st CEST pool)
Unit.dwB='ppm';
Unit.fB='';
Unit.kBA='Hz';
Unit.R1B='Hz';
Unit.R2B='Hz';

% pool C (MT pool !!!)
Unit.dwC='ppm';
Unit.fC='';
Unit.kCA='Hz';
Unit.R1C='Hz';
Unit.R2C='Hz';

% pool D (2nd CEST pool)
Unit.dwD='ppm';
Unit.fD='';
Unit.kDA='Hz';
Unit.R1D='Hz';
Unit.R2D='Hz';

% pool E (3rd CEST pool)
Unit.dwE='ppm';
Unit.fE='';
Unit.kEA='Hz';
Unit.R1E='Hz';
Unit.R2E='Hz';

% pool F (4th CEST pool)
Unit.dwF='ppm';
Unit.fF='';
Unit.kFA='Hz';
Unit.R1F='Hz';
Unit.R2F='Hz';

% pool G (5th CEST pool)
Unit.dwG='ppm';
Unit.fG='';
Unit.kGA='Hz';
Unit.R1G='Hz';
Unit.R2G='Hz';

% intramolecular exchange
Unit.kBD = 'Hz';
Unit.kBE = 'Hz';
Unit.kBF = 'Hz';
Unit.kBG = 'Hz';

Unit.kDE = 'Hz';
Unit.kDF = 'Hz';
Unit.kDG = 'Hz';

Unit.kEF = 'Hz';
Unit.kEG = 'Hz';

Unit.kFG = 'Hz';

Unit.anomer = '';

