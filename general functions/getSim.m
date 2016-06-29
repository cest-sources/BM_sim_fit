function Sim=getSim(Sim)
% comments here
% last change: 2014/04/02 by PS

%references:
% Henkelman1993  MT agar parameters
% R. M. Henkelman, X. Huang, Q. S. Xiang, G. J. Stanisz, S. D. Swanson, and M. J. Bronskill, Magn Reson Med 29, 759–766 (1993).

%Singh2011 gagcest and cartilage
% A. Singh, M. Haris, K. Cai, V. B. Kassey, F. Kogan, D. Reddy, H. Hariharan, and R. Reddy, Magnetic Resonance in Medicine: Official Journal of the Society of Magnetic Resonance in Medicine / Society of Magnetic Resonance in Medicine (2011).

% Stanisz2005   MT vivo parameters
% G. J. Stanisz, E. E. Odrobina, J. Pun, M. Escaravage, S. J. Graham, M. J. Bronskill, and R. M. Henkelman, Magnetic Resonance in Medicine 54, 507–512 (2005).

% set initial values of CEST and MT pools (required for numeric solution)
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



switch Sim.tissue
    case 'simulation' % schuenke
        T1=0.100;
        T2=0.090;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
 case 'PBS_3T' %schuenke (not correct yet)
        T1=2.92;
        T2=2;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
         
    case 'agar0.2%' %Zaiss: P. 87 notebook #III
        T1=2.99;
        T2=0.1755;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
    case 'agar0.4%' %Zaiss: P. 87 notebook #III
        T1=2.99;
        T2=0.3597;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
    case 'agar1%' %Zaiss: P. 87 notebook #III
        T1=2.93;
        T2=0.1587;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
    case 'WM'  % stanisz 3T
        T1=1.084;
        T2=0.069;
        Sim.dwA=0.0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.139;
        Sim.R2C=100000;
        Sim.kCA=23;
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=0;
        
    case 'GM' % stanisz 3T
        T1=1.820;
        T2=0.099;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.05;
        Sim.R2C=109890;  % 1/9.1µs
        Sim.kCA=40;
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=-2.6;             %asymmetric MT ?
        
    case 'cartilage3T'
        T1=1.2;       % reddy gagCEST
        T2=0.038;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;   % /reddy gagCEST
        
        Sim.R1C=1/T1;   % stanisz  3T cartilage 55°
        Sim.dwC=0;
        Sim.fC=0.182;
        Sim.R2C=1/(8.3*10^-6);  %8.3ms
        Sim.kCA=60;
        Sim.kAC=Sim.kCA*Sim.fC;        % /stanisz
        
    case 'cartilage7T'
        T1=1.5;
        T2=0.032;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;  %% until here reddy paper
        
        Sim.R1C=1/T1;   % stanisz  3T cartilage 55°
        Sim.dwC=0;
        Sim.fC=0.182;
        Sim.R2C=1/(8.3*10^-6);  %8.3ms
        Sim.kCA=60;
        Sim.kAC=Sim.kCA*Sim.fC;        % /stanisz
        
    case 'muscle' % stanisz 3T
        T1=1.412;
        T2=0.050;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.074; % M0B aus stanisz
        Sim.R2C=114942.53;  % 1/8.7µs
        Sim.kCA=66; % R aus stanisz
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=0;             %asymmetric MT ?
        
    case 'muscle_7T' % stanisz 3T extrapoliert
        T1=1.834;
        T2=0.050;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.074; % M0B aus stanisz
        Sim.R2C=114942.53;  % 1/8.7µs
        Sim.kCA=66; % R aus stanisz
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=0;             %asymmetric MT ?
        
    case 'egg_3T'  % Meißner
        T1=2.37;
        T2=1.07;
        Sim.dwA=0;
        Sim.R1A=1/T1;
        Sim.R2A=1/T2;
        
        %MT
        Sim.R1C=1/T1;
        Sim.fC=0.01;
        Sim.R2C=10000;
        Sim.kCA=20;
        Sim.kAC=Sim.kCA*Sim.fC;
        Sim.dwC=0;
        
    otherwise
       warning('NO case was selected for tissue');
        
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all R1B values are assummed to be R1B=R1A
% R2B values are often guessed or used as amine 66.66 Hz from Sun.

switch Sim.CESTagent
    
    case 'none'
        
    case 'gagOH'        %Zaiss: Singh2011
        Sim.R2B=66.66;
        Sim.dwB=1;
        Sim.fB=0.0027;
        Sim.kBA=1000;
        Sim.kAB=Sim.kBA*Sim.fB;
        
    case 'gagNH'        %Zaiss: Singh2011
        Sim.dwB=3.5;
        Sim.R2B=66.66;
        Sim.fB=0.0027;
        Sim.kBA=50;
        Sim.kAB=Sim.kBA*Sim.fB;
        
    case 'APT'    %Zaiss: van zijl und zhou review: 72mM labile protons
        Sim.R2B=66.66;
        Sim.dwB=3.5;
        Sim.fB=0.007;
        Sim.kBA=30;
        Sim.kAB=Sim.kBA*Sim.fB;
        
    case 'APT_review'  %Zaiss: wie im NMRBIOMED articel Zaiss et al.
        Sim.dwB=3.5;
        Sim.R1B=1/T1;
        Sim.R2B=66.66;
        Sim.fB=0.01;
        Sim.kBA=25;
        Sim.kAB=Sim.kBA*Sim.fB;
        
    case 'APEX_review'  %Zaiss: amine creatine like or glutamate? halt intermediate exchange
        Sim.dwB=1.9;
        Sim.R1B=1/T1;
        Sim.R2B=66.66;
        Sim.fB=0.01;
        Sim.kBA=1000;
        Sim.kAB=Sim.kBA*Sim.fB;
        
    case 'creatine_vivo'  %Rerich: T=37, pH=7.09;
        Sim.dwB=1.9;
        Sim.R1B=Sim.R1A;
        Sim.R2B=66.66;
        Sim.fB=0.007;
        Sim.kBA=1500;
        Sim.kAB=Sim.kBA*Sim.fB;
        
case 'creatine_pulsedCESL_pub'  % laborbuch III.111
        Sim.fB=0.002;
        Sim.kBA=35;           % entspricht ca. pH=6.4 T=19°
        Sim.kAB=Sim.kBA*Sim.fB;
        Sim.R1B=Sim.R1A;
        Sim.R2B=90.9;         % from  fit of suns correction of artifact paper
        Sim.dwB=1.9;          % deltaW_B in ppm
        
    case 'creatine_pH_T'     % creatine k(pH,T) from WEX paper Görke et al.
        dH0=55.84*10^3;
        kb=3.009*10^9;
        R=8.314;
        EA=32.27*10^3;
        T0=273.15;
        T25=T0+25;
        
        try
            pH=Sim.pH;
            T=T0+Sim.temperature;
        catch
            error('You have to give pH and temperature in P.pH and P.temperature.');
        end;
        
        Sim.fB=0.002;          % 55.5 mM
        Sim.kBA=kb*10.^(pH-14+ (dH0+EA)./(R.*log(10)).*(1./T25-1./T));
        Sim.kAB=Sim.kBA*Sim.fB;
        
        Sim.R1B=Sim.R1A;
        Sim.R2B=53;           % from own fit.
        Sim.dwB=1.9;
        
            
    case 'glucose'     % schuenke simulationen (100mM)
        % Pool B
        Sim.dwB=2;
        Sim.R2B=66.66;
        Sim.fB=0.0045; % 0.0045/5 bei 100mM und 5 austauschenden Protonen
        Sim.kBA=5000;
        
        % Pool D
        Sim.dwD=2.2;
        Sim.R2D=66.66;
        Sim.fD=0.0009; % 0.0045/5 bei 100mM und 5 austauschenden Protonen
        Sim.kDA=500;
        
        % Pool E
        Sim.dwE=2.8;
        Sim.R2E=66.66;
        Sim.fE=0.0009; % 0.0045/5 bei 100mM und 5 austauschenden Protonen
        Sim.kEA=1000;
        
        % Pool F
        Sim.dwF=0.6;
        Sim.R2F=66.66;
        Sim.fF=0.0009; % 0.0045/5 bei 100mM und 5 austauschenden Protonen
        Sim.kFA=5000;
        
        % Pool G
        Sim.dwG=1.2;
        Sim.R2G=66.66;
        Sim.fG=0.0009; % 0.0045/5 bei 100mM und 5 austauschenden Protonen
        Sim.kGA=5000;
        
            
    case 'amide_amine'  % Meißner
        % Pool B (Amide)
        Sim.dwB=3.5;
        Sim.R1B=1/T1;
        Sim.R2B=66.66;
        Sim.fB=0.00025;
        Sim.kBA=50;
        Sim.kAB=Sim.kBA*Sim.fB;
        
        % Pool D (Amine)
        Sim.dwD=2.0;
        Sim.R1D=1/T1;
        Sim.R2D=66.66;
        Sim.fD=0.002;
        Sim.kDA=1000;
        Sim.kAD=Sim.kDA*Sim.fD;
        
        % Pool E (NOE)
        Sim.dwE=-3.5;
        Sim.R1E=1/T1;
        Sim.R2E=66.66;
        Sim.fE=0.0075;
        Sim.kEA=1000;
        Sim.kAE=Sim.kEA*Sim.fE;
        
    
        
    case 'BarbiAcid'
        % Pool B (Barbituric Acid)
        Sim.dwB=4.9;
        Sim.R1B=1/T1;
        Sim.R2B=66.66;
        Sim.fB=0.002;
        Sim.kBA=250;
        Sim.kAB=Sim.kBA*Sim.fB;
        
        % Pool D Creatine
        Sim.dwD=1.88;
        Sim.R1D=1/T1;
        Sim.R2D=66.66;
        Sim.fD=0.0015;
        Sim.kDA=100;
        Sim.kAD=Sim.kDA*Sim.fD;
        
               
       otherwise
       warning('NO case was selected for cest pool');
end;

if Sim.n_cest_pool == 0
    Sim.kBA = 0; Sim.kDA = 0; Sim.kEA = 0; Sim.kFA = 0; Sim.kGA = 0;
end

