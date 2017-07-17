function [P, S] =get_BMrelations(P)
%this function brings in a
S=[];

if isfield(P,'pH')
    % T_fitstr='kb*exp(E/8.314*(1/298.15-1/(273.15+x)))';
    % pH_fitstr='kb.*10^7.*10.^(x-14)+c';
    if isfield(P,'T')
        eT = @(T,E) exp(E*1000/8.314*(1/298.15-1/(273.15+T))); % Arrhenius law
        
        eT = @(T,E) 1+E*(T-25)/(37-25); 
        
        P.kBA= P.kB1*10.^(P.pH-7)*eT(P.T,P.kBT1) +P.kB0*eT(P.T,P.kBT0);
        P.kDA= P.kD1*10.^(P.pH-7)*eT(P.T,P.kDT1) +P.kD0*eT(P.T,P.kDT0);
        P.kEA= P.kE1*10.^(P.pH-7)*eT(P.T,P.kET1) +P.kE0*eT(P.T,P.kDT0);
        P.kFA= P.kF1*10.^(P.pH-7)*eT(P.T,P.kFT1) +P.kF0*eT(P.T,P.kFT0);
        % P.kGA= P.kG1*10^7*10^(P.pH-14) +P.kG0;
        %P.anomer= (P.anomer1*10.^(P.pH-7) +P.anomer0)*eT(P.T,P.anomerT0);
         P.anomer= (P.anomer1*10.^(P.pH-7) +P.anomer0);
        P.R2A= P.R2A1*10.^(P.pH-7)*eT(P.T,P.R2AT1) +P.R2A0*eT(P.T,P.R2AT0);
        if P.T==37
            P.R1A=1/4.4; % 37°
        else
            P.R1A=1/3; % 21°
        end;
        P.Zi= 1 - (1-0.9)*exp(-P.R1A*P.Trec);
    else
        P.kBA= P.kB1*10.^(P.pH-7) +P.kB0;
        P.kDA= P.kD1*10.^(P.pH-7) +P.kD0;
        P.kEA= P.kE1*10.^(P.pH-7) +P.kE0;
        P.kFA= P.kF1*10.^(P.pH-7) +P.kF0;
        % P.kGA= P.kG1*10^7*10^(P.pH-14) +P.kG0;
        P.anomer= P.anomer1*(P.pH-7) +P.anomer0;
        P.R2A= P.R2A1*10.^(P.pH-7) +P.R2A0;
        
    end;
end;



% if "anomer" is given, pools F and B will be explicit dependent
if isfield(P,'anomer') % this brings poool E and F in relation ( makes sense for glucose only)
    P.fE=P.anomer * P.fB;  %P.fE must be a fraction of P.fB (pool with 1 proton), and not of P.Fe itself
    P.fF=(1-P.anomer)*P.fB; % fF is not a free variable anymore, it will be oiverwritten, not fitted and it has to be calculated like this afterwards from fE
    S.fE='P.fE=P.anomer*P.fB; P.fF=(1-P.anomer)*P.fB;';
end;


% and calculate back exchange rate
P.kAB = P.kBA.*P.fB;
P.kAC = P.kCA.*P.fC;
P.kAD = P.kDA.*P.fD;
P.kAE = P.kEA.*P.fE;
P.kAF = P.kFA.*P.fF;
P.kAG = P.kGA.*P.fG;



