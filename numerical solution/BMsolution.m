function [Z, A, Ainvb, dia] = BMsolution(x0, xZspec, tspan, z0, y0, y1)
% comments here
% last change: 2014/04/02 by PS
% BM solution solves the Bloch-McConnell equations stepwise










% set kCA=0 and fC=0, if MT is disabled
if x0(1,3) == 0;            % (x0(1,3) = Sim.MT)
   x0(4,2) = 0; y0(15) = 0;
end

% set exchange rates of unused pools to zero
switch x0(1,4)                          % (x0(1,4) = Sim.n_cest_pool)
   
    case 0
        x0(3,2) = 0; y0(14) = 1e-16;        % (x0(3,2) = Sim.kBA, y0(14) = Sim.fB)
        x0(5,2) = 0; y0(16) = 1e-16;        % (x0(5,2) = Sim.kDA, y0(16) = Sim.fD)
        x0(6,2) = 0; y0(17) = 1e-16;        % (x0(6,2) = Sim.kEA, y0(17) = Sim.fE)
        x0(7,2) = 0; y0(18) = 1e-16;        % (x0(7,2) = Sim.kFA, y0(18) = Sim.fF)
        x0(8,2) = 0; y0(19) = 1e-16;        % (x0(8,2) = Sim.kGA, y0(19) = Sim.fG)

        % set intramolecular exchange rates of unused pools to zero
        x0(3,8) = 0; x0(5,8) = 0; x0(6,8) = 0; x0(7,8) = 0;
        x0(3,7) = 0; x0(5,7) = 0; x0(6,7) = 0;
        x0(3,6) = 0; x0(5,6) = 0;
        x0(3,5) = 0;
        
    case 1
        x0(5,2) = 0; y0(16) = 1e-16;        % (x0(5,2) = Sim.kDA, y0(16) = Sim.fD)
        x0(6,2) = 0; y0(17) = 1e-16;        % (x0(6,2) = Sim.kEA, y0(17) = Sim.fE)
        x0(7,2) = 0; y0(18) = 1e-16;        % (x0(7,2) = Sim.kFA, y0(18) = Sim.fF)
        x0(8,2) = 0; y0(19) = 1e-16;        % (x0(8,2) = Sim.kGA, y0(19) = Sim.fG)

        % set intramolecular exchange rates of unused pools to zero
        x0(3,8) = 0; x0(5,8) = 0; x0(6,8) = 0; x0(7,8) = 0;
        x0(3,7) = 0; x0(5,7) = 0; x0(6,7) = 0;
        x0(3,6) = 0; x0(5,6) = 0;
        x0(3,5) = 0;
        
    case 2
        x0(6,2) = 0; y0(17) = 1e-16;        % (x0(6,2) = Sim.kEA, y0(17) = Sim.fE)
        x0(7,2) = 0; y0(18) = 1e-16;        % (x0(7,2) = Sim.kFA, y0(18) = Sim.fF)
        x0(8,2) = 0; y0(19) = 1e-16;        % (x0(8,2) = Sim.kGA, y0(19) = Sim.fG)
        
        % set intramolecular exchange rates of unused pools to zero
        x0(3,8) = 0; x0(5,8) = 0; x0(6,8) = 0; x0(7,8) = 0;
        x0(3,7) = 0; x0(5,7) = 0; x0(6,7) = 0;
        x0(3,6) = 0; x0(5,6) = 0;
        
        
    case 3
        x0(7,2) = 0; y0(18) = 1e-16;        % (x0(7,2) = Sim.kFA, y0(18) = Sim.fF)
        x0(8,2) = 0; y0(19) = 1e-16;        % (x0(8,2) = Sim.kGA, y0(19) = Sim.fG)
        
        % set intramolecular exchange rates of unused pools to zero
        x0(3,8) = 0; x0(5,8) = 0; x0(6,8) = 0; x0(7,8) = 0;
        x0(3,7) = 0; x0(5,7) = 0; x0(6,7) = 0;
        
    
    case 4
        x0(8,2) = 0; y0(19) = 1e-16;        % (x0(8,2) = Sim.kGA, y0(19) = Sim.fG)
        
        % set intramolecular exchange rates of unused pools to zero
        x0(3,8) = 0; x0(5,8) = 0; x0(6,8) = 0; x0(7,8) = 0;
        
        
end

nx = numel(xZspec);
w1 = gamma_2pi*x0(1,2);

% set y1 = y0 in the first derivation
if nargin<6             
    for k=1:1:nx
        y1(:,k) = y0;
    end;
end;

%% get intramolecular exchange rates

kbd = x0(3,5);
kdb = y0(14)/y0(16)*kbd;

kbe = x0(3,6);
keb = y0(14)/y0(17)*kbe;

kbf = x0(3,7);
kfb = y0(14)/y0(18)*kbf;

kbg = x0(3,8);
kgb = y0(14)/y0(19)*kbg;

kde = x0(5,6);
ked = y0(16)/y0(17)*kde;

kdf = x0(5,7);
kfd = y0(16)/y0(18)*kdf;

kdg = x0(5,8);
kgd = y0(16)/y0(19)*kdg;

kef = x0(6,7);
kfe = y0(17)/y0(18)*kef;

keg = x0(6,8);
kge = y0(17)/y0(19)*keg;

kfg = x0(7,8);
kgf = y0(18)/y0(19)*kfg;


%% get exchange rates with water pool (A)

kba = x0(3,2);              % (x0(3,2) = Sim.R1B)
kca = x0(4,2);              % (x0(4,2) = Sim.R1C)
kda = x0(5,2);              % (x0(5,2) = Sim.R1D)
kea = x0(6,2);              % (x0(6,2) = Sim.R1E)
kfa = x0(7,2);              % (x0(7,2) = Sim.R1F)
kga = x0(8,2);              % (x0(8,2) = Sim.R1G)

kab = y0(14)/y0(13)*kba;    % y0(14) = Sim.fB 
kac = y0(15)/y0(13)*kca;    % y0(15) = Sim.fC 
kad = y0(16)/y0(13)*kda;    % y0(16) = Sim.fD 
kae = y0(17)/y0(13)*kea;    % y0(17) = Sim.fE 
kaf = y0(18)/y0(13)*kfa;    % y0(18) = Sim.fF
kag = y0(19)/y0(13)*kga;    % y0(19) = Sim.fG


%% calculate transfer rates

k1a = x0(2,3) + kab + kac + kad + kae + kaf + kag;

k1b = x0(3,3) + kba + 0   + kbd + kbe + kbf + kbg;        % x0(3,3) = Sim.R1B
k1c = x0(4,3) + kca;                                      % x0(4,3) = Sim.R1C
k1d = x0(5,3) + kda + kdb + 0   + kde + kdf + kdg;        % x0(5,3) = Sim.R1D
k1e = x0(6,3) + kea + keb + ked + 0   + kef + keg;        % x0(6,3) = Sim.R1E
k1f = x0(7,3) + kfa + kfb + kfd + kfe + 0   + keg;        % x0(7,3) = Sim.R1F
k1g = x0(8,3) + kga + kgb + kgd + kge + kgf + 0  ;        % x0(8,3) = Sim.R1G

k2a = x0(2,4) + 0   + kab + kad + kae + kaf + kag;
k2b = x0(3,4) + kba + 0   + kbd + kbe + kbf + kbg;        % x0(3,4) = Sim.R2B
k2d = x0(5,4) + kda + kdb + 0   + kde + kdf + kdg;        % x0(5,4) = Sim.R2D
k2e = x0(6,4) + kea + keb + ked + 0   + kef + keg;        % x0(6,4) = Sim.R2E
k2f = x0(7,4) + kfa + kfb + kfd + kfe + 0   + keg;        % x0(7,4) = Sim.R2F
k2g = x0(8,4) + kga + kgb + kgd + kge + kgf + 0  ;        % x0(8,4) = Sim.R2G


%% matrix A with intramolecular exchange terms
                                                           
A = [   -k2a    kba     kda     kea     kfa     kga     0       0       0       0       0       0       0       0       0       0       0       0       0
        kab     -k2b    kdb     keb     kfb     kgb     0       0       0       0       0       0       0       0       0       0       0       0       0
        kad     kbd     -k2d    ked     kfd     kgd     0       0       0       0       0       0       0       0       0       0       0       0       0
        kae     kbe     kde     -k2e    kfe     kge     0       0       0       0       0       0       0       0       0       0       0       0       0
        kaf     kbf     kdf     kef     -k2f    kgf     0       0       0       0       0       0       0       0       0       0       0       0       0
        kag     kbg     kdg     keg     kfg     -k2g    0       0       0       0       0       0       0       0       0       0       0       0       0
        0       0       0       0       0       0       -k2a    kba     kda     kea     kfa     kga     w1      0       0       0       0       0       0
        0       0       0       0       0       0       kab     -k2b    kdb     keb     kfb     kgb     0       w1      0       0       0       0       0
        0       0       0       0       0       0       kad     kbd     -k2d    ked     kfd     kgd     0       0       0       w1      0       0       0
        0       0       0       0       0       0       kae     kbe     kde     -k2e    kfe     kge     0       0       0       0       w1      0       0
        0       0       0       0       0       0       kaf     kbf     kdf     kef     -k2f    kgf     0       0       0       0       0       w1      0
        0       0       0       0       0       0       kag     kbg     kdg     keg     kfg     -k2g    0       0       0       0       0       0       w1      
        0       0       0       0       0       0       -w1     0       0       0       0       0       -k1a    kba     kca     kda     kea     kfa     kga
        0       0       0       0       0       0       0       -w1     0       0       0       0       kab     -k1b    0       kdb     keb     kfb     kgb
        0       0       0       0       0       0       0       0       0       0       0       0       kac     0     -k1c      0       0       0       0
        0       0       0       0       0       0       0       0       -w1     0       0       0       kad     kbd     0       -k1d    ked     kfd     kgd
        0       0       0       0       0       0       0       0       0       -w1     0       0       kae     kbe     0       kde     -k1e    kfe     kge
        0       0       0       0       0       0       0       0       0       0       -w1     0       kaf     kbf     0       kdf     kef     -k1f    kgf
        0       0       0       0       0       0       0       0       0       0       0       -w1     kag     kbg     0       kdg     keg     kfg     -k1g
    ];
    
b = [ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; y0(13)*x0(2,3); y0(14)*x0(3,3); y0(15)*x0(4,3); y0(16)*x0(5,3); y0(17)*x0(6,3); y0(18)*x0(7,3); y0(19)*x0(8,3) ];
    
%%
wref=2*pi*x0(1,1);

% set 
wa = wref*x0(2,1);          % x0(2,1) = Sim.dwA
wb = wref*x0(3,1);          % x0(3,1) = Sim.dwB
wc = wref*x0(4,1);          % x0(4,1) = Sim.dwC
wd = wref*x0(5,1);          % x0(5,1) = Sim.dwD
we = wref*x0(6,1);          % x0(6,1) = Sim.dwE
wf = wref*x0(7,1);          % x0(7,1) = Sim.dwF
wg = wref*x0(8,1);          % x0(8,1) = Sim.dwG


if w1==0
    bb=0;
else
    bb=1;
end;

MT_lineshape  = z0(1);

%%

for k = 1:1:nx
    w = wref*xZspec(k);
    
    A(1,7)=-(w-wa)*bb;      % delta omega a
    A(7,1)=-(wa-w)*bb; 
    A(2,8)=-(w-wb)*bb;      % delta omega b
    A(8,2)=-(wb-w)*bb; 
    A(3,9)=-(w-wd)*bb;      % delta omega d
    A(9,3)=-(wd-w)*bb;  
    A(4,10)=-(w-we)*bb;     % delta omega e
    A(10,4)=-(we-w)*bb;     
    A(5,11)=-(w-wf)*bb;     % delta omega f
    A(11,5)=-(wf-w)*bb;  
    A(6,12)=-(w-wg)*bb;     % delta omega g
    A(12,6)=-(wg-w)*bb; 
    
    rfmt=RF_MT(1/x0(4,4),w1,w-wc,MT_lineshape);     %MT
    A(15,15) = -(k1c+rfmt);
    
    Ainvb=A\b;
    [ex diak] =expmdemo3(A*(tspan(2)-tspan(1)));
    dia(:,k)=diak/(tspan(2)-tspan(1));
    
    Z(:,k) = real(ex * (y1(:,k)+Ainvb)-Ainvb);
    
end

