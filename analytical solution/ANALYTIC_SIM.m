function [M Rex S] = ANALYTIC_SIM(Sim)
% comments here
% last change: 2014/04/03 by PS

Sim.td        = calc_td(Sim.tp,Sim.DC);     % td is always calculated form the duty-cycle and tp  


xZspec  = Sim.xZspec;
        
[S Rex] = rho_s(Sim);
       
R1rho   = S.Rho_full(:);        
theta   = S.theta(:);
Rex=Rex(:);
                
if strcmp(Sim.shape,'block')
    dEV=1;
    Pzeff=cos(theta).*dEV;
    Pz=cos(theta).*dEV;  
else
    Pzeff=1;
    Pz=1;  
end;

% calculate initial magnetisation
% to compensate pause pulse formalism relaxation time is just Sim.Trec-Sim.td
Zi = Sim.Zi;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CW SOLUTION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Sim.pulsed == 0
    
    if Sim.MT
        Z_cw = (Zi+cos(theta).*Pzeff.*S.R1obs./R1rho).*exp((R1rho*Sim.tp)) -cos(theta).*Pzeff.*S.R1obs./R1rho;
    else
        Z_cw = (Zi+cos(theta).*Pzeff.*Sim.R1A./R1rho).*exp((R1rho*Sim.tp)) -cos(theta).*Pzeff.*Sim.R1A./R1rho;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PULSED SOLUTIONS %%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else                    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% SANTYR SOLUTION %%%%%%%%
    %%%%%%%%% FOR n POOLS  %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Zss         = 1 - ((1-exp(R1rho(:)*Sim.tp)) .* (1-cos(theta(:)) .* Sim.R1A ./ (-R1rho(:)))) ./ (1-exp(R1rho(:).*Sim.tp-Sim.R1A*Sim.td));
    Z_santyr    = (Zi-Zss).*exp((R1rho(:)*Sim.tp-Sim.R1A*Sim.td)*Sim.n) +Zss;
    M.zspec_santyr = Z_santyr;
    clear Zss
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% approx pulsedCESL SOLUTION %%
    %%%%%%%%% FOR n POOLS  %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Q   = exp(-(Sim.R1A)*Sim.td);
    p   = exp(R1rho(:)*Sim.tp);
    if strcmp(Sim.shape,'block') || strcmp(Sim.shape,'SPINLOCK')
        alpha_F = alpha_f(Sim,xZspec);
    else
        Sim.alpha_fac = 0.5;
        alpha_F = Sim.alpha_fac* S.alpha_f_mean(:);    %alpha_f(P,xZspec);
    end
    Zss = (Pz.*Pzeff.*(1-Q) - (1-p.^(-1)) .* (-Pz.*cos(theta).*Sim.R1A./(R1rho))) ./ (p.^(-1)-Pz.*Pzeff.*((1-alpha_F(:)).*Q));
    xi  = ((1-alpha_F(:))*Q).*Pz.*Pzeff;   %%MZ 8.7.2013
    Z_pulsedCESL_approx = (Zi-Zss).*xi.^(Sim.n).*exp((R1rho(:)*Sim.tp)*Sim.n) +Zss;
    M.zspec_pulsedCESL_approx = Z_pulsedCESL_approx;
    clear Q p alpha_F Zss xi

    if Sim.n_cest_pool > 0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% pulsedCESL SOLUTION %%%%%% 
        %%% ONLY TWO POOL MODEL YET ! %%%
        %%% ALL OTHER POOLS = SANTYR  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Psi = (Sim.fB+Rex.minilorentz_b(:)./Sim.kBA);          %S.Rex ist negativ
        p   = exp(S.Rho_full(:)*Sim.tp);                     %S.Rho ist negativ
        r1a = Sim.R1A+Sim.kAB;
        r1b = Sim.R1B+Sim.kBA;
        a   = 1/2*(r1b+r1a+sqrt((r1b-r1a)^2+4*Sim.kBA*Sim.kAB));
        b   = 1/2*(r1b+r1a-sqrt((r1b-r1a)^2+4*Sim.kBA*Sim.kAB));
        dAA = 1/(a-b)*(-(b-r1a)*exp(-a*Sim.td)+(a-r1a)*exp(-b*Sim.td));
        dAB = Sim.kBA/(a-b)*(-exp(-a*Sim.td)+exp(-b*Sim.td));
        Zss = (1-(1-p.^(-1)).*(-Pz.*cos(theta).*Sim.R1A./(R1rho))-dAA-dAB*Sim.fB)./ (p.^(-1)-dAA-dAB*Psi);
        xi  = dAA./(1-dAB*Psi).*Pz.*Pzeff;

        Z_pulsedCESL = (Zi-Zss).*xi.^(Sim.n).*exp((R1rho*Sim.tp)*Sim.n) +Zss;
        M.zspec_pulsedCESL = Z_pulsedCESL;
        clear Psi p r1a r1b a b dAA dAB Zss xi

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% CERT SOLUTION %%%%%%%%% 
        %%% ONLY TWO POOL MODEL YET ! %%%
        %%% ALL OTHER POOLS = SANTYR  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         % pulse constants for Pool A
%         p   = exp(S.Rho_full*Sim.tp);
%         q   = (1-p).*(-Pz.*cos(theta).*Sim.R1A./(S.Rho_full));
% 
%         % pause constants
% %         r1a = Sim.R1A+Sim.kAB;
% %         r1b = Sim.R1B+Sim.kBA;
% %         a   = 1/2*(r1b+r1a+sqrt((r1b-r1a)^2+4*Sim.kBA*Sim.kAB));
% %         b   = 1/2*(r1b+r1a-sqrt((r1b-r1a)^2+4*Sim.kBA*Sim.kAB));
% 
%         % defining psi
%         alpha       = real(180/(2*pi)*(sqrt(4*(gamma_2pi*Sim.B1).^2-(Sim.R1B-Sim.R2B).^2))*Sim.tp);
%         Psi_CERT.M  = (Sim.fB+Rex.minilorentz_b./Sim.kBA); % für Moritz Lösung im Spinlock
%         Psi_CERT.A  = (Sim.fB+Rex.minilorentz_b./Sim.kBA).*(1-cosd(alpha).*exp(-(Sim.kBA+(Sim.R2B+Sim.R1B)/2)*Sim.tp)).*p;
%         Psi_CERT.B  = cosd(alpha)*exp(-(Sim.kBA+(Sim.R2B+Sim.R1B)/2)*Sim.tp);
%         Psi_CERT.C  = (Sim.fB+Rex.minilorentz_b./Sim.kBA).*(1-cosd(alpha).*exp(-(Sim.kBA+(Sim.R2B+Sim.R1B)/2)*Sim.tp)).*q;
% 
%         Mod_CERT    = calc_MTRsana_new(Sim,Pz,-S.Rho_full,Psi_CERT,theta);
%         Z_CERT      = Mod_CERT.dynamic(1,:);
%         M.zspec_CERT= Z_CERT;
%         clear p q alpha
        
    end
end; 
%%
% variables needed in PLOT_SPACE.m
M.S = S;
M.Znake=exp(R1rho*Sim.tp);        
M.R1rho=R1rho;
M.Reff=S.Reff_sincos;
M.x=xZspec;
M.Sim=Sim;
if Sim.pulsed
    M.zspec = Z_pulsedCESL_approx;
else  
    M.zspec=Z_cw;
end