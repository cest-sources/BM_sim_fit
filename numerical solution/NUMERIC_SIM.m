function A = NUMERIC_SIM(Sim)
% comments here
% last change: 2015/06/11 by PS
% NUMERIC_SIM solves the Bloch-McConnell equations stepwise for trains of
% arbitrary pulses.
% it requires:
% Sim.xZspec
% Sim.B1, Sim.n,Sim.tp,Sim.DC
% Sim.shape
% Sim.MT_lineshape
% pool concentrations
% y0 = [0 0 0 0 0 0 0 0 0 0 0 0 Sim.fA Sim.fB Sim.fC Sim.fD Sim.fE Sim.fF Sim.fG];
% pool parameters
% x0 = [  Sim.FREQ      Sim.B1        Sim.MT        Sim.n_cest_pool   0           0           0           0;
%         Sim.dwA       0             Sim.R1A       Sim.R2A           0           0           0           0;
%         Sim.dwB       Sim.kBA       Sim.R1B       Sim.R2B           Sim.kBD     Sim.kBE     Sim.kBF     Sim.kBG;
%         Sim.dwC       Sim.kCA       Sim.R1C       Sim.R2C           0           0           0           0;
%         Sim.dwD       Sim.kDA       Sim.R1D       Sim.R2D           0           Sim.kDE     Sim.kDF     Sim.kDG;
%         Sim.dwE       Sim.kEA       Sim.R1E       Sim.R2E           0           0           Sim.kEF     Sim.kEG;
%         Sim.dwF       Sim.kFA       Sim.R1F       Sim.R2F           0           0           0           Sim.kFG;
%         Sim.dwG       Sim.kGA       Sim.R1G       Sim.R2G           0           0           0           0
%         ];
% the initial magentization before saturation (M0 is always = 1)
% Sim.Zi
% the actual B1 pulse is read in from RF_pulse
% Sim.spoilf, this is the spoiling between pulses =0 full spoiling =1 no
% spoiing
%
% to model ideal readout the following parameters are introduced. However
% not used if dummies=0
% Sim.dummies
% Sim.flipangle
% Sim.readout
% Sim.TR
% Sim.linestomeasure

xZspec      = Sim.xZspec;
anz_pulses  = Sim.n;
Sim.td        = calc_td(Sim.tp,Sim.DC);     % td is always calculated form the duty-cycle and tp

% because in adiabatic Spin-Lock case the pulse-shape needs to be changed between flip- and locking-pulses,
% the initial pulse-shape (set in BATCH-file) is saved as Sim.StartShape
Sim.StartShape = Sim.shape;


if strcmp(Sim.shape,'SPINLOCK')
    SL1=1;
    SL2=1;
else
    SL1=0;
    SL2=0;
end

if strcmp(Sim.shape,'AdiaSL') || strcmp(Sim.shape,'AdiaSinCos');
    AdiaSL1=1;
    AdiaSL2=1;
else
    AdiaSL1=0;
    AdiaSL2=0;
end

if strcmp(Sim.shape,'AdiaInversion')
    Sim.B1 = Sim.AdiaB1;
    Sim.tp = Sim.AdiaDur;
end

%% manueller blochfit
z0(1) = {Sim.MT_lineshape};

y0 = [0 0 0 0 0 0 0 0 0 0 0 0 Sim.fA Sim.fB Sim.fC Sim.fD Sim.fE Sim.fF Sim.fG];

x0 = [  Sim.FREQ      Sim.B1        Sim.MT        Sim.n_cest_pool   0           0           0           0;
    Sim.dwA       0             Sim.R1A       Sim.R2A           0           0           0           0;
    Sim.dwB       Sim.kBA       Sim.R1B       Sim.R2B           Sim.kBD     Sim.kBE     Sim.kBF     Sim.kBG;
    Sim.dwC       Sim.kCA       Sim.R1C       Sim.R2C           0           0           0           0;
    Sim.dwD       Sim.kDA       Sim.R1D       Sim.R2D           0           Sim.kDE     Sim.kDF     Sim.kDG;
    Sim.dwE       Sim.kEA       Sim.R1E       Sim.R2E           0           0           Sim.kEF     Sim.kEG;
    Sim.dwF       Sim.kFA       Sim.R1F       Sim.R2F           0           0           0           Sim.kFG;
    Sim.dwG       Sim.kGA       Sim.R1G       Sim.R2G           0           0           0           0
    ];

x0(1,2) = 0;        % set B1 to zero
x00=x0;             % copy parameter

if (strcmp(Sim.shape,'block') || strcmp(Sim.shape,'SPINLOCK'))
    teile = 1;
elseif (strcmp(Sim.shape,'AdiaInversion'))
    teile = Sim.tp*Sim.AdiaSamplingFactor;
else
    teile = 200;
end;

StartInd = 0; %DKFZ PS: new running Index for adiabatic Spin-Lock simulations (starts at adiabatic flip and continues at locking and backflip)

tpulse = 0:Sim.tp/teile:Sim.tp;
%tpulse = 0:Sim.adiaDur/Sim.adiaTeile:Sim.adiaDur;

%% MV is the magnetisation vector at all offsets
MV_start    = BMsolution(x0, xZspec, [0 1], z0, y0);
MV       = MV_start.*Sim.Zi;


for dd=1:1+Sim.dummies
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %%% jetzt setze ich hier noch den readout dran
    %%% der kommt vorneweg, dann ist das wie ein erster dummy von Zi ab
    B1rf = Sim.flipangle/(360*gamma_*0.0005); % B1=alpha/(gamma*trf)
    if Sim.dummies>0 % only when dummies is >0 the readout will be played out
        if strcmp(Sim.readout,'bssfp')
            x0(1,2) = B1rf/2; % B1=alpha/(gamma*trf)  alpha halbe hier
            MV = BMsolution(x0, xZspec.*0, [0 0.0005], z0, y0, MV);
            
            for rr=1:fix(Sim.linestomeasure/2)
                x0(1,2)=-B1rf;
                MV = BMsolution(x0, xZspec.*0, [0 0.0005], z0, y0, MV);  %% flip
                x0(1,2) = 0;
                MV = BMsolution(x0, xZspec.*0, [0 Sim.TR], z0, y0, MV);  %% relax
                x0(1,2)=B1rf;
                MV = BMsolution(x0, xZspec.*0, [0 0.0005], z0, y0, MV);  %% backflip
                x0(1,2) = 0;
                MV = BMsolution(x0, xZspec.*0, [0 Sim.TR], z0, y0, MV);  %% relax
                % 2nd pulse of the TTM block
            end;
            x0(1,2) = -B1rf/2; % B1=alpha/(gamma*trf)  alpha halbe hier
            MV = BMsolution(x0, xZspec.*0, [0 0.0005], z0, y0, MV);
            
            
        elseif strcmp(Sim.readout,'gre')
            
            for rr=1:fix(Sim.linestomeasure)
                x0(1,2)=B1rf;
                MV = BMsolution(x0, xZspec.*0, [0 0.0005], z0, y0, MV);  %% flip
                MV(1:12,:) = 0;                                          %spoiling after TTM pulse block
                x0(1,2) = 0;
                MV = BMsolution(x0, xZspec.*0, [0 Sim.TR], z0, y0, MV);  %% relax
            end;
            
        elseif strcmp(Sim.readout,'se')  %% i assume that the train of 180 keeps the Z-magnetization relativley untouched, so i just leave away the TR relaxation.
            x0(1,2)=B1rf;
            MV = BMsolution(x0, xZspec.*0, [0 0.0005], z0, y0, MV);  %% flip
            
        end;
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        %Trec
        x0(1,2) = 0;
        MV          = BMsolution(x0, xZspec, [0 Sim.Trec], z0, y0, MV);
    end; % end readout phase when dummies>=1
    
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    % XXXXX calculations for different steps of the pulse train XXXXX
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    % calculate parameters for TTM pulses
    if Sim.TTM_rep >= 1;
        ttmtime = 1.*(1./8).*90.*pi./(gamma_2pi*180.*Sim.TTM_B1);
        binom_tau = 1E6 ./(2.*Sim.TTM_bandwidth.*Sim.FREQ.*1E6);
    end;
    
    % loop for the different pulses (anz_pulses = Sim.n)
    for jj=1:fix(anz_pulses)
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXX PAUSE  / TRANSFER XXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        % never start with a pause
        if jj > 1 && jj <= fix(anz_pulses)
            
            % TTM pulse block
            if Sim.TTM_rep >= 1;
                % apply TTM pulse before the current (jj) saturation pulse?
                if mod(jj, Sim.TTM_rep) == 0;
                    % how often should the TTM pulse be applied
                    for nn = 1:Sim.TTM_occurance;
                        % 1st pulse of the TTM block
                        x0(1,2) = sign(Sim.TTM_type(1))*Sim.TTM_B1;
                        MV = BMsolution(x0, xZspec.*0, [0 ttmtime*abs(Sim.TTM_type(1))], z0, y0, MV);
                        % pause between 1st and 2nd pulse
                        x0(1,2) = 0;
                        MV = BMsolution(x0, xZspec.*0, [0 (binom_tau - (abs(Sim.TTM_type(1)).*ttmtime + abs(Sim.TTM_type(2)).*ttmtime)./2)], z0, y0, MV);
                        % 2nd pulse of the TTM block
                        x0(1,2) = sign(Sim.TTM_type(2))*Sim.TTM_B1;
                        MV = BMsolution(x0, xZspec.*0, [0 ttmtime*abs(Sim.TTM_type(2))], z0, y0, MV);
                        % pause between 2nd and 3rd pulse
                        x0(1,2) = 0;
                        MV = BMsolution(x0, xZspec.*0, [0 (binom_tau - (abs(Sim.TTM_type(2)).*ttmtime + abs(Sim.TTM_type(3)).*ttmtime)./2)], z0, y0, MV);
                        % 3rd pulse of the TTM block
                        x0(1,2) = sign(Sim.TTM_type(3))*Sim.TTM_B1;
                        MV = BMsolution(x0, xZspec.*0, [0 ttmtime*abs(Sim.TTM_type(3))], z0, y0, MV);
                        % pause between 3rd and 4th pulse
                        x0(1,2) = 0;
                        MV = BMsolution(x0, xZspec.*0, [0 (binom_tau - (abs(Sim.TTM_type(3)).*ttmtime + abs(Sim.TTM_type(4)).*ttmtime)./2)], z0, y0, MV);
                        % 4th pulse of the TTM block
                        x0(1,2) = sign(Sim.TTM_type(4))*Sim.TTM_B1;
                        MV = BMsolution(x0, xZspec.*0, [0 ttmtime*abs(Sim.TTM_type(4))], z0, y0, MV);
                        
                        %spoiling after TTM pulse block
                        MV(1:12,:) = 0;
                    end
                end
            end
            
            % relaxtion during pause
            x0(1,2) = 0;
            MV = BMsolution(x0, xZspec, [tpulse(numel(w1)) tpulse(numel(w1))+Sim.td], z0, y0, MV);
            
        end;
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXX FLIP MAGNETISATION INTO EFFECTIVE FRAME (SPINLOCK CASE) XXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        if SL1
            thetA = atan(-Sim.B1./(xZspec.*Sim.FREQ/gamma_));
            thetB = thetA;
            thetC = thetA;
            thetD = thetA;
            thetE = thetA;
            thetF = thetA;
            thetG = thetA;
            
            for iter=1:numel(thetA)
                ySL = [ sin(thetA(iter)).*MV(13,iter);
                    sin(thetB(iter)).*MV(14,iter);
                    sin(thetD(iter)).*MV(16,iter);
                    sin(thetE(iter)).*MV(17,iter);
                    sin(thetF(iter)).*MV(18,iter);
                    sin(thetG(iter)).*MV(19,iter);
                    0;
                    0;
                    0;
                    0;
                    0;
                    0;
                    cos(thetA(iter)).*MV(13,iter);
                    cos(thetB(iter)).*MV(14,iter);
                    cos(thetC(iter)).*MV(15,iter);
                    cos(thetD(iter)).*MV(16,iter);
                    cos(thetE(iter)).*MV(17,iter);
                    cos(thetF(iter)).*MV(18,iter);
                    cos(thetG(iter)).*MV(19,iter)
                    ];
                MV(:,iter)=ySL;
            end
        end
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXX ADIABATIC SPINLOCK XXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXX
        
        if AdiaSL1
            x0(1,2) = Sim.AdiaB1; % set B1 to correct value for adiabatic pulse
            Sim.shape = Sim.StartShape; % set correct pulse-shape (is changed for locking-pulse later)
            Sim.tp = Sim.AdiaDur; % set duration to correct value for adiabatic pulse
            Sim.AdiaAmpFactor = 1; %increasing amplitude
            Sim.AdiaFreqFactor = 1;
            AdiaTeile = Sim.AdiaDur*Sim.AdiaSamplingFactor; % calculate new "AdiaTeile" for adiabatic pulse
            tpulse=0:Sim.tp/AdiaTeile:Sim.tp;
            [w1 dphase theta B1q] = RF_pulse(Sim, AdiaTeile);
            B1ModAdia1 = w1; % stored to check for correct amp. modulation
            FreqModAdia1 = dphase; % stored to check for correct freq. modulation
            for ii = 1:numel(w1)-1
                % set correct B1 value in every pulse part/step
                x0(1,2)=w1(ii);
                
                %frequency shift of every pulse step (adiabatic pulse has amp. and freq. modulation)
                %             x0(2,1)=x00(2,1)+dphase(ii);
                %             x0(3,1)=x00(3,1)+dphase(ii);
                %             x0(4,1)=x00(4,1)+dphase(ii);
                %             x0(5,1)=x00(5,1)+dphase(ii);
                %             x0(6,1)=x00(6,1)+dphase(ii);
                %             x0(7,1)=x00(7,1)+dphase(ii);
                %             x0(8,1)=x00(8,1)+dphase(ii);
                
                adiaxZspec=xZspec;
                adiaxZspec(xZspec>=0)=xZspec(xZspec>=0)-dphase(ii);
                adiaxZspec(xZspec<0)=xZspec(xZspec<0)+dphase(ii);
                
                
                [MV A_temp Ainvb_temp dia] = BMsolution(x0, adiaxZspec, [tpulse(ii) tpulse(ii+1)], z0, y0, MV);
                MVx(ii) = MV(1);
                MVy(ii) = MV(7);
                MVz(ii) = MV(13);
            end;
            % set parameters for the Locking-Pulse
            Sim.shape = 'block';
            Sim.tp = Sim.AdiaTSL;
            teile = Sim.AdiaTSL*Sim.AdiaSamplingFactor; % should be the same factor as for AdiaTeile before
            teile=1; % the block pulse is not discretisized
            clearvars dphase w1
            tpulse=0:Sim.tp/teile:Sim.tp;
            [w1 dphase theta B1q] = RF_pulse(Sim, teile);
            B1ModLock = ones(1,Sim.AdiaTSL*Sim.AdiaSamplingFactor).*Sim.B1;
            FreqModLock = zeros(1,Sim.AdiaTSL*Sim.AdiaSamplingFactor);
        end
        
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXX LABELING DURING PULSE XXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        [w1 dphase theta B1q] = RF_pulse(Sim, teile);
        
        if AdiaSL1
            StartInd = AdiaTeile;
        end
        for ii = 1:numel(w1)-1
            % set correct B1 value in every pulse part/step
            x0(1,2)=w1(ii);
            
            %frequency shift of every pulse step
            if numel(dphase)>1
                x0(2,1)=x00(2,1)+dphase(ii);
                x0(3,1)=x00(3,1)+dphase(ii);
                x0(4,1)=x00(4,1)+dphase(ii);
                x0(5,1)=x00(5,1)+dphase(ii);
                x0(6,1)=x00(6,1)+dphase(ii);
                x0(7,1)=x00(7,1)+dphase(ii);
                x0(8,1)=x00(8,1)+dphase(ii);
            end;
            
            [MV A_temp Ainvb_temp dia] = BMsolution(x0, xZspec, [tpulse(ii) tpulse(ii+1)], z0, y0, MV);
            
            MVx(StartInd+ii) = MV(1); 	%StartInd != 0 for adiabatic SL, else its 0
            MVy(StartInd+ii) = MV(7);	%StartInd != 0 for adiabatic SL, else its 0
            MVz(StartInd+ii) = MV(13);	%StartInd != 0 for adiabatic SL, else its 0
            
        end;
        
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXX FLIP MAGNETISATION BACK (SPINLOCK CASE) XXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        if SL2
            for iter = 1:numel(thetA)
                MV(:,iter)=[ cos(thetA(iter)).*MV(1,iter)-sin(thetA(iter)).*MV(13,iter);
                    cos(thetB(iter)).*MV(2,iter)-sin(thetB(iter)).*MV(14,iter);
                    cos(thetD(iter)).*MV(3,iter)-sin(thetD(iter)).*MV(16,iter);
                    cos(thetE(iter)).*MV(4,iter)-sin(thetE(iter)).*MV(17,iter);
                    cos(thetF(iter)).*MV(5,iter)-sin(thetF(iter)).*MV(18,iter);
                    cos(thetG(iter)).*MV(6,iter)-sin(thetG(iter)).*MV(19,iter);
                    MV(7,iter);
                    MV(8,iter);
                    MV(9,iter);
                    MV(10,iter);
                    MV(11,iter);
                    MV(12,iter);
                    sin(thetA(iter)).*MV(1,iter)+cos(thetA(iter)).*MV(13,iter);
                    sin(thetB(iter)).*MV(2,iter)+cos(thetB(iter)).*MV(14,iter);
                    MV(15,iter);
                    sin(thetD(iter)).*MV(3,iter)+cos(thetD(iter)).*MV(16,iter);
                    sin(thetE(iter)).*MV(4,iter)+cos(thetE(iter)).*MV(17,iter);
                    sin(thetF(iter)).*MV(5,iter)+cos(thetF(iter)).*MV(18,iter);
                    sin(thetG(iter)).*MV(6,iter)+cos(thetG(iter)).*MV(19,iter)
                    ];
            end
        end
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXX ADIABATIC SPINLOCK XXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXX
        
        if AdiaSL2
            StartInd = StartInd + teile; % set correct Index after locking-pulse
            x0(1,2) = Sim.AdiaB1; % set to adiabatic B1 again
            Sim.shape = Sim.StartShape; % set to initial pulse-shape
            Sim.tp = Sim.AdiaDur;
            Sim.AdiaAmpFactor = 0; % decreasing amplitude for backflip
            Sim.AdiaFreqFactor = -1;
            tpulse=0:Sim.tp/AdiaTeile:Sim.tp;
            [w1 dphase theta B1q] = RF_pulse(Sim, AdiaTeile);
            B1ModAdia2 = w1; % stored to check for correct amp. modulation
            FreqModAdia2 = dphase; % stored to check for correct freq. modulation
            for ii = 1:numel(w1)-1
                
                % set correct B1 value in every pulse part/step
                x0(1,2)=w1(ii);
                %             x0(2,1)=x00(2,1)+dphase(ii);
                %             x0(3,1)=x00(3,1)+dphase(ii);
                %             x0(4,1)=x00(4,1)+dphase(ii);
                %             x0(5,1)=x00(5,1)+dphase(ii);
                %             x0(6,1)=x00(6,1)+dphase(ii);
                %             x0(7,1)=x00(7,1)+dphase(ii);
                %             x0(8,1)=x00(8,1)+dphase(ii);
                
                adiaxZspec=xZspec;
                adiaxZspec(xZspec>=0)=xZspec(xZspec>=0)-dphase(ii);
                adiaxZspec(xZspec<0)=xZspec(xZspec<0)+dphase(ii);
                
                
                [MV A_temp Ainvb_temp dia] = BMsolution(x0, adiaxZspec, [tpulse(ii) tpulse(ii+1)], z0, y0, MV);
                
                MVx(StartInd + ii) = MV(1);
                MVy(StartInd + ii) = MV(7);
                MVz(StartInd + ii) = MV(13);
            end;
        end
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXX SPOILING XXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        % relaxtion during spoiling
        x0(1,2) = 0;
        MV = BMsolution(x0, xZspec, [0 0.0045], z0, y0, MV);
        
        MV(1:12,:)       = Sim.spoilf*MV(1:12,:);
        
        
        %% relaxation post saturation tpost
        %     if isfield(Sim,'tpost')
        %     MV          = BMsolution(x0, xZspec, [0 Sim.tpost], z0, y0, MV);
        %     end;
        
        % hier ist nach einem LTM (label transfer module)
        MV_nth(jj,:,:)  = MV;
        
        
    end; % LTMs
    
end; % dummies
A.Sim = Sim;
A.x=xZspec';
A.zspec=squeeze(MV(13,:))';
% A.zspecb=squeeze(MV(14,:))';
% A.zspecc=squeeze(MV(15,:))';
% A.zspecd=squeeze(MV(16,:))';
% A.zspece=squeeze(MV(17,:))';
% A.zspecf=squeeze(MV(18,:))';
% A.zspecg=squeeze(MV(19,:))';
%
% A.Mx=squeeze(MV(1,:))';
% A.My=squeeze(MV(7,:))';
%
% A.Mxb=squeeze(MV(2,:))';
% A.Myb=squeeze(MV(8,:))';
%
% A.Mxd=squeeze(MV(3,:))';
% A.Myd=squeeze(MV(9,:))';
%
% A.Mxe=squeeze(MV(4,:))';
% A.Mye=squeeze(MV(10,:))';
%
% A.Mxf=squeeze(MV(5,:))';
% A.Myf=squeeze(MV(11,:))';
%
% A.Mxg=squeeze(MV(6,:))';
% A.Myg=squeeze(MV(12,:))';

% parameter from RF_pulse and BMsolution
A.dia=dia;      %eigenvalues
A.theta=theta;
A.B1q=B1q;
if AdiaSL1
    try
        A.B1Mod = [B1ModAdia1 B1ModLock B1ModAdia2];
        A.FreqMod = [FreqModAdia1 FreqModLock FreqModAdia2];
    catch
    end
end
try
    A.MVx = MVx;
    A.MVy = MVy;
    A.MVz = MVz;
catch
end

A.zspec_n=squeeze(MV_nth(:,13,:));
