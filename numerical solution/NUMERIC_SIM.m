function A = NUMERIC_SIM(Sim)
% comments here
% last change: 2015/06/11 by PS

xZspec      = Sim.xZspec;
anz_pulses  = Sim.n;

% because in adiabatic Spin-Lock case the pulse-shape needs to be changed between flip- and locking-pulses,
% the initial pulse-shape (set in BATCH-file) is saved as Sim.StartShape
Sim.StartShape = Sim.shape;
AdiaSamplingFactor = 100000;

if strcmp(Sim.shape,'SPINLOCK')
    SL1=1;
    SL2=1;
else
    SL1=0;
    SL2=0;
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
    teile = Sim.tp*AdiaSamplingFactor;
else
    teile = 200;
end;

StartInd = 0; %DKFZ PS: new running Index for adiabatic Spin-Lock simulations (starts at adiabatic flip and continues at locking and backflip)

tpulse = 0:Sim.tp/teile:Sim.tp;
%tpulse = 0:Sim.adiaDur/Sim.adiaTeile:Sim.adiaDur;

%% MV is the magnetisation vector at all offsets
MV_start    = BMsolution(x0, xZspec, [0 1], z0, y0);
MV       = MV_start.*Sim.Zi;

[w1 dphase theta B1q] = RF_pulse(Sim, teile);

for dd=1:Sim.dummies
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %%% jetzt setze ich hier noch den readout dran
    %%% der kommt vorneweg, dann ist das wie ein erster dummy von Zi ab
    B1rf = Sim.flipangle/(360*gamma_*0.0005); % B1=alpha/(gamma*trf)
    
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
    
    
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    % XXXXX calculations for different steps of the pulse train XXXXX
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    
    % loop for the different pulses (anz_pulses = Sim.n)
    for jj=1:fix(anz_pulses)
        
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXX PAUSE  / TRANSFER XXXXXXXXXXXXXXXXXXXXXX
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        % never start with a pause
        if jj > 1 && jj <= fix(anz_pulses)
            
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
        
   
    
    
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    % XXXXXXXXXXXXXXXXXX LABELING DURING PULSE XXXXXXXXXXXXXXXXXXXXXX
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    
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
    
    
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    % XXXXXXXXXXXXXXXXXXXXXXXXXX SPOILING XXXXXXXXXXXXXXXXXXXXXXXXXXX
    % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    % relaxtion during spoiling
    x0(1,2) = 0;
    MV = BMsolution(x0, xZspec, [0 0.0045], z0, y0, MV);
    
    MV(1:12,:)       = Sim.spoilf*MV(1:12,:);
    
    
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

try
    A.MVx = MVx;
    A.MVy = MVy;
    A.MVz = MVz;
catch
end

A.zspec_n=squeeze(MV_nth(:,13,:));
