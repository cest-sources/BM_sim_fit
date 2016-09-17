%% this file allows to create simulated Z-spectrum data to check QUESP and QUEST formulas
%%it also allows to fit single offset QUESP and QUEST formulaes to the
%%(1) simulated data
%%(2) CESTtab structured data (go directly to 2 if you want to do so)

%% 1.1: create paraemter struct P for simulation or fit
% simulation parameters

P.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)

P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
P.MT_lineshape  = 'Gaussian';           % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
P.n_cest_pool   = 1;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'


% sequence/scanner parameters
% readout
P.TR  = 3/1000;
P.linestomeasure=1;
P.flipangle= 14;
P.readout='bssfp';
P.dummies=1;  %% or better shots?
P.TTM_rep         = 0;    

% saturation
P.FREQ          = 7*gamma_;           % frequency (=B0[T] * gamma)
P.B1            = 10;                   % B1 value in µT
P.Trec          = 3;                    % recover time in s
P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
P.Zi            = 0;                    % initial magnetisation (should be between -1 and +1)

P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation

if P.pulsed
    P.shape = 'gauss_ucl';              % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
    P.n     = 151;                      % number of saturation pulses
    P.tp    = 0.05;                      % saturation time per pulse in s
    P.DC    = 0.9091;                   % duty cycle
else
    P.tp    = 20;                       % saturation time in s
    P.n     = 1;                        % choose n=1 for cw saturation
    P.shape = 'block';                  % choose 'block' for cw saturation
    P.DC    = 1.0;                      % choose DC=1 for cw saturation
end;

%XXXXXXXXXXXXXXX DONT CHANGE ! XXXXXXXXXXXXXXXX
P.B1cwpe_quad   = -1;                     % this is the B1 power equivalent flag: -1 means that P.B1 is the average B1 over a single pulse 
P.td            = calc_td(P.tp,P.DC);     % td is always calculated form the duty-cycle and tp  
P.c             = 1;                      % this is a free parameter to play around in the solutions  
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% 1.2: set t vary and varyval (manuel)
% which variables were varied in your measurements
%% 1.2a:QUEST
      vary    = {'tp'}; % e.g. {'B1','Trec','DC'}
      clear varyval;
      varyval(1,:) = [10.0000    7.5000    5.0000    4.0000    3.0000    2.0000    1.5000    1.0000    0.7500    0.5000 ]
   
%% 1.2b:QUESP
      vary    = {'B1'};
      clear varyval;
       varyval(1,:) = [1 5 10 15 20 25 30 35]; % e.g. for B1: [0.5 0.75 1 2]
%       varyval(1,:) = [0.2:0.4:3]; % e.g. for B1: [0.5 0.75 1 2]

%% 1.3 create simulated data, adds noise and set depending variables and guess good starting values for the fit
% here the simulation data is created and starting values can be guessed
 w_x=[-80:2:-60 -59:-40 -38:2:38 40:60 62:2:80]';  % offset axis
%  w_x=[-6:0.1:6];  % offset axis
 w_xx=repmat(w_x,numel(varyval),1);                 % combined offset axis
 P.xZspec=w_x;
 
%% set start values
P.analytic=1;

tic
P.normalized=[]
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));


% get start parameters for tissue CEST-agent
P.tissue    = 'PBS_telaviv';
P.CESTagent = 'glucose_ta';
P           = getSim(P);    
expname=sprintf('Simulated data tissue:%s, agent:%s, kBA: %.2f, fB: %7f',P.tissue,P.CESTagent,P.kBA,P.fB);
 
% P.dwB=[10]
% warning(sprintf('P.dwB is at offset %.2f ppm',P.dwB));

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL A XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyA       = [   1          0           1       ];
T.dep_varsA   = {   'dwA',      'R1A',      'R2A'   };                  % {'dwA','R1A','R2A'};
T.startA      = [   P.dwA       P.R1A       P.R2A   ];                  % [P.dwA P.R1A P.R2A];
T.lowerA      = [   -1        1/4*P.R1A     1/3*P.R2A     ];
T.upperA      = [   +1        2*P.R1A       20000*P.R2A     ];

[T.dep_varsA, T.startA, T.lowerA, T.upperA] = selectVars( T.varyA, T.dep_varsA, T.startA, T.lowerA, T.upperA );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL B XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyB       = [   1           1           1           1           0           0           0           0           ];
T.dep_varsB   = {   'dwB',      'fB',       'kBA',      'R2B',      'kBD',      'kBE',      'kBF',      'kBG'       };     
T.startB      = [   P.dwB       P.fB        P.kBA       P.R2B       P.kBD       P.kBE       P.kBF       P.kBG       ];
T.lowerB      = [   P.dwB-15         P.fB*0.1     50          1/20        0           0           0           0           ];
T.upperB      = [   P.dwB+15         P.fB*10        1500000        66         10000       10000       10000       10000       ];

[T.dep_varsB, T.startB, T.lowerB, T.upperB] = selectVars( T.varyB, T.dep_varsB, T.startB, T.lowerB, T.upperB );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL C XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% 
% T.varyC       = [ 0         0           0           0           0           ];
% T.dep_varsC   = {'dwC',     'fC',       'kCA',      'R1C',      'R2C'       };     
% T.startC      = [P.dwC      P.fC        P.kCA       P.R1C       P.R2C       ];
% T.lowerC      = [0          0.0001      10          1/20        1/20        ];
% T.upperC      = [9.9        0.1         10000       1/0.1       1/0.1       ];
% 
% [T.dep_varsC, T.startC, T.lowerC, T.upperC] = selectVars( T.varyC, T.dep_varsC, T.startC, T.lowerC, T.upperC );
% 
% 

% % XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% % XXXXXXXXXX POOL D XXXXXXXXXX
% % XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% 
% % 
T.varyD       = [ 1         1           1           1           0           0           0           ];
T.dep_varsD   = {'dwD',     'fD',       'kDA',      'R2D',      'kDE',      'kDF',      'kDG'       };     
T.startD      = [P.dwD      P.fD        P.kDA       P.R2D       P.kDE       P.kDF       P.kDG       ];
T.lowerD      = [0        0.0000001     10          1/20        0           0           0           ];
T.upperD      = [3        0.01        50000         20         10000       10000       10000       ];
 
[T.dep_varsD, T.startD, T.lowerD, T.upperD] = selectVars( T.varyD, T.dep_varsD, T.startD, T.lowerD, T.upperD )
% 


% determine which case (--> parameters) is used
[dep_vars startValue lowerbounds upperbounds] = casedetermination(P,T);


% plot your measured and simulated data (uses starting point values)
if P.analytic
        GUESS = conv_ana(startValue,w_x,P,dep_vars,vary,varyval);
else
        GUESS = conv_num(startValue,w_x,P,dep_vars,vary,varyval);
end

% plot "1st-guess" with deviations from data points
figure(2001),

GUESSmat=reshape(GUESS,numel(w_x),numel(varyval));
plot(w_x,GUESSmat,'.-');
clear leg;

for ii=1:numel(startValue)
    lab{ii} = (sprintf('%s=%.4f',dep_vars{ii},startValue(ii)));
end; 

for ii=1:numel(varyval)
     leg{ii} = (sprintf('%s=%.2f',vary{1},varyval(ii)));
end; 
legend(leg);
xlabel(lab);
title(expname);

set(gca,'XDir','reverse');

%% 1.4 here the simulation data  is noised and is used as input data (overwrites Z)
Z_xx = ricernd(GUESS',0.001); % std value is 0.005
Z_xx=Z_xx';                                     % Z-spectrum of combined offset axis
Z_x=reshape(Z_xx,numel(w_x),numel(varyval));    % Z-spectrum matrix of offset axis (Thats teh one you can plot by plot(w_x,Z_x)
figure('Name','Noisy data'); plot(w_x,Z_x); title(sprintf('Noisy data, varying %s',vary{1}));
%% 2 QUEST and QUESP on experimental data  
%%2.1: pick QUEST and QUESP data from CESTtab structure
%first load a Ztab
Ztab=Ztab_fullglint;

rowname='ta2';                                 % define  the row you want to evaluate    
expname=Ztab{rowname,'exp'}{1};  

P.normalized=[-4.6]
 Ztab(rowname,:) = norm_run(Ztab(rowname,:),'B1_run',P.normalized) 
 Ztab(rowname,:) = exclude_run(Ztab(rowname,:),'B1_run',-5) 
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));

[w_x, Z_x, w_xx, Z_xx,  varyval, vary, Ptab]= plot_tab(Ztab,rowname,'B1_run');

P = catstruct(P,Ptab);
%% 3 fitting of the Z_x data (wherever it comes from)

%% 3.1 RUN full BM OPTIMIZATION 
% you need a startvalue, run 1.3 first!
P.analytic=1;  % set this to 1 if analytic fit should be used, numeric =0 can take forever
P.n_cest_pool=2;

try close 2001 % close "1st-guess"
end
if P.analytic == 1 % use analytical solution
        f = @(x,xdata) conv_ana(x,w_x,P,dep_vars,vary,varyval);
else % use numerical solution
        f = @(x,xdata) conv_num(x,w_x,P,dep_vars,vary,varyval);
end

% fit-options
options = optimset('PlotFcns',{@optimplotfval},'TolFun',1e-8,'MaxIter',400,'Display','on');

    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,w_xx,Z_xx,lowerbounds,upperbounds,options);

[ci, varb, corrb, varinf] = nlparci(x,RES,JACOBIAN,0.95);
ci(:,1) = ci(:,1)-x';
ci(:,2) = ci(:,2)-x';
figure, imagesc(abs(corrb))
xlabel(dep_vars)
ylabel(dep_vars)

FIT=f(x,w_x);

for i=1:numel(dep_vars)
    fprintf('%s: real=%f  fit=%f±%.5f  start=%f\n',dep_vars{i},P.(dep_vars{i}),x(i),ci(i,2), startValue(i));
    end;

clear leg FITres;
figure('Name',expname)

Z_x=reshape(Z_xx,numel(w_x),numel(varyval));
FITmat=reshape(FIT,numel(w_x),numel(varyval));
plot(w_x,Z_x,'.',w_x,FITmat,w_x, Z_x-FITmat);
title(expname);
[Label, Unit]=getSimLabelUnit()
for ii=1:numel(startValue) 
    leg{ii} = (sprintf('%s=%.7f±%.7f %s (start: %.7f)',Label.(dep_vars{ii}),x(ii),ci(ii,2),Unit.(dep_vars{ii}),startValue(ii)));
    leg{ii} = (sprintf('%s=%s %s (start: %.7f)',Label.(dep_vars{ii}),errorbar_str(x(ii),ci(ii,2)),Unit.(dep_vars{ii}),startValue(ii)));
            
    FITres.(dep_vars{ii})=x(ii)+j*ci(ii,2);
end; 
xlabel(leg);
set(gca,'XDir','reverse');

FITres.P=P;
FITres.T=T;
FITres.dep_vars=dep_vars;
FITres.startValue=startValue;
FITres.lowerbounds=lowerbounds;
FITres.upperbounds=upperbounds;
FITres.w_xx=w_xx;
FITres.Z_xx=Z_xx;
FITres.w_x=w_x;
FITres.Z_x=Z_x;
Fitres.varyval=varyval;
Fitres.vary=vary;

%% save fitresult in Ztable
Ztab(rowname,'FITres')={{FITres}};
%check
Ztab{rowname,'FITres'}{1}.kBA


%% 3.2 single offset preparation

single_offset=x(strcmp(dep_vars,'dwB')>0);

single_offset=3;


clear leg
figure(42), plot(w_x,Z_x); title(sprintf('check data again! The chosen offset is %.2f ppm',single_offset));
set(gca,'XDir','reverse');
for ii=1:numel(varyval)
     leg{ii} = (sprintf('%s=%.2f',vary{1},varyval(ii)));
end; 
legend(leg);

ind_lab=find_nearest(w_x,single_offset);
ind_ref=find_nearest(w_x,-single_offset);
Zlab=Z_x(ind_lab,:);
Zref=Z_x(ind_ref,:);

 [S Rex] = rho_s(P);
-Rex.minilorentz_b(ind_lab)
-Rex.Hyper_full(ind_lab)

alpha=-Rex.minilorentz_b(ind_lab)/(P.fB*P.kBA)

%% 3.3 single offset QUEST ( using MTRasym) 
% MTR for that single_offset
if strcmp(vary{1},'tp')==0
    error('It seems you varied %s, you need to vary tp for QUEST, or try QUESP',vary{1});
end;
for ii=1:4
    if ii==1
MTR=Zref-Zlab;
% original QUEST

modelstr= @(fb,kb,R1,w1,x) fb*kb/(R1+fb*kb)*(w1)^2/(w1^2+kb^2)*(1- exp(-(R1+fb*kb)*x) );

    elseif ii==2
MTR=Zref-Zlab;
       
% revised QUEST
modelstr= @(fb,kb,R1,w1,x) fb*kb*w1^2/(w1^2+kb^2)/(R1+fb*kb)*(1- exp(-(R1+fb*kb*w1^2/(w1^2+kb^2))*x) );


    elseif ii==3       
% full revised QUEST
MTR=Zref-Zlab;
modelstr= @(fb,kb,R1,w1,x) fb*kb*w1^2/(w1^2+kb^2)/(R1+fb*kb*w1^2/(w1^2+kb^2))*(1- exp(-(R1+fb*kb*w1^2/(w1^2+kb^2))*x) );
 
    elseif ii==4    
        % full revised QUEST with arbitrary initial conditions

        MTR=Zref-Zlab;     
        if isempty(P.normalized)
modelstr= @(fb,kb,R1,w1,x) fb*kb*w1^2/(w1^2+kb^2)/(R1+fb*kb*w1^2/(w1^2+kb^2))-(P.Zi- R1/(R1+fb*kb*w1^2/(w1^2+kb^2)))*exp(-(R1+fb*kb*w1^2/(w1^2+kb^2))*x)+(P.Zi-1)*exp(-R1*x);
        else
            modelstr= @(fb,kb,R1,w1,x) (fb*kb*w1^2/(w1^2+kb^2)/(R1+fb*kb*w1^2/(w1^2+kb^2))-(P.Zi- R1/(R1+fb*kb*w1^2/(w1^2+kb^2)))*exp(-(R1+fb*kb*w1^2/(w1^2+kb^2))*x)+(P.Zi-1)*exp(-R1*x))./(1+(P.Zi-1)*exp(-R1*x));
         end;
    end;
    
 
[xData, yData] = prepareCurveData( varyval, MTR );

% Set up fittype and options.
ft = fittype(modelstr, 'problem', {'R1','w1'}, 'independent', {'x'}, 'dependent', 'y' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.000135 0];
opts.StartPoint = [0.000135 4000];
opts.Upper = [0.000135 15000];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A P.B1*gamma_2pi} );

% Plot fit with data.
figure(2);
subplot(4,1,ii);

h = plot( fitresult, xData, yData ); hold on;
ci = confint(fitresult);
ci = ci(2,:)-coeffvalues(fitresult);
legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,vary{1}), sprintf('QUEST fit 1, fb=%.6f±%.9f, kb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );
title(func2str(modelstr));
% Label axes
xlabel( vary );
ylabel( 'MTR' );
grid on


end;
%% 3.4 single offset QUESP(using MTRasym) 
if strcmp(vary{1},'B1')==0
    error('It seems you varied %s, you need to vary B1 for QUEST',vary{1});
end;

for ii=1:6
    if ii==1
MTR=Zref-Zlab;
% original QUESP

modelstr= @(fb,kb,R1,x,w1)  fb.*kb./(R1+fb.*kb).*w1.^2./(w1.^2+kb.^2).*(1- exp(-(R1+fb.*kb).*x) );
         
    elseif ii==2
MTR=Zref-Zlab;
%fully revised
            
modelstr= @(fb,kb,R1,x,w1)  fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*(1- exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x) );

    
    elseif ii==3
        % general QUESP(t)
        if isempty(P.normalized)
            modelstr= @(fb,kb,R1,x,w1) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x);
        else
              modelstr= @(fb,kb,R1,x,w1) (fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x))./(1+(P.Zi-1).*exp(-R1.*x));
        end;
        
    elseif ii==4
        % general QUESP(t) with Rexmean
        
        if isempty(P.normalized)
            modelstr= @(fb,kb,R1,x,w1) Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi)./(R1+Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi))-(P.Zi- R1./(R1+Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi))).*exp(-(R1+Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi)).*x)+(P.Zi-1).*exp(-R1.*x);
        else
              modelstr= @(fb,kb,R1,x,w1) (Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi)./(R1+Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi))-(P.Zi- R1./(R1+Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi))).*exp(-(R1+Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi)).*x)+(P.Zi-1).*exp(-R1.*x))./(1+(P.Zi-1).*exp(-R1.*x));
        end;
        
    elseif ii==5
        
        % inverse QUESP
        MTR=1./Zlab-1./Zref;
        
        modelstr= @(fb,kb,R1,x,w1) fb.*kb.*w1.^2./(w1.^2+kb.^2)./R1*x/x;
        if P.pulsed
            c1= 1/(2*2.92)*sqrt(2*pi);
            c22=(c1*sqrt(sqrt(2))).^2;
            %  c1=0.5
            %  c22=0.35;
            modelstr= @(fb,kb,R1,x,w1) P.DC.*c1.*fb.*kb.*w1.^2./(w1.^2+kb.^2.*c22)./R1*x/x;
        end;
        % .*(1-exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x));
        
    elseif ii==6
        MTR=1./Zlab-1./Zref;
       % WIP use fully integrated Rex
        modelstr= @(fb,kb,R1,x,w1) x./x*P.DC./R1.*Rex_mean(P,P.R2B,fb,kb,w1./gamma_2pi) ;
        
    end;

[xData, yData] = prepareCurveData( varyval*gamma_2pi, MTR );

% Set up fittype and options.
ft = fittype(modelstr, 'problem', {'R1','x'}, 'independent', {'w1'}, 'dependent', 'y' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.0000135 0];
opts.StartPoint = [0.000135 4000];
opts.Upper = [0.0135 150000];

% Fit model to data.
if P.pulsed
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A (P.tp*P.n/P.DC)} );
else
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A P.tp} );
end;

% Plot fit with data.
figure(1);
subplot(6,1,ii);
h = plot( fitresult, xData, yData ); hold on;
ci = confint(fitresult);
ci = ci(2,:)-coeffvalues(fitresult);
legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,vary{1}), sprintf('QUESP fit 1, fb=%.2e±%.1e, kb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );
title(func2str(modelstr));
% Label axes
xlabel( vary );
ylabel( 'MTR' );
grid on

fb=0.0009;
end;

%% add omegaplot

% mine
for ii=1:3
    
if ii==1
MTR=1./Zlab-1./Zref;
ylab='1/MTR_{Rex}';
modelstr= @(fb,kb,R1,xx) R1*(1./(fb.*kb) + kb./fb.*xx )
[xData, yData] = prepareCurveData( 1./(varyval*gamma_2pi).^2, 1./MTR );

elseif ii==2
   % original from dixons paper
    MTR=Zlab./(1-Zlab);
    ylab='MTR_{Dixon}';
    modelstr= @(fb,kb,R1,xx) R1*(1./(fb.*kb) + kb./fb.*xx )
[xData, yData] = prepareCurveData( 1./(varyval*gamma_2pi).^2, MTR );

else 
   % pulsed form factors
    MTR=1./Zlab-1./Zref;
    ylab='MTR_{Rex,pulsed}';
    c1= 1/(2*2.92)*sqrt(2*pi);
    c22=(c1*sqrt(sqrt(2))).^2;
    
    modelstr= @(fb,kb,R1,xx) R1./(P.DC.*c1)*(1./(fb.*kb) + c22.*kb./fb.*xx );
[xData, yData] = prepareCurveData( 1./(varyval*gamma_2pi).^2, 1./MTR );

end;


% Set up fittype and options.
ft = fittype(modelstr, 'problem', {'R1'}, 'independent', {'xx'}, 'dependent', 'y' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.0000135 0];
opts.StartPoint = [0.000135 2000];
opts.Upper = [0.00135 150000];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A} );

% Plot fit with data.
figure(11);
subplot(3,1,ii);
h = plot( fitresult, xData, yData ); hold on;
ci = confint(fitresult);
ci = ci(2,:)-coeffvalues(fitresult);
legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,vary{1}), sprintf('QUESP fit 1, fb=%.2e±%.1e, kb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );
title(func2str(modelstr));
% Label axes
xlabel( vary );
ylabel( ylab );
grid on

fb=0.0009;

end;


%% analytic formula comparison
% create a fast standard P
P.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
P.asym_fit      = 0;                    % which spectrum u want 2 fit? - cases : Z-spectrum(0) , Asym-Spectrum(1)
P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
P.n_cest_pool   = 1;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
P.MT_lineshape  = 'Gaussian';         % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
P.TR  = 3/1000;
P.linestomeasure=1;
P.flipangle= 14;
P.readout='bssfp';
P.dummies=1;  %% or better shots?
P.TTM_rep         = 0;    
P.FREQ          = 7*gamma_;       % frequency (=B0[T] * gamma)
P.B1            = 100;                  % B1 value in µT
P.Trec          = 0;                    % recover time in s
P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
P.Zi            = 1;                    % initial magnetisation (should be between -1 and +1)
P.shape         = 'block';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation
P.tp    = 3;                       % saturation time in s
P.n     = 1;                        % choose n=1 for cw saturation
P.shape = 'block';                  % choose 'block' for cw saturation
P.DC    = 1.0;                      % choose DC=1 for cw saturation
P.B1cwpe_quad   = -1;                     
P.td            = calc_td(P.tp,P.DC);       
P.c             = 1;                        

P.normalized=[]
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));

P.tissue    = 'PBS_PARA';
P.CESTagent = 'PARACEST';
P           = getSim(P);    

%% create 
Pso=P; %single offset
Pso.dwB=100;
Pso.analytic=1;
Pso.xZspec=[-Pso.dwB Pso.dwB];

%% QUEST

% original QUESP
QUEST1= @(fb,kb,R1,w1,x) fb.*kb./(R1+fb.*kb).*(w1).^2./(w1.^2+kb.^2).*(1- exp(-(R1+fb.*kb).*x) );

%fully revised            
QUEST2= @(fb,kb,R1,w1,x) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*(1- exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x) );

% general QUESP(t) with Zi
QUEST3= @(fb,kb,R1,w1,x) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x);

% inverse QUESP
QUESTi1= @(fb,kb,R1,w1,x) fb.*kb.*w1.^2./(w1.^2+kb.^2)./R1*x/x;

tp=2;
k=1:10000;
for B1=[ 1 5 25]
figure(2), plot(k,QUEST1(Pso.fB,k,P.R1A,B1*gamma_2pi,tp)); hold on
figure(2), plot(k,QUEST2(Pso.fB,k,P.R1A,B1*gamma_2pi,tp)); hold on
% figure(2), plot(k,QUEST3(Pso.fB,k,P.R1A,B1*gamma_2pi,tp)); hold on

 varyval=linspace(1,k(end),50);
MTR=fast_BM_MTR(Pso,{'kBA'},varyval,B1,tp);
figure(2), plot(varyval,MTR,'rx'); hold on

end;

%% QUESP

QUEST1= @(fb,kb,R1,w1,x) fb.*kb./(R1+fb.*kb).*(w1).^2./(w1.^2+kb.^2).*(1- exp(-(R1+fb.*kb).*x) );

%fully revised            
QUEST2= @(fb,kb,R1,w1,x) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*(1- exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x) );

% general QUESP(t) with Zi
QUEST3= @(fb,kb,R1,w1,x) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x);

% inverse QUESP
QUESTi1= @(fb,kb,R1,w1,x) fb.*kb.*w1.^2./(w1.^2+kb.^2)./R1*x/x;

tp=100;
k=1:10000;
for B1=[ 1 5 25 ]

figure(3), plot(k,QUEST1(Pso.fB,k,P.R1A,B1*gamma_2pi,tp)); hold on
figure(3), plot(k,QUEST2(Pso.fB,k,P.R1A,B1*gamma_2pi,tp)); hold on
% figure(3), plot(k,QUEST3(Pso.fB,k,P.R1A,B1*gamma_2pi,tp)); hold on

 varyval=linspace(1,k(end),50);
MTR=fast_BM_MTR(Pso,{'kBA'},varyval,B1,tp);
figure(3), plot(varyval,MTR,'rx'); hold on

end;

