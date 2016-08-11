%% this file allows to create simulated Z-spectrum data to check QUESP and QUEST formulas
%%it also allows to fit single offset QUESP and QUEST formulaes to the
%%(1) simulated data
%%(2) CESTtab structured data (go directly to 2 if you want to do so)

%% 1.1: create paraemter struct P for simulation or fit
% simulation parameters

P.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
P.asym_fit      = 0;                    % which spectrum u want 2 fit? - cases : Z-spectrum(0) , Asym-Spectrum(1)

P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
P.n_cest_pool   = 1;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
P.MT_lineshape  = 'Gaussian';         % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian

% sequence/scanner parameters

P.TR  = 3/1000;
P.linestomeasure=1;
P.flipangle= 14;
P.readout='bssfp';
P.dummies=1;  %% or better shots?
P.TTM_rep         = 0;    

P.offset        = 4;                    % offset in ppm

P.FREQ          = 7*gamma_;       % frequency (=B0[T] * gamma)
P.B1            = 10;                  % B1 value in µT
P.Trec          = 0;                    % recover time in s
P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
P.Zi            = 1;                    % initial magnetisation (should be between -1 and +1)

P.shape         = 'block';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation

if P.pulsed
    P.n     = 30;                       % number of saturation pulses
    P.tp    = 0.1;                      % saturation time per pulse in s
    P.DC    = 0.9;                      % duty cycle
else
    P.tp    = 20;                       % saturation time in s
    P.n     = 1;                        % choose n=1 for cw saturation
    P.shape = 'block';                  % choose 'block' for cw saturation
    P.DC    = 1.0;                      % choose DC=1 for cw saturation
end;

%XXXXXXXXXXXXXXX DONT CHANGE ! XXXXXXXXXXXXXXXX
P.B1cwpe_quad   = -1;                       %XX
P.td            = calc_td(P.tp,P.DC);       %XX
P.c             = 1;                        %XX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% 1.2: set t vary and varyval (manuel)
% which variables were varied in your measurements
%%1.2a:QUEST
      vary    = {'tp'}; % e.g. {'B1','Trec','DC'}
      clear varyval;
      varyval(1,:) = [10.0000    7.5000    5.0000    4.0000    3.0000    2.0000    1.5000    1.0000    0.7500    0.5000 ]
   
%% 1.2b:QUESP
      vary    = {'B1'};
      clear varyval;
      varyval(1,:) = [10 15 17.5 20 22.5 25 30 35 ]; % e.g. for B1: [0.5 0.75 1 2]
      
%% 1.3 create simulated data, adds noise and set depending variables and guess good starting values for the fit

P.normalized=[]
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));


% get start parameters for tissue CEST-agent
P.tissue    = 'PBS_PARA';
P.CESTagent = 'PARACEST';
P           = getSim(P);    

 
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
T.lowerB      = [   40         P.fB*0.9     50          1/20        0           0           0           0           ];
T.upperB      = [   60         P.fB*1.1        15000        10000         10000       10000       10000       10000       ];

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

% determine which case (--> parameters) is used
[dep_vars startValue lowerbounds upperbounds] = casedetermination(P,T);

% here the  simulation data is created
xxx=-80:1:80;
xZspec=repmat(xxx,numel(varyval));
P.xZspec=xxx;

% plot your measured and simulated data (uses starting point values)
if P.analytic
    if P.asym_fit
        Z = conv_ana_asym(startValue,xxx,P,dep_vars,vary,varyval);
    else
        Z = conv_ana(startValue,xxx,P,dep_vars,vary,varyval);
    end
else
    if P.asym_fit
        Z = conv_num_asym(startValue,xxx,P,dep_vars,vary,varyval);
    else
        Z = conv_num(startValue,xxx,P,dep_vars,vary,varyval);
    end
end
%here the simulation data is noised
Z = ricernd(Z',0.00001); % std value is 0.005
Z=Z';

% plot "1st-guess" with deviations from data points
figure(2001),

Zmat=reshape(Z,numel(xxx),numel(varyval));
plot(xxx,Zmat,'.-');
clear leg;

for ii=1:numel(startValue)
    lab{ii} = (sprintf('%s=%.4f',dep_vars{ii},startValue(ii)));
end; 

for ii=1:numel(varyval)
     leg{ii} = (sprintf('%s=%.2f',vary{1},varyval(ii)));
end; 
legend(leg);
xlabel(lab);


%% 2 QUEST and QUESP on experimental data  

%%2.1: pick QUEST and QUESP data from CESTstab structure
[xxx, xZspec, Z, varyval, vary, P]= plot_tab(Ztab,'goran1','B1_run');
% Zmat=reshape(Z,numel(xxx),numel(varyval));


%% 3 fitting of the Zmat data (wherever it comes from)

%% 3.1 RUN full BM OPTIMIZATION 

try close 2001 % close "1st-guess"
end
if P.analytic == 1 % use analytical solution
    if P.asym_fit % fit Asym-spectrum
        f = @(x,xdata) conv_ana_asym(x,xxx,P,dep_vars,vary,varyval);
    else % fit Z-spectrum
        f = @(x,xdata) conv_ana(x,xxx,P,dep_vars,vary,varyval);
    end
else % use numerical solution
    if P.asym_fit % fit Asym-spectrum
        f = @(x,xdata) conv_num_asym(x,xxx,P,dep_vars,vary,varyval);
    else % fit Z-spectrum
        f = @(x,xdata) conv_num(x,xxx,P,dep_vars,vary,varyval);
    end
end

%  ff= @(x,xdata) f(x,xdata)./(f(x,-80)+0.00001);

% fit-options
options = optimset('PlotFcns',{@optimplotfval},'TolFun',1e-8,'MaxIter',400,'Display','on');

if P.asym_fit
    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,x_asym_all,y_asym_all,lowerbounds,upperbounds,options);
else
    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,xZspec,Z,lowerbounds,upperbounds,options);
end

[ci, varb, corrb, varinf] = nlparci(x,RES,JACOBIAN,0.5);
ci(:,1) = ci(:,1)-x';
ci(:,2) = ci(:,2)-x';
figure, surf(imresize(abs(corrb), 1,'method','nearest'))
figure, imagesc(abs(corrb))
xlabel(dep_vars)
ylabel(dep_vars)


FIT=f(x,xxx);

for i=1:numel(dep_vars)
    fprintf('%s: real=%f  fit=%f+-%.5f  start=%f\n',dep_vars{i},P.(dep_vars{i}),x(i),ci(i,2), startValue(i));
end;

clear leg;
figure(2000)
if P.asym_fit
    plot(x_asym_all,y_asym_all,'x',x_asym_all,FIT);
else
    Zmat=reshape(Z,numel(xxx),numel(varyval));
    FITmat=reshape(FIT,numel(xxx),numel(varyval));
    plot(xxx,Zmat,'.',xxx,FITmat,xxx, Zmat-FITmat);
end

for ii=1:numel(startValue) 
    leg{ii} = (sprintf('%s=%.7f+-%.7f  (start: %.7f)',dep_vars{ii},x(ii),ci(ii,2),startValue(ii)));
end; 
xlabel(leg);
set(gca,'XDir','reverse');



%% 3.2 single offset preparation

Zmat=reshape(Z,numel(xxx),numel(varyval));  % 
clear leg
figure(42), plot(xxx,Zmat); title('check data again!');
for ii=1:numel(varyval)
     leg{ii} = (sprintf('%s=%.2f',vary{1},varyval(ii)));
end; 
legend(leg);

single_offset=50;
ind_lab=find_nearest(xxx,single_offset);
ind_ref=find_nearest(xxx,-single_offset);
Zlab=Zmat(ind_lab,:);
Zref=Zmat(ind_ref,:);

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
legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,vary{1}), sprintf('QUEST fit 1, fb=%.8f, kb=%.4f',fitresult.fb,fitresult.kb), 'Location', 'NorthEast' );
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

for ii=1:4
    if ii==1
MTR=Zref-Zlab;
% original QUESP

modelstr= @(fb,kb,R1,x,w1) fb.*kb./(R1+fb.*kb).*(w1).^2./(w1.^2+kb.^2).*(1- exp(-(R1+fb.*kb).*x) );
         
    elseif ii==2
MTR=Zref-Zlab;
%fully revised
            
modelstr= @(fb,kb,R1,x,w1) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*(1- exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x) );

    
    elseif ii==3
        % general QUESP(t)
        if isempty(P.normalized)
            modelstr= @(fb,kb,R1,x,w1) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x);
        else
              modelstr= @(fb,kb,R1,x,w1) (fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x))./(1+(P.Zi-1).*exp(-R1.*x));
        end;
        
        
    elseif ii==4
        
% inverse QUESP
MTR=1./Zlab-1./Zref;

modelstr= @(fb,kb,R1,x,w1) fb.*kb.*w1.^2./(w1.^2+kb.^2)./R1*x/x;

% .*(1-exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x));


    end;

[xData, yData] = prepareCurveData( varyval*gamma_2pi, MTR );

% Set up fittype and options.
ft = fittype(modelstr, 'problem', {'R1','x'}, 'independent', {'w1'}, 'dependent', 'y' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.0000135 0];
opts.StartPoint = [0.000135 4000];
opts.Upper = [0.00135 15000];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A P.tp} );

% Plot fit with data.
figure(1);
subplot(4,1,ii);
h = plot( fitresult, xData, yData ); hold on;
legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,vary{1}), sprintf('QUESP fit 1, fb=%.6f, kb=%.4f',fitresult.fb,fitresult.kb), 'Location', 'NorthEast' );
title(func2str(modelstr));
% Label axes
xlabel( vary );
ylabel( 'MTR' );
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
P.offset        = 4;                    % offset in ppm
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

tp=Pso.tp;
for B1=[ 1 10 100]
figure(2), plot(1:10000,QUEST1(Pso.fB,1:10000,P.R1A,B1*gamma_2pi,tp)); hold on
figure(2), plot(1:10000,QUEST2(Pso.fB,1:10000,P.R1A,B1*gamma_2pi,tp)); hold on
figure(2), plot(1:10000,QUEST3(Pso.fB,1:10000,P.R1A,B1*gamma_2pi,tp)); hold on

varyval=linspace(1,10000,50);
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
for B1=[ 1 10 100 ]

figure(3), plot(1:10000,QUEST1(Pso.fB,1:10000,P.R1A,B1*gamma_2pi,tp)); hold on
figure(3), plot(1:10000,QUEST2(Pso.fB,1:10000,P.R1A,B1*gamma_2pi,tp)); hold on
figure(3), plot(1:10000,QUEST3(Pso.fB,1:10000,P.R1A,B1*gamma_2pi,tp)); hold on

 varyval=linspace(1,10000,50);
MTR=fast_BM_MTR(Pso,{'kBA'},varyval,B1,tp);
figure(3), plot(varyval,MTR,'rx'); hold on

end;

