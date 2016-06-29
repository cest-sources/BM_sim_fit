%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BATCH_fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
%            !!! MATLAB OPTIMIZATION TOOLBOX REQUIRED !!!
%            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
%            !!!!!!!!!! RUN SECTION BY SECTION !!!!!!!!!!
%            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% - choose analytical or numerical optimization, analytical is MUCH faster
%   NOTE: numeric optimization with pulsed saturation and 
%   shaped pulses (e.g. gauss, seq_gauss...) can take up to 48h !!!!
% - analytic MT modulation is not very accurate yet
% - if u run analytical optimization in pulsed case, please make sure the
%   correct pause modulation is used. Therefore check at the end of
%   mod_semi_pulsed.m which case is saved in M.zspec 
%   (e.g. if P.pulsed:  M.zspec = Z_pulsedCESL_approx)
%
% last change: 2014/04/09 by PS


%% clear all
clear all
clc

%% simulation parameters

P.analytic      = 0;                    % Optimization type - cases: analytical(1), numerical(0)
P.asym_fit      = 1;                    % which spectrum u want 2 fit? - cases : Z-spectrum(0) , Asym-Spectrum(1)

P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
P.n_cest_pool   = 2;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
P.MT_lineshape  = 'Gaussian';         % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian


%% sequence/scanner parameters

P.offset        = 4;                    % offset in ppm

P.FREQ          = 7*gamma_;       % frequency (=B0[T] * gamma)
P.B1            = 2;                  % B1 value in µT
P.Trec          = 4;                    % recover time in s
P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
P.Zi            = 1;                    % initial magnetisation (should be between -1 and +1)

P.shape         = 'SPINLOCK';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
P.pulsed        = 1;                    % 0 = cw saturation, 1 = pulsed saturation

if P.pulsed
    P.n     = 30;                       % number of saturation pulses
    P.tp    = 0.1;                      % saturation time per pulse in s
    P.DC    = 0.9;                      % duty cycle
else
    P.tp    = 12;                       % saturation time in s
    P.n     = 1;                        % choose n=1 for cw saturation
    P.shape = 'block';                  % choose 'block' for cw saturation
    P.DC    = 1.0;                      % choose DC=1 for cw saturation
end;

%XXXXXXXXXXXXXXX DONT CHANGE ! XXXXXXXXXXXXXXXX
P.B1cwpe_quad   = -1;                       %XX
P.td            = calc_td(P.tp,P.DC);       %XX
P.c             = 1;                        %XX
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% READ DATA FROM ZRUN STRUCT (CREATED WITH BATCH_ZRUN)

% reset data 
Z=[];
xZspec=[];
x_asym_all=[];
y_asym_all=[];
clear xxx zzz vary varyval

n_measure  = numel(ZRUN);
n_ROI      = numel(ZRUN{1});

% which ROI do u want 2 fit?
which_ROI = 1;

for ii=1:n_measure

    % read x- and Z-values for current measurement
    R       = ZRUN{ii}{which_ROI};
    xxx     = R.x([1:end]);
    zzz     = R.Zmean([1:end]);
    
    if isrow(xxx)
        xxx = xxx';
    end
    if isrow(zzz)
        zzz = zzz';
    end
    
    % save current values into one large multi-parameter spectrum
    xZspec=[xZspec ; xxx];
    Z=[Z ; zzz];
    
    if P.asym_fit
        [x_asym, y_asym] = asym_PS(xxx,zzz);
        x_asym = x_asym';
        y_asym = y_asym';
        x_asym_all = [x_asym_all ; x_asym];
        y_asym_all = [y_asym_all ; y_asym];
    end
end;

P.xZspec = xxx;

%% READ DATA FROM OPEN FIGURE

%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  !!! to make sure the script uses the correct figure, close all others !!
%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  !!!!!!! last plot will be read first, adapt your varyval vector !!!!!!!!
%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% reset data 
Z=[];
xZspec=[];
x_asym_all=[];
y_asym_all=[];
clear xxx zzz vary varyval

obj = get(gca,'Children');         

for i = 1:numel(obj);
    xxx=get(obj(i),'XDATA');               
    zzz=get(obj(i),'YDATA');
    
    if isrow(xxx)
        xxx = xxx';
    end
    if isrow(zzz)
        zzz = zzz';
    end
    
    % save current values into one large multi-parameter spectrum
    Z       = [Z ; zzz];
    xZspec  = [xZspec ; xxx];
    
    if P.asym_fit
        [x_asym, y_asym] = asym_PS(xxx,zzz);
        x_asym = x_asym';
        y_asym = y_asym';
        x_asym_all = [x_asym_all ; x_asym];
        y_asym_all = [y_asym_all ; y_asym];
    end
end
P.xZspec = xxx;
clear obj i zzz gX gZ

%% add rice/rician noise to simulated data
Z = ricernd(Z',0.001); % std value is 0.005
Z=Z';


%% set vary and varyval (manuel)
% which variables vary in your measurements
vary    = {'B1'}; % e.g. {'B1','Trec','DC'}

% set values (THE ORDER IS IMPORTANT)
varyval(1,:) = [6.4 3.2 1.6 0.8 0.4 0.2]; % e.g. for B1: [0.5 0.75 1 2]
% varyval(2,:) = [.8 .5 .8]; % e.g. for DC: [0.8 0.2 0.1 0.05]
% varyval(3,:) = [5 10 5]; % e.g. for n : [50 30 20 10]

%% set vary and varyval (semi automatic from maps)
% which variables vary in your measurements (and/or should be read from a
% map (e.g. R1A, R1B... from T1-map)
vary    = {'B1','R1A','R1B','R1C','R1D','R1E'}; % e.g. {'B1','Trec','DC'}

% which values does your variables take?
for ii=1:n_measure
    varyval(1,ii) = ZRUN{ii}{which_ROI}.B1; % e.g. B1: ZRUN{n_measure}{which_ROI}.B1
    varyval(2,ii) = 1/ZRUN{ii}{which_ROI}.T1; % e.g. R1A: 1/ZRUN{n_measure}{which_ROI}.T1
    varyval(3,ii) = 1/ZRUN{ii}{which_ROI}.T1; % e.g. R1B: 1/ZRUN{n_measure}{which_ROI}.T1
    varyval(4,ii) = 1/ZRUN{ii}{which_ROI}.T1; % e.g. R1B: 1/ZRUN{n_measure}{which_ROI}.T1
    varyval(5,ii) = 1/ZRUN{ii}{which_ROI}.T1; % e.g. R1B: 1/ZRUN{n_measure}{which_ROI}.T1
    varyval(6,ii) = 1/ZRUN{ii}{which_ROI}.T1; % e.g. R1B: 1/ZRUN{n_measure}{which_ROI}.T1
end

%% set depending variables and guess good starting values

% get start parameters for tissue CEST-agent
P.tissue    = 'simulation';
P.CESTagent = 'simulation';
P           = getSim(P);        

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL A XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyA       = [   1           1           1       ];
T.dep_varsA   = {   'dwA',      'R1A',      'R2A'   };                  % {'dwA','R1A','R2A'};
T.startA      = [   P.dwA       P.R1A       P.R2A   ];                  % [P.dwA P.R1A P.R2A];
T.lowerA      = [   -0.1        1/4         1/3     ];
T.upperA      = [   +0.1        1/2         200     ];

[T.dep_varsA, T.startA, T.lowerA, T.upperA] = selectVars( T.varyA, T.dep_varsA, T.startA, T.lowerA, T.upperA );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL B XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyB       = [   1           1           1           1           0           0           0           0           ];
T.dep_varsB   = {   'dwB',      'fB',       'kBA',      'R2B',      'kBD',      'kBE',      'kBF',      'kBG'       };     
T.startB      = [   P.dwB       P.fB/2        P.kBA/2       P.R2B       P.kBD       P.kBE       P.kBF       P.kBG       ];
T.lowerB      = [   1.0         0.00001     50          1/20        0           0           0           0           ];
T.upperB      = [   3.0         0.01        1000        200         10000       10000       10000       10000       ];

[T.dep_varsB, T.startB, T.lowerB, T.upperB] = selectVars( T.varyB, T.dep_varsB, T.startB, T.lowerB, T.upperB );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL C XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyC       = [ 0         0           0           0           0           ];
T.dep_varsC   = {'dwC',     'fC',       'kCA',      'R1C',      'R2C'       };     
T.startC      = [P.dwC      P.fC        P.kCA       P.R1C       P.R2C       ];
T.lowerC      = [0          0.0001      10          1/20        1/20        ];
T.upperC      = [9.9        0.1         10000       1/0.1       1/0.1       ];

[T.dep_varsC, T.startC, T.lowerC, T.upperC] = selectVars( T.varyC, T.dep_varsC, T.startC, T.lowerC, T.upperC );


% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL D XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX


T.varyD       = [ 1         1           1           1           0           0           0           ];
T.dep_varsD   = {'dwD',     'fD',       'kDA',      'R2D',      'kDE',      'kDF',      'kDG'       };     
T.startD      = [P.dwD      P.fD*2        P.kDA*2       P.R2D       P.kDE       P.kDF       P.kDG       ];
T.lowerD      = [-4.0        0.00001     10          1/20        0           0           0           ];
T.upperD      = [-1.0        0.01        500         200         10000       10000       10000       ];

[T.dep_varsD, T.startD, T.lowerD, T.upperD] = selectVars( T.varyD, T.dep_varsD, T.startD, T.lowerD, T.upperD );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL E XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX


T.varyE       = [ 1         1           1           1           1           1           ];
T.dep_varsE   = {'dwE',     'fE',       'kEA',      'R2E',      'kEF',      'kEG',      };     
T.startE      = [P.dwE      P.fE        P.kEA       P.R2E       P.kEF       P.kEG       ];
T.lowerE      = [2.7        0.00001     10          1/20        0           0           ];
T.upperE      = [3.0        0.01        500         200         10000       10000       ];

[T.dep_varsE, T.startE, T.lowerE, T.upperE] = selectVars( T.varyE, T.dep_varsE, T.startE, T.lowerE, T.upperE );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL F XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyF       = [ 1         1           1           1           1           ];
T.dep_varsF   = {'dwF',     'fF',       'kFA',      'R2F',      'kFG',      };     
T.startF      = [P.dwF      P.fF        P.kFA       P.R2F       P.kFG       ];
T.lowerF      = [0.3        0.00001     500         1/20        0           ];
T.upperF      = [1.0        0.01        10000       200         10000       ];

[T.dep_varsF, T.startF, T.lowerF, T.upperF] = selectVars( T.varyF, T.dep_varsF, T.startF, T.lowerF, T.upperF );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL G XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyG       = [ 1         1           1           1           ];
T.dep_varsG   = {'dwG',     'fG',       'kGA',      'R2G',      };     
T.startG      = [P.dwF      P.fF        P.kFA       P.R2F       ];
T.lowerG      = [0          0.00001     10          1/20        ];
T.upperG      = [3          0.01        10000       200         ];

[T.dep_varsG, T.startG, T.lowerG, T.upperG] = selectVars( T.varyG, T.dep_varsG, T.startG, T.lowerG, T.upperG );


% determine which case (--> parameters) is used
[dep_vars startValue lowerbounds upperbounds] = casedetermination(P,T);


% plot your measured and simulated data (uses starting point values)
if P.analytic
    if P.asym_fit
        FIT = conv_ana_asym(startValue,xxx,P,dep_vars,vary,varyval);
    else
        FIT = conv_ana(startValue,xxx,P,dep_vars,vary,varyval);
    end
else
    if P.asym_fit
        FIT = conv_num_asym(startValue,xxx,P,dep_vars,vary,varyval);
    else
        FIT = conv_num(startValue,xxx,P,dep_vars,vary,varyval);
    end
end

% plot "1st-guess" with deviations from data points
figure(2001),
if P.asym_fit % if u fit an Asym-spectrum
    plot(x_asym_all,y_asym_all,'x',x_asym_all,FIT,x_asym_all,y_asym_all-FIT); hold on;
    clear leg;
else % if u fit a Z-spectrum
    plot(xZspec,Z,'x',xZspec,FIT,xZspec,Z-FIT); hold on;
    clear leg;
end
for ii=1:numel(startValue)
    leg{ii} = (sprintf('%s=%.4f',dep_vars{ii},startValue(ii)));
end; 
xlabel(leg);



%% RUN OPTIMIZATION

close 2001 % close "1st-guess"

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

% fit-options
options = optimset('PlotFcns',{@optimplotfval},'TolFun',1e-6,'MaxIter',400,'Display','on');

if P.asym_fit
    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,x_asym_all,y_asym_all,lowerbounds,upperbounds,options);
else
    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,xZspec,Z,lowerbounds,upperbounds,options);
end

[ci, varb, corrb, varinf] = nlparci(x,resnorm,JACOBIAN,0.5);
ci(:,1) = ci(:,1)-x';
ci(:,2) = ci(:,2)-x';
figure, surf(imresize(abs(corrb), 1,'method','nearest'))
figure, imagesc(abs(corrb))
xlabel(dep_vars)
ylabel(dep_vars)

% @optimplotx plots the current point.
% @optimplotfunccount plots the function count.
% @optimplotfval plots the function value.
% @optimplotresnorm plots the norm of the residuals.
% @optimplotstepsize plots the step size.
% @optimplotfirstorderopt plots the first-order optimality measure.

FIT=f(x,xxx);

for i=1:numel(dep_vars)
    fprintf('%s: real=%f  fit=%f+-%.5f  start=%f\n',dep_vars{i},P.(dep_vars{i}),x(i),ci(i,2), startValue(i));
end;

clear leg;
figure(2000)
if P.asym_fit
    plot(x_asym_all,y_asym_all,'x',x_asym_all,FIT);
else
    plot(xZspec,Z,'x',xZspec,FIT);
end

for ii=1:numel(startValue) 
    leg{ii} = (sprintf('%s=%.5f+-%.5f  (start: %.5f)',dep_vars{ii},x(ii),ci(ii,2),startValue(ii)));
end; 
xlabel(leg);
set(gca,'XDir','reverse');


% cc = hsv(numel(ZRUN));
% figure (2001)
% hname = {'B_1 = 0.5µT','B_1 = 1.0µT','B_1 = 1.50µT','B_1 = 2.00µT','B_1 = 2.50µT'};
% for jj = 1:numel(ZRUN);
%     xxx = xxx;
%     yyy = FIT(((jj-1)*numel(xxx)+1 : jj*numel(xxx)));
%     h = plot(xxx,yyy,'color',cc(jj,:));
%     hvec(jj) = h;
%     hold on
%     plot(xxx,Z(((jj-1)*numel(xxx)+1 : jj*numel(xxx))),'x','color',cc(jj,:));   
% end
% xlabel('\Delta\omega [ppm]','FontSize',18);
% ylabel('M/M_0','FontSize',18);
% ylim([0.0 1.0]);
% xlim([-4 4]);
% grid on
% set(gca,'xdir','reverse');
% legend(hvec,hname,'Location','NorthWest','FontSize',14);   % legend
% hold off





%% SAVE OPTIM RESULTS
% get current folder
oldPath = pwd;
% choose folder
folder_name = uigetdir;
cd(folder_name);
% save workspace
save('matlab_optim_results.mat');
% create and switch into figure folder
mkdir('figures');
cd('.\figures');
% save figures
% saveas(figure(3),'figure3','fig');
% saveas(figure(4),'figure4','fig');
saveas(figure(2000),'figure2000','fig');
% switch back to old folder
cd(oldPath)




