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

P.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
P.asym_fit      = 0;                    % which spectrum u want 2 fit? - cases : Z-spectrum(0) , Asym-Spectrum(1)

P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
P.n_cest_pool   = 4;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
P.MT_lineshape  = 'Gaussian';         % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian

% sequence/scanner parameters

P.offset        = 4;                    % offset in ppm

P.FREQ          = 7*gamma_;       % frequency (=B0[T] * gamma)
P.B1            = 2;                  % B1 value in µT
P.Trec          = 3;                    % recover time in s
P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
P.Zi            = 1;                    % initial magnetisation (should be between -1 and +1)

P.shape         = 'block';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation

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

%% modify Ztab
% to exclude offsets use this code
Ztab=Ztab_fullglint;

rowname='ta2';                                 % define  the row you want to evaluate    
expname=Ztab{rowname,'exp'}{1}; 

Ztab(rowname,:) = norm_run(Ztab(rowname,:),'B1_run',P.normalized) 
Ztab(rowname,:) = exclude_run(Ztab(rowname,:),'B1_run',-5) 


%% READ DATA FROM Ztab
clear vary varyval
[w_x, Z_x, w_xx, Z_xx, varyval, vary, Ptab]= plot_tab(Ztab,rowname,'B1_run');

P.xZspec=w_x;
P.tp=Ptab.tsat;
P.Trec=Ptab.Trec;
P.Zi=cos(Ptab.readout_flipangle*pi/180);
P.B1=Ptab.B1;
P.FREQ=Ptab.B0*gamma_;

P.normalized=-4.6;

P.shape         = 'block';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation

if P.pulsed
    P.n     = Ptab.n;                       % number of saturation pulses
    P.tp    = Ptab.tp;                      % saturation time per pulse in s
    P.DC    = Ptab.DC;                      % duty cycle
end;

%%
P.shape         = 'SPINLOCK';
P.TTM_rep=0;
P.dummies=1;
P.analytic=0;
P.flipangle=10;
P.readout='bssfp';
P.linestomeasure=1;



%% set depending variables and guess good starting values
clear T
P.n_cest_pool   = 2;  
% get start parameters for tissue CEST-agent
P.tissue    = 'PBS_telaviv';
P.CESTagent = 'glucose_ta';
P           = getSim(P);        

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL A XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyA       = [   1           0           1       ];
T.dep_varsA   = {   'dwA',      'R1A',      'R2A'   };                  % {'dwA','R1A','R2A'};
T.startA      = [   P.dwA       P.R1A       P.R2A   ];                  % [P.dwA P.R1A P.R2A];
T.lowerA      = [   -0.1        1/4         1/3     ];
T.upperA      = [   +0.1        1/2         20     ];

[T.dep_varsA, T.startA, T.lowerA, T.upperA] = selectVars( T.varyA, T.dep_varsA, T.startA, T.lowerA, T.upperA );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL B XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyB       = [   0          0           1           1           0           0           0           0           ];
T.dep_varsB   = {   'dwB',      'fB',       'kBA',      'R2B',      'kBD',      'kBE',      'kBF',      'kBG'       };     
T.startB      = [   P.dwB       P.fB        P.kBA       P.R2B       P.kBD       P.kBE       P.kBF       P.kBG       ];
T.lowerB      = [   0.5        0.00001     5          1/20        0           0           0           0           ];
T.upperB      = [   2          0.01          100000        66         10000       10000       10000       10000       ];

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
% 
T.varyD       = [ 0         1           1           1           0           0           0           ];
T.dep_varsD   = {'dwD',     'fD',       'kDA',      'R2D',      'kDE',      'kDF',      'kDG'       };     
T.startD      = [P.dwD      P.fD        P.kDA       P.R2D       P.kDE       P.kDF       P.kDG       ];
T.lowerD      = [2.0        0.00001     1000         1/20        0           0           0           ];
T.upperD      = [3        0.01        100000         66         10000       10000       10000       ];
 
[T.dep_varsD, T.startD, T.lowerD, T.upperD] = selectVars( T.varyD, T.dep_varsD, T.startD, T.lowerD, T.upperD );
% 
% % XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% % XXXXXXXXXX POOL E XXXXXXXXXX
% % XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% 
% 
T.varyE       = [ 1         1           1           1           1           1           ];
T.dep_varsE   = {'dwE',     'fE',       'kEA',      'R2E',      'kEF',      'kEG',      };     
T.startE      = [P.dwE      P.fE        P.kEA       P.R2E       P.kEF       P.kEG       ];
T.lowerE      = [0.5        0.00001     10          1/20        0           0           ];
T.upperE      = [5       0.01        10000         20         10000       10000       ];

[T.dep_varsE, T.startE, T.lowerE, T.upperE] = selectVars( T.varyE, T.dep_varsE, T.startE, T.lowerE, T.upperE );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL F XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyF       = [ 1         1           1           1           1           ];
T.dep_varsF   = {'dwF',     'fF',       'kFA',      'R2F',      'kFG',      };     
T.startF      = [P.dwF      P.fF        P.kFA       P.R2F       P.kFG       ];
T.lowerF      = [0.3        0.00001     500         1/20        0           ];
T.upperF      = [5        0.01        10000       20         10000       ];

[T.dep_varsF, T.startF, T.lowerF, T.upperF] = selectVars( T.varyF, T.dep_varsF, T.startF, T.lowerF, T.upperF );

% XXXXXXXXXXXXXXXXXXXXXXXXXXXX
% XXXXXXXXXX POOL G XXXXXXXXXX
% XXXXXXXXXXXXXXXXXXXXXXXXXXXX

T.varyG       = [ 1         1           1           1           ];
T.dep_varsG   = {'dwG',     'fG',       'kGA',      'R2G',      };     
T.startG      = [P.dwF      P.fF        P.kFA       P.R2F       ];
T.lowerG      = [0          0.00001     10          1/20        ];
T.upperG      = [5          0.01        10000       200         ];

% [T.dep_varsG, T.startG, T.lowerG, T.upperG] = selectVars( T.varyG, T.dep_varsG, T.startG, T.lowerG, T.upperG );

% determine which case (--> parameters) is used
[dep_vars startValue lowerbounds upperbounds] = casedetermination(P,T);


% plot your measured and simulated data (uses starting point values)
if P.analytic
    if P.asym_fit
        FIT = conv_ana_asym(startValue,w_x,P,dep_vars,vary,varyval);
    else
        FIT = conv_ana(startValue,w_x,P,dep_vars,vary,varyval);
    end
else
    if P.asym_fit
        FIT = conv_num_asym(startValue,w_x,P,dep_vars,vary,varyval);
    else
        FIT = conv_num(startValue,w_x,P,dep_vars,vary,varyval);
    end
end

% plot "1st-guess" with deviations from data points
figure(2001),
if P.asym_fit % if u fit an Asym-spectrum
    plot(x_asym_all,y_asym_all,'x',x_asym_all,FIT,x_asym_all,y_asym_all-FIT); hold on;
    clear leg;
else % if u fit a Z-spectrum
    plot(w_xx,Z_xx,'x',w_xx,FIT,w_xx,Z_xx-FIT); hold on;
    clear leg;
    set(gca,'XDir','reverse');
end
for ii=1:numel(startValue)
    leg{ii} = (sprintf('%s=%.4f',dep_vars{ii},startValue(ii)));
end; 
xlabel(leg);



%% RUN OPTIMIZATION

try    close 2001 % close "1st-guess" 
end

if P.analytic == 1 % use analytical solution
    if P.asym_fit % fit Asym-spectrum
        f = @(x,xdata) conv_ana_asym(x,w_x,P,dep_vars,vary,varyval);
    else % fit Z-spectrum
        f = @(x,xdata) conv_ana(x,w_x,P,dep_vars,vary,varyval);
    end
else % use numerical solution
    if P.asym_fit % fit Asym-spectrum
        f = @(x,xdata) conv_num_asym(x,w_x,P,dep_vars,vary,varyval);
    else % fit Z-spectrum
        f = @(x,xdata) conv_num(x,w_x,P,dep_vars,vary,varyval);
    end
end

% fit-options
options = optimset('PlotFcns',{@optimplotfval},'TolFun',1e-8,'MaxIter',400,'Display','on');

if P.asym_fit
    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,x_asym_all,y_asym_all,lowerbounds,upperbounds,options);
else
    [x resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,w_xx,Z_xx,lowerbounds,upperbounds,options);
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

FIT=f(x,w_x);

for i=1:numel(dep_vars)
    fprintf('%s: real=%f  fit=%f+-%.5f  start=%f\n',dep_vars{i},P.(dep_vars{i}),x(i),ci(i,2), startValue(i));
end;

clear leg;
figure(2000)
if P.asym_fit
    plot(x_asym_all,y_asym_all,'x',x_asym_all,FIT);
else
    plot(w_xx,Z_xx,'x',w_xx,FIT,w_xx,Z_xx-FIT);
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




