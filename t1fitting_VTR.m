function [M0, T1map,c,S_fit] = t1fitting_VTR(S, TRs)
% =========================================================================
% T1 fitting using Variable TR.
% =========================================================================

%% Error checking
if ~nargin
 error('T1fitting:Inputs','There are no inputs.')
else if nargin <2
        error('T1fitting: Missing Input.');
    end
end
if size(S) ~= size (TRs)
    error('The dimensions of the signal vector and the TR vector should be equal.');
end

%% Default values
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective'...
    ,'Display','off','MaxFunEvals',100000,'MaxIter',100000,'TolFun',...
    1e-8,'TolX',1e-8);

lb = []; % lower bound
ub = []; % upper bound 
x0 = [400 1000 0]; % Initial values - hard coded

%% Equation
modelfun = @(p,x) (p(1).*(1-exp(-x./p(2)))+p(3));
% modelfun = @(p,x) (p(1).*(1-exp(-x./p(2))));

%% Curve fitting
p = lsqcurvefit(modelfun,x0,TRs,S,lb,ub,options);
M0 = p(1);
T1map = p(2);
c=p(3);
S_fit = modelfun(p, TRs);
    
%% Plot
     % times = linspace(TRs(1),TRs(end));
     % figure;
     % plot(TRs,S,'ko');
end

% =========================================================================
% 20191107 AK: 1st version.
% =========================================================================
