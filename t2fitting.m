function [M0_2, T2map, c_2, S_fit] = t2fitting(S, TEs)
% =========================================================================
%                               T2 fitting.
% =========================================================================

%% Error checking:
if ~nargin
 error('T2fitting:Inputs','There are no inputs.')
else if nargin <2
        error('T2fitting: Missing Input.');
    end
end
if size(S) ~= size (TEs)
    error('The dimensions of the signal vector and the TE vector should be equal.');
end

%% Default values
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective'...
    ,'Display','off','MaxFunEvals',100000,'MaxIter',100000,'TolFun',...
    1e-8,'TolX',1e-8);

lb = [];  % lower bound
ub = [];  % upper bound
%x0 = [400 10000 0];  % Initial values - hard coded
x0 = [3000 50 0];  % Initial values - hard coded

%% Equation
modelfun = @(p,x) ((p(1).*exp(-x./p(2)))+p(3));

%% Curve fitting
p = lsqcurvefit(modelfun,x0,TEs,S,lb,ub,options);
M0_2 = p(1);
T2map = p(2);
c_2 = p(3);
S_fit = modelfun(p,TEs);
    
%% Plot
%     times = linspace(TEs(1),TEs(end));
%     figure;
%     plot(TEs,S,'ko',times,modelfun(p,times),'b-');
end

% =========================================================================
%                           20191107 AK: 1st version.
% =========================================================================