function  [ci,varb,corrb,varinf]  = nlparci(beta,resid,J,quantil)
%NLPARCI Confidence intervals on parameters of nonlinear models.
%  %   and the variance and correlation matrix for the parameters b
%   VARINF are the variance inflation factors (if>10 ==> Multicollinearity)  = NLPARCI(BETA,RESID,J) returns the 95% confidence interval CI
%   and the variance and correlation matrix for the parameters beta
%   VARINF are the variance inflation factors (if>10 ==> Multicollinearity) 
%   on the nonlinear least squares parameter estimate BETA, given the 
%   residuals RESID, and the Jacobian matrix J at the solution.
%
%   The confidence interval calculation is valid for systems where 
%   the number of rows of J exceeds the length of BETA. 
%
%   NLPARCI uses the outputs of NLINFIT for its inputs.
%   Example:
%      [beta,resid,J]=nlinfit(input,output,'f',betainit);
%      [ci,varb,corrb,varinf] = nlparci(beta,resid,J);
%
%   See also NLINFIT.
%

%   Bradley Jones 1-28-94
%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.11 $  $Date: 2000/05/26 18:53:21 $
%   Modified by aj 08/01 to provide Var and Corr info also

%initialization
if nargin < 3
   error('Requires three inputs.');
end;

resid = resid(:);
[m,n] = size(J);
if m <= n
   error('The number of observations must exceed the number of parameters.');
end;

if length(beta) ~= n
   error('The length of beta must equal the number of columns in J.')
end

% approximation when a column is zero vector
temp = find(max(abs(J)) == 0);
if ~isempty(temp)
   J(temp,:) = J(temp,:) + sqrt(eps);
end;

%calculate covariance
[Q R] = qr(J,0);
Rinv = R\eye(size(R));
diag_info = sum((Rinv.*Rinv)')';

v = m-n;
rmse = sqrt(sum(resid.*resid)/v);

% calculate confidence interval
delta = sqrt(diag_info) .* rmse*tinv(quantil,v);
ci = [(beta(:) - delta) (beta(:) + delta)];

% Var and Corr Calculation for parameters --aj 08/01
xtxinv=Rinv*Rinv';% inv(X'X) see regstat
varb=xtxinv*rmse*rmse;
%calc correlation matrix corrb
varinf = diag(xtxinv); % variance inflation factors
corrb = xtxinv./sqrt(varinf*varinf');

%--end of nlparci.m---