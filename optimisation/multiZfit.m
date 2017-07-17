
function [FIT] =multiZfit(P,Sim,T,w_x,Z_x)
%MULTIZFIT simultaneous multiple Z-spectra Bloch-McConnell fit, e.g. Z(B1)
%
%   FIT = multiZfit(P,Sim,T,w_x,Z_x) yields the simultaneous Bloch-McConnell fit of multiple Z-spectra Z_x. 
% INPUT:
% P measurmenet parameters
% Sim simulation parameters
% T fit starting values and boundaries
% w_x   frequency axis vector in ppm
% Z_x   multi-B1 Z-spectrum stack Z_x(w_x,B1)
%
% OUTPUT:
% FIT.f=f(popt,w_x);            % fitted function
% FIT.popt=popt;                % optimal parameters
% FIT.pci=pci;                  % parameter confidence interval
% FIT.corrb=corrb;              % correlation matrix   
% FIT.R2=1-resnorm./mean(Z_xx); % goodness of fit parametr R^2
% FIT.T=T;                      % input startvalues and bopundaries T
% FIT.(dep_vars{ii});           % fitted parameters as field names as
% complex number, e.g. FIT.kBA= 320 +23i  means kBA= (320+- 23) Hz
%
%   Examples:
%   use BATCH files QUESP_and_QUEST_evaluations.m to set up 
%   (P,Sim,T,w_x,Z_x)
%
%   See also MULTIZPLOT, QUESP_and_QUEST_evaluations, norm_run,
%   exclude_run, b0correct_run

%   Version: see git repository www.cest-sources.org
%   Authors: Moritz Zaiss  - moritz.zaiss@tuebingen.mpg.de
%   Magentic Resonance Renter
%   Max-Planck institute of biological cybernetics Tübingen, Germany, 
%   http://www.kyb.tuebingen.mpg.de/nc/de/mitarbeiter/details/mzaiss.html
%   CEST sources - Copyright (C) 2017  Moritz Zaiss
%   **********************************
%   This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or(at your option) any later version.
%    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%    You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   **********************************


FIT.P=P;
FIT.Sim=Sim;

P = catstruct(Sim,P);

% this translates multi dim data to a single vector ready for teh
% optimization
w_xx=repmat(w_x,size(Z_x,2),1);
Z_xx= reshape(Z_x,size(Z_x,1)*size(Z_x,2),1);


try close 2001 % close "1st-guess"
end

% determine which case (--> parameters) is used
[dep_vars, startValue, lowerbounds, upperbounds] = casedetermination(P,T);

if P.analytic == 1 % use analytical solution
        f = @(x,xdata) conv_ana(x,w_x,P,dep_vars,P.vary,P.varyval);
else % use numerical solution
        f = @(x,xdata) conv_num(x,w_x,P,dep_vars,P.vary,P.varyval);
end

if isfield(T,'extOpt')
    options=T.extOpt;
else
    options = optimset('PlotFcns',{@optimplotfval},'TolFun',1e-8,'MaxIter',400,'Display','on');
end;
[popt resnorm RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,startValue,w_xx,Z_xx,lowerbounds,upperbounds,options);

try    
[ci, varb, corrb, varinf] = nlparci(popt,RES,JACOBIAN,0.95);
ci(:,1) = ci(:,1)-popt';
pci = ci(:,2)-popt';

catch 
    warning('nlparci could not be evaluated, ');
    pci=popt*NaN;
    corrb=NaN;
end;


FIT.f=f(popt,w_x);

FIT.popt=popt;
FIT.pci=pci;
FIT.corrb=corrb;
FIT.R2=1-resnorm./mean(Z_xx);
FIT.T=T;

for ii=1:numel(startValue)  
    FIT.(dep_vars{ii})(1)=popt(ii);
    FIT.(dep_vars{ii})(2)=pci(ii);
end; 

