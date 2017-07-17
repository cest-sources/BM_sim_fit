function QUEST(Zlab,Zref,single_offset,P,ext_opts)
if iscell(P.vary)
    P.vary=P.vary{1};
end;

if strcmp(P.vary,'tp')==0
    error('It seems you varied %s, you need to vary tp for QUEST, or try QUESP',P.vary{1});
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

if nargin>4
opts.Lower  = ext_opts.Lower  ;
opts.StartPoint = ext_opts.StartPoint;
opts.Upper      = ext_opts.Upper ;

else
        warning('Internal Startvalue and boundaries used');
end;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A P.B1*gamma_2pi} );

% Plot fit with data.
figure(2);
subplot(4,1,ii);

h = plot( fitresult, xData, yData ); hold on;
ci = confint(fitresult);
ci = ci(2,:)-coeffvalues(fitresult);
legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,P.vary), sprintf('QUEST fit 1, fb=%.6f±%.9f, kb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );
title(func2str(modelstr));
% Label axes
xlabel( vary );
ylabel( 'MTR' );
grid on


end;