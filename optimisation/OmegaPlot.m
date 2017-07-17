function [ome] = OmegaPlot(Zlab,Zref,single_offset,P,ext_opts,ind)
if iscell(P.vary)
    P.vary=P.vary{1};
end;

if nargin<6
    ind=1
end;
% mine
for ii=ind
    
if ii==1
MTR=1./Zlab-1./Zref;
ylab='(MTR_{Rex})^-1';


  modelstr= @(fb,kb,R1,xx) R1*(1./(fb.*kb) + kb./fb.*xx );
   


[xData, yData] = prepareCurveData( 1./(P.varyval*gamma_2pi).^2, 1./MTR );

elseif ii==2
   % original from dixons paper
    MTR=Zlab./(1-Zlab);
    ylab='MTR_{Dixon}';
    modelstr= @(fb,kb,R1,xx) R1*(1./(fb.*kb) + kb./fb.*xx )
[xData, yData] = prepareCurveData( 1./(P.varyval*gamma_2pi).^2, MTR );

else 
   % pulsed form factors
    MTR=1./Zlab-1./Zref;
    ylab='MTR_{Rex,pulsed}';
    c1= 1/(2*2.92)*sqrt(2*pi);
    c22=(c1*sqrt(sqrt(2))).^2;
    
    modelstr= @(fb,kb,R1,xx) R1./(P.DC.*c1)*(1./(fb.*kb) + c22.*kb./fb.*xx );
[xData, yData] = prepareCurveData( 1./(P.varyval*gamma_2pi).^2, 1./MTR );

end;


% Set up fittype and options.
ft = fittype(modelstr, 'problem', {'R1'}, 'independent', {'xx'}, 'dependent', 'y' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.0000135 0];
opts.StartPoint = [0.000135 2000];
opts.Upper = [0.00135 150000];

if nargin>4
opts.Lower  = ext_opts.Lower  ;
opts.StartPoint = ext_opts.StartPoint;
opts.Upper      = ext_opts.Upper ;

else
        warning('Internal Startvalue and boundaries used');
end;

opts.Weights=1./yData;
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A} );

% Plot fit with data.

h = plot( fitresult, xData, yData,'o' ); hold on;
ci = confint(fitresult);
ci = ci(2,:)-coeffvalues(fitresult);
legend( h, sprintf('MTR_{Rex}^-1(%.2f ppm) vs. %s',single_offset,P.vary), sprintf('Omega-plot fit, \nfb=%.2e±%.1e, \nkb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );
title(func2str(modelstr));
% Label axes
xlabel( P.vary );
ylabel( ylab );
grid on

fb=0.0009;

if ii==1
    ome.kBA=[fitresult.kb ci(2)];
    ome.fB=[fitresult.fb ci(1)];
    ome.MTR=MTR;
    ome.w_x=xData;
    ome.fit=fitresult;
end;

end;
