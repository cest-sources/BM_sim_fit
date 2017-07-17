function [quesp, quesp_inv] =  QUESP(Zlab,Zref,single_offset,P, ext_opts)
if iscell(P.vary)
    P.vary=P.vary{1};
end;


if strcmp(P.vary,'B1')==0
    error('It seems you varied %s, you need to vary B1 for QUEST',P.vary{1});
end;

for ii=1:2
    if ii==1
        MTR=Zref-Zlab;
        if isempty(P.normalized)
            modelstr= @(fb,kb,R1,x,w1) fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x);
        else
              modelstr= @(fb,kb,R1,x,w1) (fb.*kb.*w1.^2./(w1.^2+kb.^2)./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))-(P.Zi- R1./(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2))).*exp(-(R1+fb.*kb.*w1.^2./(w1.^2+kb.^2)).*x)+(P.Zi-1).*exp(-R1.*x))./(1+(P.Zi-1).*exp(-R1.*x));
        end;

        
    elseif ii==2
        
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
        
    end;

[xData, yData] = prepareCurveData( P.varyval*gamma_2pi, MTR );

% Set up fittype and options.
ft = fittype(modelstr, 'problem', {'R1','x'}, 'independent', {'w1'}, 'dependent', 'y' ); 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower      = [0.0000135 0];
opts.StartPoint = [0.000135 4000];
opts.Upper      = [0.0135 150000];

if nargin>4  % read in extern options
opts.Lower  = ext_opts.Lower  ;
opts.StartPoint = ext_opts.StartPoint;
opts.Upper      = ext_opts.Upper ;

else
        warning('Internal Startvalue and boundaries used');
end;
    

% Fit model to data.
if P.pulsed
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A (P.tp*P.n/P.DC)} );
else
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {P.R1A P.tp} );
end;

% Plot fit with data.
subplot(2,3,3+ii);
h = plot( fitresult, xData, yData,'o' ); hold on;
ci = confint(fitresult);
ci = ci(2,:)-coeffvalues(fitresult);


if ii==1 % MTRASYM
    quesp.kBA=[fitresult.kb ci(2)];
    quesp.fB=[fitresult.fb ci(1)];
    quesp.MTR=MTR;
    quesp.w_x=xData;
    quesp.fit=fitresult;
    
    legend( h, sprintf('MTR_{asym}(%.2f ppm) vs. %s',single_offset,P.vary), sprintf('QUESP fit, \nfb=%.2e±%.1e, \nkb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );

% Label axes
xlabel( P.vary );
ylabel( 'MTR_{asym}' );
grid on


end;

if ii==2 % MTRREX
    quesp_inv.kBA=[fitresult.kb ci(2)];
    quesp_inv.fB=[fitresult.fb ci(1)];
    quesp_inv.MTR=MTR;
    quesp_inv.w_x=xData;
    quesp_inv.fit=fitresult;
    ylabel( 'MTR_{Rex}' );
       legend( h, sprintf('MTR_{Rex}(%.2f ppm) vs. %s',single_offset,P.vary), sprintf('QUESP fit, \nfb=%.2e±%.1e, \nkb=%.2f±%.2f',fitresult.fb,ci(1),fitresult.kb,ci(2)), 'Location', 'NorthEast' );

% Label axes
xlabel( P.vary );

grid on
end;

end;