function [Rex_MT R1obs] = Rex_MT(P,xZspec) 

w_ref   = 2*pi*P.FREQ;
w1      = P.B1*gamma_2pi;
da      = (xZspec-P.dwA)*w_ref;
dc      = (xZspec-P.dwC)*w_ref;
theta   = atan(w1./da);

if strcmp(P.MT_sol_type,'Rex_MT')
    
    Reff    = -(P.R1A*cos(theta).^2 + P.R2A*sin(theta).^2);
    r1a     = P.R1A+Reff;
    r2a     = P.R2A+Reff;
    r1c     = P.R1C+Reff;
    R2c     = P.R2C;
    rfmt    = RF_MT(1/R2c,w1,dc,P.MT_lineshape);
    
    Rex_MT =    -(((da.^2 + r2a.^2).*(P.kCA.*r1a + (P.kAC + r1a).*(r1c + rfmt)) + ...
                r2a.*(P.kCA + r1c + rfmt).*w1.^2)./(da.^2.*(P.kAC + P.kCA + r1a + r1c + rfmt) + ...
                r2a.*(P.kCA.*(2.*r1a + r2a) + r2a.*(r1c + rfmt) + ...
                P.kAC.*(2.*r1c + r2a + 2.*rfmt) + ...
                r1a.*(2.*r1c + r2a + 2.*rfmt)) + (P.kCA + r1c + r2a + rfmt).*w1.^2));
    
    R1obs =    0.5*( P.kAC + P.kCA + P.R1A+ P.R1C - sqrt(( P.kAC + P.kCA + P.R1A + P.R1C )^2 - 4*( P.kCA*P.R1A + P.kAC*P.R1C + P.R1A*P.R1C )));

end;


        

        
        
 

        
