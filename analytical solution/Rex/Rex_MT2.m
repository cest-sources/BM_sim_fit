function [Rex_MT rfmt] = Rex_MT2(da,w1,dc,P)

    theta   = atan(w1./da);
    Reff    = -(P.R1A*cos(theta).^2 + P.R2A*sin(theta).^2);
    r1a     = P.R1A+Reff;
    r2a     = P.R2A+Reff;
    r1c     = P.R1C+Reff;
    R2C     = P.R2C;
    rfmt    = RF_MT(1/R2C,w1,dc,P.MT_lineshape);
    
    kCA     = P.kCA;
    kAC     = P.kCA*P.fC;

    Rex_MT =    -(((da.^2 + r2a.^2).*(kCA.*r1a + (kAC + r1a).*(r1c + rfmt)) + ...
                r2a.*(kCA + r1c + rfmt).*w1.^2)./(da.^2.*(kAC + kCA + r1a + r1c + rfmt) + ...
                r2a.*(kCA.*(2.*r1a + r2a) + r2a.*(r1c + rfmt) + ...
                kAC.*(2.*r1c + r2a + 2.*rfmt) + ...
                r1a.*(2.*r1c + r2a + 2.*rfmt)) + (kCA + r1c + r2a + rfmt).*w1.^2));
end


