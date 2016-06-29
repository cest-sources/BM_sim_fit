function Rex_Lorentz = Rex_Lorentz(da,w1,di,fi,ki,r2i)
% calculates Rex_Lorentz for one pool (i)

ka=ki*fi;

REXMAX  = -ka.*w1^2./(da.^2+w1.^2).*((da-di).^2 +(da.^2+w1.^2).*r2i./ki + r2i.*(ki+r2i));
GAMMA   = 2*sqrt( (ki+r2i)./ki.*w1.^2 + (ki+r2i).^2);
Rex_Lorentz = REXMAX./((GAMMA./2).^2+di.^2);

end