function Rex_minilorentz = Rex_minilorentz(da,w1,di,fi,ki,r2i)
% calculates Rex_minilorentz for one pool (i)

ka=ki*fi;

REXMAX  = -ka.*w1^2./(w1.^2+ki.*(ki+r2i));
GAMMA   = 2*sqrt( (ki+r2i)./ki.*w1.^2 + (ki+r2i).^2);
Rex_minilorentz = REXMAX*(GAMMA./2).^2./((GAMMA./2).^2+di.^2);

end