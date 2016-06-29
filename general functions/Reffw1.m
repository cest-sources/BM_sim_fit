function Reff = Reffw1(da,w1,r2a,r1a)
% calculates Reffw1 for pulsed gaussian solution
% last change: 2014/04/03 by PS

Reff = -(r1a.*da.^2./(da.^2+w1.^2) +(r2a).*w1.^2./(da.^2+w1.^2));

end

