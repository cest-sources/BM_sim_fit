function [G2 G1]= superLorentzian_extrap_DG(t2b,w1,delta, cutoff)


for j=1:length(delta)
    if abs(delta(j)) >= cutoff
          % see Morrsion and Henkelman 1995.  units = s.  Seems weird to
          % me.  Need to multiply by w1^2 * pi to get saturation rate.
        integrand= @(u)sqrt(2/pi)*t2b./abs(3*u.^2-1) .* exp(-2*(delta(j).*t2b./(3*u.^2-1)).^2);
        G1(j)=w1.^2.*pi.*quad(integrand, 0,1);
        
        du=.0001;
        u=0.3:du:1;
        integrand2=sqrt(2/pi)*t2b./abs(3*u.^2-1) .* exp(-2*(delta(j)*t2b./abs(3*u.^2-1)).^2);
        G2(j)=w1.^2.*pi.*sum(integrand2)*du;
    else
        X = [-1.1*cutoff -cutoff cutoff 1.1*cutoff];
        Y = superLorentzian_extrap_DG(t2b,w1,X,cutoff);
        G2(j)=interp1(X,Y,delta(j),'spline');
    end
end
