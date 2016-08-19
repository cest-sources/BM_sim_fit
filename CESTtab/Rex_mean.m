function Rex_m=Rex_mean(P,r2b,fb,kb,B1in)

teile=200;

for ii=1:numel(B1in)
    P.B1=B1in(ii);
    [B1]=RF_pulse(P,teile,P.B1);
    Rex_m(ii) = -trapz(Rex_Hyper(P.dwB*2*pi*P.FREQ,gamma_2pi*B1,0,fb,kb,r2b))/teile ;
end;

if iscolumn(B1in)
    Rex_m=Rex_m';
end;

end


function Rex=Rex_Hyper(da,w1,di,fi,ki,r2i)
% calculate Rex HyperCEST solution for pool i
ka  =   ki*fi;
Rex =   -((ka.*ki.*w1.^2.*((-da+di).^2 + (r2i.*(da.^2 + (ka + ki).^2 + ki.*r2i + w1.^2))./ki))./...
        ((ka + ki).*(di.^2.*w1.^2 + ka.*r2i.*w1.^2) + ...
        (ka + ki).*((da.*di - ka.*r2i).^2 + (di.*ka + da.*(ki + r2i)).^2 + ...
        (ka + ki + r2i).^2.*w1.^2) + (ka + ki + r2i).*(da.^2.*w1.^2 + w1.^4)));
end
