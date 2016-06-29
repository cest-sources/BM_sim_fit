function zspec=conv_ana_asym(startValue,xxx,P,dep_vars,vary,val)

if nargin<5
    multiple = 0;
else
    multiple = 1;
end;

for i=1:numel(dep_vars)
    P.(dep_vars{i}) = startValue(i);
end;

P.dwB = P.dwA+P.dwB;
P.kAB = P.kBA*P.fB;
P.dwC = P.dwA+P.dwC;
P.kAC = P.kCA*P.fC;
P.dwD = P.dwA+P.dwD;
P.kAD = P.kDA*P.fD;
P.dwE = P.dwA+P.dwE;
P.kAE = P.kEA*P.fE;
P.dwF = P.dwA+P.dwF;
P.kAF = P.kFA*P.fF;
P.dwG = P.dwA+P.dwG;
P.kAG = P.kGA*P.fG;

if multiple
    zspec=[];
    xZspec=[];
    for ii=1:numel(val(1,:))
        for jj = 1:numel(vary)
            P.(vary{jj}) = val(jj,ii)*P.c;
            P.td=calc_td(P.tp,P.DC);
        end
            Msim    = ANALYTIC_SIM(P);
            [x_t, y_t] = asym_PS(Msim.x, Msim.zspec);
            zspec   = [zspec ; y_t'];
            xZspec  = [xZspec ; x_t];
    end;
else
    Msim    = ANALYTIC_SIM(P);
    zspec   = Msim.zspec;
end;