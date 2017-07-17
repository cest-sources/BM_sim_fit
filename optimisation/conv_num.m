function zspec=conv_num(startValue,xxx,P,dep_vars,vary,val)

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

   if ~iscell(vary) % make sure both old an new def of vary work
            temp=vary;
            clear vary
            vary{1}=temp; clear temp;
           %warning('please use new notation for vary: ist a struct now P.vary{1}=B1');
      
        end;

if multiple
    zspec=[];
    xZspec=[];
    for ii=1:numel(val(1,:))
        
            
        for jj = 1:numel(vary)
            try
                P.(vary{jj}) = val(jj,ii)*P.c;
            catch
                P.(vary) = val(jj,ii)*P.c;
            end;
            P.td=calc_td(P.tp,P.DC);
        end
            Msim    = NUMERIC_SIM(P);
            
            if isfield(P,'normalized')
                if ~isempty(P.normalized)
                    
                    Msim.zspec=Msim.zspec./Msim.zspec(find_nearest(xxx, P.normalized));
                    
                end;
            end;
            
            zspec   = [zspec ; abs(Msim.zspec)];
            xZspec  = [xZspec ; xxx];
    end;
else
    Msim    = NUMERIC_SIM(P);
    
    if isfield(P,'normalized')
        if P.normalized==1
            Msim.zspec=Msim.zspec./Msim.zspec(end);
            
        end;
    end;
    zspec   = Msim.zspec;
end;

