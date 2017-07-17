function [zspec] = conv_ana(startValue,xxx,P,dep_vars,vary,val)

if nargin<5
    multiple = 0;
else
    multiple = 1;
end;

for i=1:numel(dep_vars)
    P.(dep_vars{i}) = startValue(i);
end;


% calculate offset always from the water minimum (once)
P.dwB = P.dwA+P.dwB;
P.dwC = P.dwA+P.dwC;
P.dwD = P.dwA+P.dwD;
P.dwE = P.dwA+P.dwE;
P.dwF = P.dwA+P.dwF;
P.dwG = P.dwA+P.dwG;



if multiple
    zspec=[];
    xZspec=[];
    for ii=1:numel(val(1,:))
        
        if ~iscell(vary) % make sure both old an new def of vary work
            temp=vary;
            clear vary
            vary{1}=temp; clear temp;
           warning('please use new notation for vary: ist a struct now P.vary{1}=B1');
      
        end;
        
        for jj = 1:numel(vary) % thsi realizes multiple varying parameters
            P.(vary{jj}) = val(jj,ii)*P.c;
        end
        
% add function here for relations between variables
[P S] =get_BMrelations(P);
        
            [Msim Rex S]    = ANALYTIC_SIM(P);
            if isfield(P,'normalized')
                if ~isempty(P.normalized)
                    
                    Msim.zspec=Msim.zspec./Msim.zspec(find_nearest(xxx, P.normalized));
                    
                end;
            end;
                
            zspec   = [zspec ; abs(Msim.zspec)];
            xZspec  = [xZspec ; xxx];
%             assignin('base','Rex',Rex);  if Rex shall be written in the  workspace
%             assignin('base','S',S);
    end;
else
    [Msim Rex S]    = ANALYTIC_SIM(P);
     if isfield(P,'normalized')
        if P.normalized==1
            Msim.zspec=Msim.zspec./Msim.zspec(end);
            
        end;
    end;
    
    zspec   = Msim.zspec;
    assignin('base','Rex',Rex);
end;