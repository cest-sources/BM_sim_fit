function [w_x, Z_x, w_xx, Z_xx, varyval, vary,P]= plot_tab(Ztab,row,col)
% row = row name 
% col= col name


w_xx =[];Z_xx =[];varyval= [];  % initialize this for multi param fit            

Zrow=Ztab({sprintf('%s',row)},:);
try
    P=Zrow.P{1};
catch
    P=Zrow.P;
end;



try % this is for the old notation with teh additional p row
    Zrowp=Ztab({sprintf('%s%s',row,'p')},:);
    for ii=1:numel(Zrowp.B1_run)
        if ~isempty(Zrowp.B1_run{ii})
            B1(ii)=Zrowp.B1_run{ii};
        end
    end
    
    for ii=1:numel(Zrowp.tsat_run)
        if ~isempty(Zrowp.tsat_run{ii})
            tsat(ii)=Zrowp.tsat_run{ii};
        end
    end
    
catch % this is for the new notation where var B1/tsat is stored in  P
    B1=P.varyval;
    tsat=P.varyval;
end;





if strcmp(col,'B1_run')
    
    for ii=1:numel(B1)  %%choose
        
        [w_x, wi]=sort(Zrow.offset{1});
        
        Z=Zrow.B1_run(1,:);
        if ~isempty(Z{ii})
            Zi=Z{ii}(wi);
            figure(42), plot(w_x,Zi,'DisplayName',sprintf('B1=%.2fµT',B1(ii))); hold on;    
            
            % this is for multi param fit    
            w_xx =[w_xx w_x(:)'];
            Z_xx =[Z_xx Zi(:)'];
            vary={'B1'};
            varyval= [varyval B1(ii)];
        end;
        clear Z Zi wi
               
    end;
     legend('show');
     title(Zrow.exp);
    
end

if strcmp(col,'tsat_run')
    
    for ii=1:numel(tsat)  %%choos
        
        [w_x, wi]=sort(Zrow.offset{1});
        
        Z=Zrow.tsat_run(1,:);
        if ~isempty(Z{ii})
            Zi=Z{ii}(wi);
            figure(42), plot(w_x,Zi,'DisplayName',sprintf('tsat=%.2f s',tsat(ii))); hold on;  
            
             % this is for multi param fit    
            w_xx =[w_xx w_x(:)'];
            Z_xx =[Z_xx Zi(:)'];
            vary={'tp'};
            varyval= [varyval tsat(ii)];
            
        end;
               
    end;
     legend('show');
    
end
if isrow(w_x)
w_x=w_x';
end;
Z_xx=Z_xx';

Z_x=reshape(Z_xx,numel(w_x),numel(varyval));
P.xZspec=w_x;

end

