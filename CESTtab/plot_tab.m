function [w_x, Z_x, w_xx, Z_xx, varyval, vary,P]= plot_tab(Ztab,row,col)
% row = row name 
% col= col name


w_xx =[];Z_xx =[];varyval= [];  % initialize this for multi param fit            
            
Zrow=Ztab({sprintf('%s%s',row,'p'),sprintf('%s',row)},:)
P=Zrow.P{2};
if strcmp(col,'B1_run')
    
    for ii=1:numel(Zrow.B1_run(2,:))  %%choose
        
        [w_x, wi]=sort(Zrow.offset{2});
        B1=Zrow.B1_run{1,ii};
        Z=Zrow.B1_run(2,:);
        if ~isempty(Z{ii})
            Zi=Z{ii}(wi);
            figure(42), plot(w_x,Zi,'DisplayName',sprintf('B1=%.2fµT',B1)); hold on;    
            
            % this is for multi param fit    
            w_xx =[w_xx w_x(:)'];
            Z_xx =[Z_xx Zi(:)']
            vary={'B1'};
            varyval= [varyval B1];
        end;
        clear B1 Z Zi wi B1
               
    end;
     legend('show');
     title(Zrow.exp{2});
    
end

if strcmp(col,'tsat_run')
    
    for ii=1:numel(Zrow.tsat_run(2,:))  %%choos
        
        [w_x, wi]=sort(Zrow.offset{2});
        tsat=Zrow.tsat_run{1,ii};
        Z=Zrow.tsat_run(2,:);
        if ~isempty(Z{ii})
            Zi=Z{ii}(wi);
            figure(42), plot(w_x,Zi,'DisplayName',sprintf('tsat=%.2f s',tsat)); hold on;  
            
             % this is for multi param fit    
            w_xx =[w_xx w_x(:)'];
            Z_xx =[Z_xx Zi(:)']
            vary={'tp'};
            varyval= [varyval tsat];
            
        end;
               
    end;
     legend('show');
    
end
if isrow(w_x)
w_x=w_x';
end;
Z_xx=Z_xx';

Z_x=reshape(Z_xx,numel(w_x),numel(varyval));

end

