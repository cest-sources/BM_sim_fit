function [w_x, Z_x, w_xx, Z_xx, varyval, vary,P]= plot_tab(Ztab,row,col,run_ind)
% row = row name 
% col= col name
% B1_ind  indices of B1 values that shoudl be chosen


subplot(2,1,1);

w_xx =[];Z_xx =[];varyval= [];  % initialize this for multi param fit            

if isstr(row)
Zrow=Ztab({sprintf('%s',row)},:);
else
   Zrow=Ztab(row,:);
end;
try
    P=Zrow.P{1};
catch
    P=Zrow.P;
end;

if nargin<4
    run_ind=1:size(P.varyval,2);
end;


try % this is for the old notation with teh additional p row
    Zrowp=Ztab({sprintf('%s%s',row,'p')},:);
    for ii=run_ind
        if ~isempty(Zrowp.B1_run{ii})
            B1(ii)=Zrowp.B1_run{ii};
        end
    end
    
    for ii=run_ind
        if ~isempty(Zrowp.tsat_run{ii})
            tsat(ii)=Zrowp.tsat_run{ii};
        end
    end
    
catch % this is for the new notation where var B1/tsat is stored in  P
    B1=P.varyval(run_ind);
    tsat=P.varyval(run_ind);
end;




    for ii=run_ind %%choose
        
        [w_x, wi]=sort(Zrow.offset{1});
        
        Z=Zrow.B1_run(1,:);
        if ~iscell(P.vary)
            vary{1}=P.vary;
        else
            vary=P.vary;
        end;
        
        
        if ~isempty(Z{ii})
            Zi=Z{ii}(wi);
             plot(w_x,Zi,'.','DisplayName',sprintf('%s=%.2f ',vary{1},P.varyval(:,ii))); hold on;    
             set(gca,'XDir','reverse');
            % this is for multi param fit    
            w_xx =[w_xx w_x(:)'];
            Z_xx =[Z_xx Zi(:)'];
            
        end;
        
        clear Z Zi wi
               
    end;
     legend('show');
     th=title(Zrow.exp);
     set(th,'interpreter','none');
    

if isrow(w_x)
w_x=w_x';
end;
Z_xx=Z_xx';

Z_x=reshape(Z_xx,numel(w_x),numel(run_ind));
P.xZspec=w_x;
ax = gca; ax.ColorOrderIndex = 1;
varyval= P.varyval(:,run_ind);



end



