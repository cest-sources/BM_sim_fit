function Zrow = exclude_run(Zrow,runstr,offset)
%  Ztab('ucl1') = norm_run(Ztab('ucl1'),1,0)
if ~isempty(offset)
    
w=Zrow.offset{1};
subtab=Zrow{:,runstr};

offseti=[];

    for ii=1:numel(offset)
    offseti=[offseti find_nearest(w,offset(ii))];
    end;



% normalize
II=1:numel(subtab); % for all run entries

for ii=II
    if ~isempty(subtab{ii})
       
        subtab{ii}(offseti)=[]; % delet from Z-spectrum
        
    end;
end;

Zrow.offset{1}(offseti)=[]; % delete from offsetlist

Zrow{:,runstr}=subtab;
else
    warning('offset is empty, nothing was done!');
end;






