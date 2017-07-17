function [Zrow, Mnorm] = norm_run(Zrow,runstr,offset)
%  Ztab('ucl1') = norm_run(Ztab('ucl1'),1,0)
if ~isempty(offset)
    
w=Zrow.offset{1};
subtab=Zrow{:,runstr};


offseti=find_nearest(w,offset);


% normalize
II=1:numel(subtab); % for all run entries

for ii=II
    if ~isempty(subtab{ii})
        Mnorm=subtab{ii}(offseti);
        subtab{ii}=subtab{ii}./Mnorm;

    end;
end;

Zrow{:,runstr}=subtab;
else
    warning('offset is empty, nothing was done!');
    Mnorm=1;
end;






