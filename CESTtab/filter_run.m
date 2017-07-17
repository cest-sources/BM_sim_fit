function Zrow = filter_run(Zrow,runstr,ntimes)
%  Ztab('ucl1') = norm_run(Ztab('ucl1'),1,0)
    if nargin<3
        ntimes=1;
    end;
    
    subtab=Zrow{:,runstr};
    
    % normalize
    II=1:numel(subtab); % for all run entries
    
    for ii=II
        if ~isempty(subtab{ii})
            Z=subtab{ii};
                for jj=1:ntimes
                 Z = bfilt(Z);
                end;           
            subtab{ii}= Z;

        end;
    end;
    
    Zrow{:,runstr}=subtab;







