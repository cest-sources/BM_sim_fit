function subtab = norm_run(subtab,offseti,excludi)
%  Ztab('ucl1','B1_run') = norm_run(Ztab('ucl1','B1_run'),1,0)

if nargin<3
    excludi=[];
end;

% normalize
II=1:numel(subtab{1,1});

for ii=II
    if ~isempty(subtab{1,1}{ii})
        subtab{1,1}{ii}=subtab{1,1}{ii}./subtab{1,1}{ii}(offseti);
        % exclude
        subtab{1,1}{ii}(excludi)=[];
    end;
end;

% exclude





