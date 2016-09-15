function [ind value]=find_nearest(vec,val,howmany)
%[ind value]=find_nearest(vec,val,howmany)
% gives back the index and value in vec nearest to the given val

[vec VECi] =sort(vec);



rough_index=find((vec>=val),1);
if isempty(rough_index)
    rough_index=numel(vec);
end;

if rough_index>1
rough_index_pre=rough_index-1;
else
    rough_index_pre=rough_index;
end;

if abs(vec(rough_index)-val)<abs(vec(rough_index_pre)-val)
ind=rough_index;
else
    ind=rough_index_pre;
end;
    value=vec(ind);
    
if isempty(ind)
    
    try
        [C,I]=min(abs(vec-val*ones(numel(vec),1)));
    catch
        [C,I]=min(abs(vec-val*ones(1,numel(vec))));
    end
    ind=I;
    value=vec(ind)
end

vec=vec(VECi);
ind=VECi(ind);
