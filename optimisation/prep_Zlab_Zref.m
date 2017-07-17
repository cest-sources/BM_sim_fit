function [Zlab, Zref]=prep_Zlab_Zref(P,single_offset,w_x,Z_x)
% prep_Zlab_Zref creates label and reference scan (opposite frequency) 
% from multi-Z-spectrum datahe
%
% [Zlab, Zref]=prep_Zlab_Zref(P,single_offset,w_x,Z_x)
% P is a parameter struct, this function reqires P.vary (string) and
% P.varyval (array of values), numel(P.varyval) must be equal to size(Z_x,2)
% single_offset is the offset to be evaluated, e.g 3.5
% w_x is the frequency axis numel(w_x)=size(Z_x,1)
% Z_x is teh multi-Z-spectrum stack with dimensions (numel(w_x),numel(P.varyval)

if iscell(P.vary)
    P.vary=P.vary{1};
end;

clear leg
plot(w_x,Z_x); hold on;
ax = gca; ax.ColorOrderIndex = 1;
for ii=1:numel(P.varyval)
     leg{ii} = (sprintf('%s=%.2f',P.vary,P.varyval(ii)));
end; 
legend(leg);

ind_lab=find_nearest(w_x,single_offset);
ind_ref=find_nearest(w_x,-single_offset);
Zlab=Z_x(ind_lab,:);
Zref=Z_x(ind_ref,:);

plot(w_x([ind_ref ind_lab]),[Zref; Zlab],'o');
title(sprintf('Single-offset evaluation\ncheck data, label and reference scan again! The chosen offsets %.2f ppm (lab) %.2f (ref)',w_x(ind_lab),w_x(ind_lab))); 
set(gca,'XDir','reverse');