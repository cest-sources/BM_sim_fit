function val=tinv(percent,dgrfreed)

load stud_t


vecx=stud_t(:,1);
vec= stud_t(:,[0 0.5 0.75 0.8 0.9 0.95 0.98 0.99 0.998]==percent);

vecii=interp1(vecx,vec,1:10000,'linear','extrap');

 val=squeeze(vecii(dgrfreed));

    
