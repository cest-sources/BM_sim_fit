function Rex = Rex_minilorentz_full(P,xZspec) 
% calculates Rex minilorentz solution for n pools
% pool A is the water pool
% pools B,D,E,F are CEST pools
% pool C is reserved for MT ... Rex_MT.m

w_ref   = 2*pi*P.FREQ;
w1      = P.B1*gamma_2pi;

switch P.n_cest_pool
    
    case 0
        Rex.minilorentz_full = 0;
    
    case 1
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        
        Rex.minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        
        Rex.minilorentz_full = Rex.minilorentz_b;
                
    case 2
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;

        Rex.minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        
        Rex.minilorentz_full = Rex.minilorentz_b + Rex.minilorentz_d;
        
    case 3
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        
        Rex.minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex.minilorentz_e = Rex_minilorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        
        Rex.minilorentz_full = Rex.minilorentz_b + Rex.minilorentz_d + Rex.minilorentz_e;
          
     case 4
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        df  = (xZspec-P.dwF)*w_ref;
        
        Rex.minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex.minilorentz_e = Rex_minilorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        Rex.minilorentz_f = Rex_minilorentz(da,w1,df,P.fF,P.kFA,P.R2F);
        
        Rex.minilorentz_full = Rex.minilorentz_b + Rex.minilorentz_d + Rex.minilorentz_e + Rex.minilorentz_f;
        
    case 5
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        df  = (xZspec-P.dwF)*w_ref;
        dg  = (xZspec-P.dwG)*w_ref;
        
        Rex.minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex.minilorentz_e = Rex_minilorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        Rex.minilorentz_f = Rex_minilorentz(da,w1,df,P.fF,P.kFA,P.R2F);
        Rex.minilorentz_g = Rex_minilorentz(da,w1,dg,P.fG,P.kGA,P.R2G);
        
        Rex.minilorentz_full = Rex.minilorentz_b + Rex.minilorentz_d + Rex.minilorentz_e + Rex.minilorentz_f + Rex.minilorentz_g;
        
end