function Rex = Rex_Lorentz_full(P,xZspec) 
% calculates Rex Lorentz solution for n pools
% pool A is the water pool
% pools B,D,E,F are CEST pools
% pool C is reserved for MT ... Rex_MT.m

w_ref   = 2*pi*P.FREQ;
w1      = P.B1*gamma_2pi;

switch P.n_cest_pool
    
    case 0
        Rex.Lorentz_full = 0;

    case 1
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        
        Rex.Lorentz_b = Rex_Lorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        
        Rex.Lorentz_full = Rex.Lorentz_b;
                
    case 2
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;

        Rex.Lorentz_b = Rex_Lorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.Lorentz_d = Rex_Lorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        
        Rex.Lorentz_full = Rex.Lorentz_b + Rex.Lorentz_d;
        
    case 3
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        
        Rex.Lorentz_b = Rex_Lorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.Lorentz_d = Rex_Lorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex.Lorentz_e = Rex_Lorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        
        Rex.Lorentz_full = Rex.Lorentz_b + Rex.Lorentz_d + Rex.Lorentz_e;
          
     case 4
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        df  = (xZspec-P.dwF)*w_ref;
        
        Rex.Lorentz_b = Rex_Lorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.Lorentz_d = Rex_Lorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex.Lorentz_e = Rex_Lorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        Rex.Lorentz_f = Rex_Lorentz(da,w1,df,P.fF,P.kFA,P.R2F);
        
        Rex.Lorentz_full = Rex.Lorentz_b + Rex.Lorentz_d + Rex.Lorentz_e + Rex.Lorentz_f;
        
    case 5
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        df  = (xZspec-P.dwF)*w_ref;
        dg  = (xZspec-P.dwG)*w_ref;
        
        Rex.Lorentz_b = Rex_Lorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex.Lorentz_d = Rex_Lorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex.Lorentz_e = Rex_Lorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        Rex.Lorentz_f = Rex_Lorentz(da,w1,df,P.fF,P.kFA,P.R2F);
        Rex.Lorentz_g = Rex_Lorentz(da,w1,dg,P.fG,P.kGA,P.R2G);
        
        Rex.Lorentz_full = Rex.Lorentz_b + Rex.Lorentz_d + Rex.Lorentz_e + Rex.Lorentz_f + Rex.Lorentz_g;
        
end