function alpha_f = alpha_f(P,xZspec)
% calculates alpha_f for approx pulsedCESL solution
% last change: 2014/04/02 by PS

w_ref=2*pi*P.FREQ;
w1 = P.B1*gamma_2pi;

switch P.n_cest_pool
    
    case 0
        alpha_f = 0;
    
    case 1
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        
        Rex_minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        
        alpha_f_b = -Rex_minilorentz_b./P.kBA;
        
        alpha_f = alpha_f_b;
        
    case 2
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        
        Rex_minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex_minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        
        alpha_f_b = -Rex_minilorentz_b./P.kBA;
        alpha_f_d = -Rex_minilorentz_d./P.kDA;
        
        alpha_f = alpha_f_b + alpha_f_d;
        
    case 3
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        
        Rex_minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex_minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex_minilorentz_e = Rex_minilorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        
        alpha_f_b = -Rex_minilorentz_b./P.kBA;
        alpha_f_d = -Rex_minilorentz_d./P.kDA;
        alpha_f_e = -Rex_minilorentz_e./P.kEA;
        
        alpha_f = alpha_f_b + alpha_f_d + alpha_f_e;
       
    case 4
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        df  = (xZspec-P.dwF)*w_ref;
        
        Rex_minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex_minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex_minilorentz_e = Rex_minilorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        Rex_minilorentz_f = Rex_minilorentz(da,w1,df,P.fF,P.kFA,P.R2F);
        
        alpha_f_b = -Rex_minilorentz_b./P.kBA;
        alpha_f_d = -Rex_minilorentz_d./P.kDA;
        alpha_f_e = -Rex_minilorentz_e./P.kEA;
        alpha_f_f = -Rex_minilorentz_f./P.kFA;
        
        alpha_f = alpha_f_b + alpha_f_d + alpha_f_e + alpha_f_f;
        
    case 5
        da  = (xZspec-P.dwA)*w_ref;
        db  = (xZspec-P.dwB)*w_ref;
        dd  = (xZspec-P.dwD)*w_ref;
        de  = (xZspec-P.dwE)*w_ref;
        df  = (xZspec-P.dwF)*w_ref;
        dg  = (xZspec-P.dwG)*w_ref;
        
        Rex_minilorentz_b = Rex_minilorentz(da,w1,db,P.fB,P.kBA,P.R2B);
        Rex_minilorentz_d = Rex_minilorentz(da,w1,dd,P.fD,P.kDA,P.R2D);
        Rex_minilorentz_e = Rex_minilorentz(da,w1,de,P.fE,P.kEA,P.R2E);
        Rex_minilorentz_f = Rex_minilorentz(da,w1,df,P.fF,P.kFA,P.R2F);
        Rex_minilorentz_g = Rex_minilorentz(da,w1,dg,P.fG,P.kGA,P.R2G);
        
        alpha_f_b = -Rex_minilorentz_b./P.kBA;
        alpha_f_d = -Rex_minilorentz_d./P.kDA;
        alpha_f_e = -Rex_minilorentz_e./P.kEA;
        alpha_f_f = -Rex_minilorentz_f./P.kFA;
        alpha_f_g = -Rex_minilorentz_g./P.kGA;
        
        alpha_f = alpha_f_b + alpha_f_d + alpha_f_e + alpha_f_f + alpha_f_g;
        
        
end

if P.MT == 1;
    da  = (xZspec-P.dwA)*w_ref;
    dc  = (xZspec-P.dwC)*w_ref;
    Rex_minilorentz_MT = Rex_minilorentz(da,w1,dc,P.fC,P.kCA,P.R2C);
    alpha_f_MT = -Rex_minilorentz_MT./P.kCA;
    alpha_f = alpha_f_MT + alpha_f;
end

end

