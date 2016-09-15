function [S Rex] = rho_s(Sim)
% comments here
% last change: 2014/04/03 by PS

P = Sim;
xZspec = Sim.xZspec;

w_ref           = 2*pi*P.FREQ;
w1              = P.B1*gamma_2pi;
da              = (xZspec-P.dwA)*w_ref;
theta           = atan(w1./da);
S.theta         = theta;
S.Reff_sincos   = -P.R1A*cos(theta).^2 -(P.R2A)*sin(theta).^2;

if (strcmp(P.shape,'SPINLOCK')) || strcmp(P.shape,'block')
    
    switch P.Rex_sol    
        case 'Hyper' 
            if P.MT == 0;
                Rex = Rex_Hyper_full(P,xZspec);
                S.Rex_Hyper_full = Rex.Hyper_full;
                S.Rho_full = S.Reff_sincos + S.Rex_Hyper_full;
            else 
                Rex = Rex_Hyper_full(P,xZspec);
                S.Rex_Hyper_full = Rex.Hyper_full;
                [S.Rex_MT S.R1obs] = Rex_MT(P,xZspec);
                S.Rho_full = S.Reff_sincos + S.Rex_Hyper_full + S.Rex_MT;    
            end

        case 'Lorentz'
             if P.MT == 0;
                Rex = Rex_Lorentz_full(P,xZspec);
                S.Rex_Lorentz_full = Rex.Lorentz_full;
                S.Rho_full = S.Reff_sincos + S.Rex_Lorentz_full;
             else 
                Rex = Rex_Lorentz_full(P,xZspec);
                S.Rex_Lorentz_full = Rex.Lorentz_full;
                [S.Rex_MT S.R1obs] = Rex_MT(P,xZspec);
                S.Rho_full = S.Reff_sincos + S.Rex_Lorentz_full + S.Rex_MT;
            end

         case 'minilorentz'
             if P.MT == 0;
                Rex = Rex_minilorentz_full(P,xZspec);
                S.Rex_minilorentz_full = Rex.minilorentz_full;
                S.Rho_full = S.Reff_sincos + S.Rex_minilorentz_full;
             else 
                Rex = Rex_minilorentz_full(P,xZspec);
                S.Rex_minilorentz_full = Rex.minilorentz_full;
                [S.Rex_MT S.R1obs] = Rex_MT(P,xZspec);
                S.Rho_full = S.Reff_sincos + S.Rex_minilorentz_full + S.Rex_MT;
            end
    end
    
% all other pulses (e.g. seq_gauss, gauss, sinc...)    
else
    
    % preallocating for speed
    n = numel(da);
    teile=200;
    [B1]=RF_pulse(P,teile);
    
    Reff_mean   = ones(1,n);
    
    switch P.n_cest_pool
        
        case 0 
            da  = (xZspec-P.dwA)*w_ref;
            for i=1:n
                Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.R2A,P.R1A))/teile ;
            end
            
            S.R1rho_hyper_integration = Reff_mean;
            S.alpha_f_mean = 0;
            
            S.Rex_hyper_integration = 0;
                
        case 1
            db  = (xZspec-P.dwB)*w_ref;
            
            for i=1:n
                Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.R2A,P.R1A))/teile ;
                Rex_mean_b(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,db(i),P.fB,P.kBA,P.R2B))/teile ;
            end

            Rex_mean = Rex_mean_b;
            S.alpha_f_mean = -(Rex_mean_b/P.kBA);
            
            S.R1rho_hyper_integration=Reff_mean+Rex_mean; % R1p bar
            S.Rex_hyper_integration=Rex_mean;


        case 2
            db  = (xZspec-P.dwB)*w_ref;
            dd  = (xZspec-P.dwD)*w_ref;

            for i=1:n
                Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.R2A,P.R1A))/teile ;
                Rex_mean_b(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,db(i),P.fB,P.kBA,P.R2B))/teile ;
                Rex_mean_d(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,dd(i),P.fD,P.kDA,P.R2D))/teile ;
            end

            Rex_mean = Rex_mean_b + Rex_mean_d;
            S.alpha_f_mean = -(Rex_mean_b/P.kBA + Rex_mean_d/P.kDA);
            
            S.R1rho_hyper_integration=Reff_mean+Rex_mean;
            S.Rex_hyper_integration=Rex_mean;

        case 3
            db  = (xZspec-P.dwB)*w_ref;
            dd  = (xZspec-P.dwD)*w_ref;
            de  = (xZspec-P.dwE)*w_ref;

            for i=1:n
                Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.R2A,P.R1A))/teile ;
                Rex_mean_b(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,db(i),P.fB,P.kBA,P.R2B))/teile ;
                Rex_mean_d(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,dd(i),P.fD,P.kDA,P.R2D))/teile ;
                Rex_mean_e(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,de(i),P.fE,P.kEA,P.R2E))/teile ;
            end

            Rex_mean = Rex_mean_b + Rex_mean_d + Rex_mean_e;
            S.alpha_f_mean = -(Rex_mean_b/P.kBA + Rex_mean_d/P.kDA + Rex_mean_e/P.kEA);
            
            S.R1rho_hyper_integration=Reff_mean+Rex_mean;
            S.Rex_hyper_integration=Rex_mean;

        case 4
            db  = (xZspec-P.dwB)*w_ref;
            dd  = (xZspec-P.dwD)*w_ref;
            de  = (xZspec-P.dwE)*w_ref;
            df  = (xZspec-P.dwF)*w_ref;
            
            for i=1:n
                Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.R2A,P.R1A))/teile ;
                Rex_mean_b(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,db(i),P.fB,P.kBA,P.R2B))/teile ;
                Rex_mean_d(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,dd(i),P.fD,P.kDA,P.R2D))/teile ;
                Rex_mean_e(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,de(i),P.fE,P.kEA,P.R2E))/teile ;
                Rex_mean_f(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,df(i),P.fF,P.kFA,P.R2F))/teile ;
            end

            Rex_mean = Rex_mean_b + Rex_mean_d + Rex_mean_e + Rex_mean_f;
            S.alpha_f_mean = -(Rex_mean_b/P.kBA + Rex_mean_d/P.kDA + Rex_mean_e/P.kEA + Rex_mean_f/P.kFA );
            
            S.R1rho_hyper_integration=Reff_mean+Rex_mean;
            S.Rex_hyper_integration=Rex_mean; 
            
        case 5
            db  = (xZspec-P.dwB)*w_ref;
            dd  = (xZspec-P.dwD)*w_ref;
            de  = (xZspec-P.dwE)*w_ref;
            df  = (xZspec-P.dwF)*w_ref;
            dg  = (xZspec-P.dwG)*w_ref;
            
            for i=1:n
                Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.R2A,P.R1A))/teile ;
                Rex_mean_b(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,db(i),P.fB,P.kBA,P.R2B))/teile ;
                Rex_mean_d(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,dd(i),P.fD,P.kDA,P.R2D))/teile ;
                Rex_mean_e(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,de(i),P.fE,P.kEA,P.R2E))/teile ;
                Rex_mean_f(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,df(i),P.fF,P.kFA,P.R2F))/teile ;
                Rex_mean_g(i) = trapz(Rex_Hyper(da(i),gamma_2pi*B1,dg(i),P.fG,P.kGA,P.R2G))/teile ;
            end

            Rex_mean = Rex_mean_b + Rex_mean_d + Rex_mean_e + Rex_mean_f + Rex_mean_g;
            S.alpha_f_mean = -(Rex_mean_b/P.kBA + Rex_mean_d/P.kDA + Rex_mean_e/P.kEA + Rex_mean_f/P.kFA + Rex_mean_g/P.kGA);
            S.R1rho_hyper_integration=Reff_mean+Rex_mean;
            S.Rex_hyper_integration=Rex_mean;    
    end
    
    if P.MT == 1
        if (strcmp(P.shape,'SPINLOCK')) || strcmp(P.shape,'block')
            [S.Rex_MT S.R1obs] = Rex_MT(P,xZspec);  
            S.Rho_full = S.R1rho_hyper_integration + S.Rex_MT;
        else
            dc          = (xZspec-P.dwC)*w_ref;
            S.R1obs     = 0.5*( P.kAC + P.kCA + P.R1A+ P.R1C - sqrt(( P.kAC + P.kCA + P.R1A + P.R1C )^2 - 4*( P.kCA*P.R1A + P.kAC*P.R1C + P.R1A*P.R1C )));
            for i=1:n
                Rex_mean_MT(i) = trapz(Rex_MT2(da(i),gamma_2pi*B1,dc(i),P))/teile;
            end
            S.Rex_mean_MT   = Rex_mean_MT;
            S.Rho_full      = S.R1rho_hyper_integration + S.Rex_mean_MT;
            S.alpha_f_mean = S.alpha_f_mean - Rex_mean_MT/P.kCA;
            
        end
    else
        S.Rho_full = S.R1rho_hyper_integration;
    end

end

% minilorentz solution is always needed for pulsedCESL_approx
if P.n_cest_pool > 0
    aaa = Rex_minilorentz_full(P,xZspec);
    Rex.minilorentz_b = aaa.minilorentz_b;
else
    Rex = 0;
end

end
