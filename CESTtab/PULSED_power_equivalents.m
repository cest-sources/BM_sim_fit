

%% analytic formula comparison
% create a fast standard P
P.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
P.asym_fit      = 0;                    % which spectrum u want 2 fit? - cases : Z-spectrum(0) , Asym-Spectrum(1)
P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
P.n_cest_pool   = 1;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
P.MT_lineshape  = 'Gaussian';         % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
P.TR  = 3/1000;
P.linestomeasure=1;
P.flipangle= 14;
P.readout='bssfp';
P.dummies=1;  %% or better shots?
P.TTM_rep         = 0;    
P.offset        = 4;                    % offset in ppm
P.FREQ          = 7*gamma_;       % frequency (=B0[T] * gamma)
P.B1            = 100;                  % B1 value in µT
P.Trec          = 0;                    % recover time in s
P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
P.Zi            = 1;                    % initial magnetisation (should be between -1 and +1)
P.shape         = 'block';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation
P.tp    = 3;                       % saturation time in s
P.n     = 1;                        % choose n=1 for cw saturation
P.shape = 'gauss_ucl';                  % choose 'block' for cw saturation
P.DC    = 1.0;                      % choose DC=1 for cw saturation
P.B1cwpe_quad   = -1;                     
P.td            = calc_td(P.tp,P.DC);       
P.c             = 1;                        

P.normalized=[]
warning(sprintf('P.normalized is at offset %.2f ppm',P.normalized));

P.tissue    = 'PBS_PARA';
P.CESTagent = 'PARACEST';
P           = getSim(P);    

P.xZspec=[-P.dwB P.dwB];

%%
clear II Rex_integrated Rex_cwpe Rex_cwae;

P.xZspec=[-P.dwB P.dwB];

II=1:20;
P.B1=50;
P.tp    = 0.1; 

P.kBA=5000;
P.kAB=P.fB*P.kBA;
P.DC=0.9091;

B1=linspace(0.1,10,numel(II));

for ii=II
ii*267.5153
 P.B1=B1(ii);
% P.kBA=ii;
% P.kAB=P.fB*P.kBA;
P.B1cwpe_quad   = -1;   
P.shape = 'gauss_ucl';                  % choose 'block' for cw saturation
[S Rex] = rho_s(P);
Rex_integrated(ii)=-S.Rex_hyper_integration(2)*P.DC;
Rex_integrated(ii)=Rex_mean(P,P.R2B,P.fB,P.kBA,P.B1)*P.DC;
                    
[w1 dphase theta, B1q]=RF_pulse(P,200);

P.B1cwpe_quad   = -1;   

P.B1= B1q; % power eqiv
P.shape = 'block';                  % choose 'block' for cw saturation
[S Rex] = rho_s(P);
Rex_cwpe(ii)=-S.Rex_Hyper_full(2);

P.B1= theta/gamma_2pi./((P.tp/P.DC)); % amplitide eqiv
P.shape = 'block';                  % choose 'block' for cw saturation
[S Rex] = rho_s(P);
Rex_cwae(ii)=-S.Rex_Hyper_full(2);

end;

figure(8),
c1= 1/(2*2.92)*sqrt(2*pi);
c22=(c1*sqrt(sqrt(2))).^2;
% c1=1;
% c22=1;
Rex_form=P.DC*c1.*P.fB.*P.kBA.*(B1*267.5153).^2./((B1*267.5153).^2+P.kBA*(P.kBA+P.R2B).*c22);
plot(II,Rex_integrated,II,Rex_cwpe,II,Rex_cwae,II,Rex_form); hold on;
legend({'Rex-integrated','Rex-cwpe','Rex-cwae','Rex+formfactors'});


