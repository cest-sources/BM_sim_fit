function [MTRasym,  MTRRex, Zlab, Zref]= fast_BM_MTR(P,vary,varyval,B1,tp)
%fast MTR calculation,  example: fast_BM_MTR(P,{'kBA'},[1:1000],1,1);
%calculates MTR for B1=1, tp=1 (overides) and a parameter struct P like below: 
% P.analytic      = 1;                    % Optimization type - cases: analytical(1), numerical(0)
% P.asym_fit      = 0;                    % which spectrum u want 2 fit? - cases : Z-spectrum(0) , Asym-Spectrum(1)
% P.MT            = 0;                    % 1 = with MT, 0 = no MT pool (MT is always pool C)
% P.n_cest_pool   = 1;                    % number of CEST/NOE pools (CEST pools: B,D,E,F,G)
% P.Rex_sol       = 'Hyper';              % cases: 'Hyper', 'Lorentz' , 'minilorentz'
% P.MT_lineshape  = 'Gaussian';         % MT lineshape - cases: SuperLorentzian, Gaussian, Lorentzian
% P.TR  = 3/1000;
% P.linestomeasure=1;
% P.flipangle= 14;
% P.readout='bssfp';
% P.dummies=1;  %% or better shots?
% P.TTM_rep         = 0;    
% P.offset        = 4;                    % offset in ppm
% P.FREQ          = 7*gamma_;       % frequency (=B0[T] * gamma)
% P.B1            = 100;                  % B1 value in µT
% P.Trec          = 0;                    % recover time in s
% P.spoilf        = 0;                    % spoilingfactor (0 = full spoiling, 1= no spoiling)
% P.Zi            = 1;                    % initial magnetisation (should be between -1 and +1)
% P.shape         = 'block';           % cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4
% P.pulsed        = 0;                    % 0 = cw saturation, 1 = pulsed saturation
% P.tp    = 3;                       % saturation time in s
% P.n     = 1;                        % choose n=1 for cw saturation
% P.shape = 'block';                  % choose 'block' for cw saturation
% P.DC    = 1.0;                      % choose DC=1 for cw saturation
% P.B1cwpe_quad   = -1;                     
% P.td            = calc_td(P.tp,P.DC);       
% P.c             = 1;                        
% % P.normalized=[]
% P.tissue    = 'PBS_PARA';
% P.CESTagent = 'PARACEST';
% P           = getSim(P);    

P.kAB = P.kBA*P.fB;
single_offset=P.dwB;

P.B1=B1;
P.tp=tp;

if P.analytic
    Z = conv_ana([],P.xZspec,P,[],vary,varyval);
else
    Z = conv_num([],P.xZspec,P,[],vary,varyval);
end

Zmat=reshape(Z,numel(P.xZspec),numel(varyval));

ind_lab=find_nearest(P.xZspec,single_offset);
ind_ref=find_nearest(P.xZspec,-single_offset);
Zlab=Zmat(ind_lab,:);
Zref=Zmat(ind_ref,:);

MTRasym=Zref-Zlab;

MTRRex=1./Zlab-1./Zref