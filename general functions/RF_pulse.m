function [w1 dphase theta B1q]=RF_pulse(P,teile,B1overide)
% comments here
% last change: 2015/06/11 by PS
if nargin>2
    P.B1=B1overide;
end;

tpulse=0:P.tp/teile:P.tp;
dphase=0;

if any(strcmp(P.shape,'SPINLOCK'))
    % lineare erhöhung der pulse bei td ungl. 0
    % w1 =tpulse*0+ (P.n*P.tp+(P.n-1)*P.td)/(P.n*P.tp)* P.B1;
    w1 =tpulse*0+ P.B1;
    %B1t = P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);

    if P.B1cwpe_quad==1
        %B1t = P.B1; % quadratische erhöhung der pulse bei interpulsdelay ungl. 0
        w1 =tpulse*0+sqrt((P.n*P.tp+(P.n-1)*P.td)/(P.n*P.tp))* P.B1;

    elseif P.B1cwpe_quad==0
        w1=w1./P.DC; 
    end;

    
    
elseif any(strcmp(P.shape,'block'))
    % lineare erhöhung der pulse bei td ungl. 0
    % w1 =tpulse*0+ (P.n*P.tp+(P.n-1)*P.td)/(P.n*P.tp)* P.B1;
    w1 =tpulse*0+ 1/P.DC* P.B1;
    %B1t = P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);

    if P.B1cwpe_quad==1
        %B1t = P.B1; % quadratische erhöhung der pulse bei interpulsdelay ungl. 0
        w1 =tpulse*0+sqrt((P.n*P.tp+(P.n-1)*P.td)/(P.n*P.tp))* P.B1;

    elseif P.B1cwpe_quad==-1
        w1=w1.*P.DC; 
    end;

    
    
elseif any(strcmp(P.shape,'block_trap'))
    % lineare erhöhung der pulse bei td ungl. 0
    % w1 =tpulse*0+ (P.n*P.tp+(P.n-1)*P.td)/(P.n*P.tp)* P.B1;
    w1 =tpulse*0+ 1/P.DC* P.B1;
    helpx=0.1:0.1:1;
    helpxx=1:-0.1:0.1;
    w1(1:10)=helpx*1/P.DC* P.B1
    w1(end-9:end)=helpxx*1/P.DC* P.B1
    %B1t = P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);

    if P.B1cwpe_quad==1
        %B1t = P.B1; % quadratische erhöhung der pulse bei interpulsdelay ungl. 0
        w1 =tpulse*0+sqrt((P.n*P.tp+(P.n-1)*P.td)/(P.n*P.tp))* P.B1;

    elseif P.B1cwpe_quad==-1
        w1=w1.*P.DC; 
    end;

    
    
elseif any(strcmp(P.shape,'seq_gauss'))
    nsig=2.229;
    epsilon=0.00001;
    %b1 linear
    gauss=@(t,t0,w1max,sig) 1./sqrt(2*pi*sig^2) *w1max.*exp(-(t-t0).^2./(2.*sig.^2));   
    w1=gauss(tpulse,P.tp/2,(P.tp+P.td.*( (P.n-1)./P.n ))*P.B1,P.tp./(2*nsig));
    w1=w1-w1(1)+epsilon;
    %b1 linear
    %norm=trapz(tpulse,abs(w1))./((P.tp+P.td*( (P.n-1)/P.n )));
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1./norm*P.B1;

    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp/P.DC)));
        w1=w1./norm*P.B1;
    elseif P.B1cwpe_quad==-1
        w1=w1.*P.DC; 
    end;
    
elseif any(strcmp(P.shape,'seq_siemens'))  %% got from scanner VD (20ms, 0.6µT)
   
    gauss=@(x,w1max) w1max*(-25.88*(x/P.tp).^6 + 76.88*(x/P.tp).^5  -67.47*(x/P.tp).^4 + 8.011*(x/P.tp).^3 + 8.034*(x/P.tp).^2 +  0.4235*(x/P.tp)  -0.0002965);
    w1=gauss(tpulse,1);
    
    %b1 linear
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1./norm*P.B1;

    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp/P.DC)));
        w1=w1./norm*P.B1;
    elseif P.B1cwpe_quad==-1
        w1=w1.*P.DC; 
    end;
    
    
     
  
    
    
elseif any(strcmp(P.shape,'gauss'))    
    %b1 linear
    gauss=@(t,t0,w1max,sig) 1./sqrt(2*pi*sig^2) *w1max.*exp(-(t-t0).^2./(2.*sig.^2));
    %b1 linear
    w1=gauss(tpulse,P.tp/2,(P.tp/P.DC)*P.B1,P.tp./(2*3));
    theta=trapz(tpulse,w1*gamma_2pi);
    
    if P.B1cwpe_quad==1
        %b1 quadratisch
        gauss=@(t,t0,w1max,sig) 1/sqrt(sqrt(pi)*sig)*w1max.*exp(-(t-t0).^2./(2.*sig^2));
        %b1 quadratisch
        w1=gauss(tpulse,P.tp/2,sqrt((P.tp/P.DC) )*P.B1,P.tp/(2*3));
    end;
    
elseif any(strcmp(P.shape,'gauss_ucl'))
    nsig=2.92;
    epsilon=0.00001;
    %b1 linear
    gauss=@(t,t0,w1max,sig) 1./sqrt(2*pi*sig^2) *w1max.*exp(-(t-t0).^2./(2.*sig.^2)); 
    
    w1=gauss(tpulse,P.tp/2,P.B1,P.tp./(2*nsig));
    
    %b1 linear
    %norm=trapz(tpulse,abs(w1))./((P.tp+P.td*( (P.n-1)/P.n )));
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1./norm*P.B1;

    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp/P.DC)));
        w1=w1./norm*P.B1;
    elseif P.B1cwpe_quad==-1
        w1=w1.*P.DC; 
    end;
    

 
    
elseif any(strcmp(P.shape,'sech'))   
    f_sech=@(t,t0,w1max,sig) w1max.*sech((t-t0).*sig);
    %B1 linear
    w1=f_sech(tpulse,P.tp/2,(P.tp+P.td*( (P.n-1)/P.n ))*P.B1,P.nsig);
    %Normiere auf P.B1
    norm=trapz(tpulse,abs(w1))/((P.tp+P.td*( (P.n-1)/P.n )));
    w1=w1/norm*P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);
    %dphase
    f_phase = @(t,t0,mu,sig)  mu*(-tanh((t-t0).*sig));
    dphase= f_phase(tpulse,P.tp/2,P.mu,P.beta);
    plot(tpulse,w1./max(w1),tpulse,dphase)
  
elseif any(strcmp(P.shape,'AdiaSL'))
    if P.AdiaAmpFactor == 1
        t0 = tpulse(end);
    elseif P.AdiaAmpFactor == 0
        t0 = 0;
    end
    f_Abs = @ (t,t0,w1max,mu,bw) w1max./cosh( (bw*pi()/mu).*(t-t0));
    w1 = f_Abs(tpulse,t0,P.AdiaB1,P.AdiaMu,P.AdiaBw);
    
    %Phase
    f_Pha = @ (t,t0,mu,bw,FreqFactor) FreqFactor*bw*pi().*tanh((bw*pi()/mu).*(t-t0));
    dfreq = f_Pha(tpulse,t0,P.AdiaMu,P.AdiaBw,P.AdiaFreqFactor);
    dphase = (dfreq)/(2*pi()*P.FREQ); %convert from rad to ppm
    
elseif any(strcmp(P.shape,'AdiaSinCos'))
    
    if P.AdiaAmpFactor == 1
        t0 = 0;
    elseif P.AdiaAmpFactor == 0
        t0 = tpulse(end);
    end
    
    f_Abs = @ (t,t0,w1max,mu,bw) w1max.*sin( bw*pi()/mu.*(t+t0));
    w1 = f_Abs(tpulse,t0,P.AdiaB1,P.AdiaMu,P.AdiaBw);
    
    %Phase
    f_Pha = @ (t,t0,mu,bw,FreqFactor) FreqFactor*bw*pi().*cos( bw*pi()/mu.*(t+t0));
    dfreq = f_Pha(tpulse,t0,P.AdiaMu,P.AdiaBw,P.AdiaFreqFactor);
    dphase = (dfreq)/(2*pi()*P.FREQ); %convert from rad to ppm
    
    figure
    subplot(1,3,1)
    plot(tpulse,w1)
    title('B1-modulation')
    subplot(1,3,2)
    plot(tpulse,dfreq)
    title('Freq-Modulation')
    subplot(1,3,3)
    plot(tpulse,dphase)
    title('Phase-Modulation')
    

elseif any(strcmp(P.shape,'AdiaInversion'))
%     if P.AdiaAmpFactor == 1
        t0 = tpulse(end);
%     elseif P.AdiaAmpFactor == 0
%         t0 = 0;
%     end
    f_Abs = @ (t,t0,w1max,mu,bw) w1max./cosh( (bw*pi()/mu).*(2.*t-t0));
    w1 = f_Abs(tpulse,t0,P.AdiaB1,P.AdiaMu,P.AdiaBw);
    
    %Phase
    f_Pha = @ (t,t0,mu,bw,FreqFactor) FreqFactor*bw*pi().*tanh((bw*pi()/mu).*(2.*t-t0));
    dfreq = f_Pha(tpulse,t0,P.AdiaMu,P.AdiaBw,P.AdiaFreqFactor);
    dphase = (dfreq)/(2*pi()*P.FREQ); %convert from rad to ppm
    
    figure
    subplot(1,3,1)
    plot(tpulse,w1)
    title('B1-modulation')
    subplot(1,3,2)
    plot(tpulse,dfreq)
    title('Freq-Modulation')
    subplot(1,3,3)
    plot(tpulse,dphase)
    title('Phase-Modulation')
      
elseif any(strcmp(P.shape,'sinc_1'))
    epsilon=0.00001;
    wiggles = 1;  
    sinc=@(t,t0,w1max,sig) epsilon+w1max.*(epsilon+sin((t-t0)*pi./sig))./(epsilon+(t-t0)*pi./sig);
    w1=sinc(tpulse,P.tp/2,(P.tp+P.td*( (P.n-1)/P.n ))*P.B1,P.tp/(2*wiggles));
    % norm=trapz(tpulse,abs(w1))./((P.tp+P.td*( (P.n-1)/P.n )));
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1/norm*P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);
    
    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp+P.td*( (P.n-1)/P.n ))));
        w1=w1./norm*P.B1;
    end

    
    
elseif any(strcmp(P.shape,'sinc_2'))
    epsilon=0.00001;
    wiggles = 2;    
    sinc=@(t,t0,w1max,sig) epsilon+w1max.*(epsilon+sin((t-t0)*pi./sig))./(epsilon+(t-t0)*pi./sig);
    w1=sinc(tpulse,P.tp/2,(P.tp+P.td*( (P.n-1)/P.n ))*P.B1,P.tp/(2*wiggles));
    % norm=trapz(tpulse,abs(w1))./((P.tp+P.td*( (P.n-1)/P.n )));
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1/norm*P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);
    
    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp+P.td*( (P.n-1)/P.n ))));
        w1=w1./norm*P.B1;
    end;

elseif any(strcmp(P.shape,'sinc_3'))
    epsilon=0.00001;
    wiggles = 3;   
    sinc=@(t,t0,w1max,sig) epsilon+w1max.*(epsilon+sin((t-t0)*pi./sig))./(epsilon+(t-t0)*pi./sig);
    w1=sinc(tpulse,P.tp/2,(P.tp+P.td*( (P.n-1)/P.n ))*P.B1,P.tp/(2*wiggles));
    % norm=trapz(tpulse,abs(w1))./((P.tp+P.td*( (P.n-1)/P.n )));
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1/norm*P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);
    
    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp+P.td*( (P.n-1)/P.n ))));
        w1=w1./norm*P.B1;
    end;
    
    
    
elseif any(strcmp(P.shape,'sinc_4'))
    epsilon=0.00001;
    wiggles = 6;    
    sinc=@(t,t0,w1max,sig) epsilon+w1max.*(epsilon+sin((t-t0)*pi./sig))./(epsilon+(t-t0)*pi./sig);
    w1=sinc(tpulse,P.tp/2,(P.tp+P.td*( (P.n-1)/P.n ))*P.B1,P.tp/(2*wiggles));
    % norm=trapz(tpulse,abs(w1))./((P.tp+P.td*( (P.n-1)/P.n )));
    norm=trapz(tpulse,abs(w1))./(P.tp/P.DC);
    w1=w1/norm*P.B1;
    theta=trapz(tpulse,w1*gamma_2pi);
    
    if P.B1cwpe_quad==1
        %b1 quadratisch
        norm=sqrt(trapz(tpulse,abs(w1.^2))./((P.tp+P.td*( (P.n-1)/P.n ))));
        w1=w1./norm*P.B1;
    end;
end;

theta=trapz(tpulse,w1*gamma_2pi);
B1q=sqrt( trapz(tpulse,w1.^2)./((P.tp/P.DC)) );

%  fprintf('mittleres B1squared =%f\n',B1q);
%  fprintf('mittleres B1lin =%.4f, flipwinkel=%.2f (%.2f) \n',trapz(tpulse,abs(w1))./(P.tp+P.td),trapz(tpulse,abs(w1))*gamma_*360,trapz(tpulse,abs(w1))/P.tp );
% fprintf('innerpulse duty-cycle: %.2f',sum( tpulse((abs(w1)>P.B1)) )/sum(tpulse));