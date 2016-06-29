function M=calc_Mag_ana(P,S,Pz,Pzeff,theta,Psi)
    %constants
    r1a=P.R1A+P.kAB;
    r1b=P.R1B+P.kBA;
    a=1/2*(r1b+r1a+sqrt((r1b-r1a)^2+4*P.kBA*P.kAB));
    b=1/2*(r1b+r1a-sqrt((r1b-r1a)^2+4*P.kBA*P.kAB));
    dAA=1/(a-b)*(-(b-r1a)*exp(-a*P.td)+(a-r1a)*exp(-b*P.td));
    dAB=P.kBA/(a-b)*(-exp(-a*P.td)+exp(-b*P.td));
    dAC=1-dAA-dAB*P.fB;
    dBA=P.fB*dAB;
    dBB=1/(a-b)*((a-r1a)*exp(-a*P.td)-(b-r1a)*exp(-b*P.td));
    dBC=(P.fB-dBA-dBB*P.fB);


    % pulse constants for Pool A
    M=zeros(2,length(S.Rho_Zaiss_Hyper_full));
    for i=1:length(S.Rho_Zaiss_Hyper_full)
        if strcmp(P.shape,'SPINLOCK')
            p=exp(S.Rho_Zaiss_Hyper_full(i)*P.tp); %Zerfall
            q=(1-p).*(Pz.*cos(theta(i)).*P.R1A./(-S.Rho_Zaiss_Hyper_full(i)));
            A=[dAA dAB; dBA dBB]*[p*Pz*Pzeff 0; Psi.A(i) Psi.B(i)];
            C=[q+p*Pz*Pzeff*dAC; Psi.A(i)*dAC+Psi.B(i)*dBC+Psi.C(i)];
        else 
            p=exp(S.Rho_Zaiss_Hyper_full(i)*P.tp); %Zerfall
            q=(1-p).*(Pz(i).*cos(theta(i)).*P.R1A./(-S.Rho_Zaiss_Hyper_full(i)));
            A=[p*Pz(i)*Pzeff(i) 0; Psi.A(i) Psi.B(i)]*[dAA dAB; dBA dBB];
            C=[q+p*Pz(i)*Pzeff(i)*dAC; Psi.A(i)*dAC+Psi.B(i)*dBC+Psi.C(i)];
        end
        %calculation of the steady state
        alpha=1./(1-A(1,1)-A(1,2).*A(2,1)-A(2,2)+A(1,1).*A(2,2));
        Mss =alpha.*[((1-A(2,2)).*C(1)+A(1,2).*C(2));(A(2,1)).*C(1)+(1-A(1,1)).*C(2)]; 
        
        M(:,i)=A^(P.n)*([1;P.fB]-Mss)+Mss;
    end
end
    
    