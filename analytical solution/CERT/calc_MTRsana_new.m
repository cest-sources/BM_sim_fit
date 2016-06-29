function M=calc_MTRsana_new(P,Pz,R1rho,Psi,theta)

% pause constants
r1a=P.R1A+P.kAB;
r1b=P.R1B+P.kBA;
a=1/2*(r1b+r1a+sqrt((r1b-r1a)^2+4*P.kBA*P.kAB));
b=1/2*(r1b+r1a-sqrt((r1b-r1a)^2+4*P.kBA*P.kAB));
dAA=1/(a-b)*(-(b-r1a)*exp(-a*P.td)+(a-r1a)*exp(-b*P.td));
dAB=P.kBA/(a-b)*(-exp(-a*P.td)+exp(-b*P.td));
dBA=P.fB*dAB;
dBB=1/(a-b)*((a-r1a)*exp(-a*P.td)-(b-r1a)*exp(-b*P.td));


% pulse constants for Pool A
p=exp(-R1rho*P.tp); %Zerfall
q=(1-p).*(Pz.*cos(theta).*P.R1A./(R1rho)-1);

% pulse constants for Pool B
Phi.A=Psi.A+zeros(1,length(theta)); %factor for dMza
Phi.B=Psi.B+zeros(1,length(theta)); %factor for dMzb
Phi.C=Psi.C+Psi.A+Psi.B.*P.fB-P.fB; %constant factor

%calculating the matrix constants
a11=p.*dAA;
a12=p.*dAB;
a21=Phi.A.*dAA+Phi.B.*dBA;
a22=Phi.A.*dAB+Phi.B.*dBB;

%calculation of the steady state
alpha=1./(1-a11-a12.*a21-a22+a11.*a22);
M.ss(1,:) =alpha.*((1-a22).*q+a12.*Phi.C); %steady state Pool A
M.ss(2,:) =alpha.*((a21).*q+(1-a11).*Phi.C); %steady state Pool B

%calculation of the A^n matrix
for i=1:length(theta)
   A=[a11(i) a12(i);
      a21(i) a22(i)];
   An=A^(P.n);
   M.dynamic(:,i)=([1 0; 0 1]-An)*M.ss(:,i);
end

%rescaling
M.ss=M.ss+[1+zeros(1,length(theta));P.fB+zeros(1,length(theta))];
M.dynamic=M.dynamic+[1+zeros(1,length(theta));P.fB+zeros(1,length(theta))];
end