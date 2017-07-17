function y=bfilt(x)

% Blaise Filter
x=x(:)';
l=length(x);
% x=x';

m(1,:)=[2*x(1)-x(3) 2*x(1)-x(2) x(1:l-2)];
m(2,:)=[2*x(1)-x(2) x(1:l-1)];
m(3,:)=x;
m(4,:)=[x(2:l) 2*x(l)-x(l-1)];
m(5,:)=[x(3:l) 2*x(l)-x(l-1) 2*x(l)-x(l-2)];

b=[0 .2 0.6 .2 0];
y=b*m;

