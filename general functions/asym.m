function [xasym, yasym] = asym(xx,ZZ)
% calculates Asym-spectrum from z-spectrum
% input: x-values (xx), z-values (ZZ)
% output: x-values (xasym), asym-values (yasym)
% last change: 2014/04/02 by PS


if ~mod(numel(xx),2) %if even number of elements = no zero
    int1=numel(xx)/2:-1:1;
    int2=numel(xx)/2+1:numel(xx);   
else
    int1=fix(numel(xx)/2):-1:1;
    int2=fix(numel(xx)/2)+2:numel(xx);
end;

xasym(1)=0;
yasym(1)=0;
xasym(2:numel(int2)+1)=xx(int2);
yasym(2:numel(int2)+1)=(ZZ(int1)-ZZ(int2));





