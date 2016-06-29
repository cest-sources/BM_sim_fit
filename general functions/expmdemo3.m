function [E dia]= expmdemo3(A)
%EXPMDEMO3  Matrix exponential via eigenvalues and eigenvectors.
%   E = EXPMDEMO3(A) illustrates one possible way to compute the matrix
%   exponential.  As a practical numerical method, the accuracy
%   is determined by the condition of the eigenvector matrix.
%
%   See also EXPM, EXPMDEMO1, EXPMDEMO2.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2004/08/16 01:37:28 $

[V,D] = eig(A);
dia=sort(diag(D));

E = V * diag(exp(diag(D))) / V;
