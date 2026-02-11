function [M0_2, T2map, c_2, S_fit] = t2fitting_new(S, TEs)
% =========================================================================
%                               T2 fitting (alias).
% =========================================================================
% This is an alias function that calls t2fitting.m
% Created for compatibility with phantom_T2.m

[M0_2, T2map, c_2, S_fit] = t2fitting(S, TEs);

end
