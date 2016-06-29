function [dep_vars, start, lower, upper  ] = selectVars( vary, dep_vars, start, lower, upper )
% comments here
% last change: 2014/04/09 by PS

dep_vars = dep_vars(vary==1);
start    = start(vary==1);
lower    = lower(vary==1);
upper    = upper(vary==1);

end

