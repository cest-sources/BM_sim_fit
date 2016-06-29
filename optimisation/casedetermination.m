function [dep_vars, startValue, lowerbounds, upperbounds] = casedetermination(P, T)
% determine which case (--> parameters) is used
% last change: 2014/04/04 by PS

if P.MT
    switch P.n_cest_pool
        case 0
            dep_vars    = [T.dep_varsA T.dep_varsC];
            startValue  = [T.startA T.startC];
            lowerbounds = [T.lowerA T.lowerC];
            upperbounds = [T.upperA T.upperC];
        case 1
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsC];
            startValue  = [T.startA T.startB T.startC];
            lowerbounds = [T.lowerA T.lowerB T.lowerC];
            upperbounds = [T.upperA T.upperB T.upperC];
        case 2
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsC T.dep_varsD];
            startValue  = [T.startA T.startB T.startC T.startD];
            lowerbounds = [T.lowerA T.lowerB T.lowerC T.lowerD];
            upperbounds = [T.upperA T.upperB T.upperC T.upperD];
        case 3
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsC T.dep_varsD T.dep_varsE];
            startValue  = [T.startA T.startB T.startC T.startD T.startE];
            lowerbounds = [T.lowerA T.lowerB T.lowerC T.lowerD T.lowerE];
            upperbounds = [T.upperA T.upperB T.upperC T.upperD T.upperE];
        case 4
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsC T.dep_varsD T.dep_varsE T.dep_varsF];
            startValue  = [T.startA T.startB T.startC T.startD T.startE T.startF];
            lowerbounds = [T.lowerA T.lowerB T.lowerC T.lowerD T.lowerE T.lowerF];
            upperbounds = [T.upperA T.upperB T.upperC T.upperD T.upperE T.upperF];
        case 5
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsC T.dep_varsD T.dep_varsE T.dep_varsF T.dep_varsG];
            startValue  = [T.startA T.startB T.startC T.startD T.startE T.startF T.startG];
            lowerbounds = [T.lowerA T.lowerB T.lowerC T.lowerD T.lowerE T.lowerF T.lowerG];
            upperbounds = [T.upperA T.upperB T.upperC T.upperD T.upperE T.upperF T.upperG];
    end
else
    switch P.n_cest_pool
        case 0
            dep_vars    = [T.dep_varsA];
            startValue  = [T.startA];
            lowerbounds = [T.lowerA];
            upperbounds = [T.upperA];
        case 1
            dep_vars    = [T.dep_varsA T.dep_varsB];
            startValue  = [T.startA T.startB];
            lowerbounds = [T.lowerA T.lowerB];
            upperbounds = [T.upperA T.upperB];
        case 2
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsD];
            startValue  = [T.startA T.startB T.startD];
            lowerbounds = [T.lowerA T.lowerB T.lowerD];
            upperbounds = [T.upperA T.upperB T.upperD];
        case 3
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsD T.dep_varsE];
            startValue  = [T.startA T.startB T.startD T.startE];
            lowerbounds = [T.lowerA T.lowerB T.lowerD T.lowerE];
            upperbounds = [T.upperA T.upperB T.upperD T.upperE];
        case 4
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsD T.dep_varsE T.dep_varsF];
            startValue  = [T.startA T.startB T.startD T.startE T.startF];
            lowerbounds = [T.lowerA T.lowerB T.lowerD T.lowerE T.lowerF];
            upperbounds = [T.upperA T.upperB T.upperD T.upperE T.upperF];
        case 5
            dep_vars    = [T.dep_varsA T.dep_varsB T.dep_varsD T.dep_varsE T.dep_varsF T.dep_varsG];
            startValue  = [T.startA T.startB T.startD T.startE T.startF T.startG];
            lowerbounds = [T.lowerA T.lowerB T.lowerD T.lowerE T.lowerF T.lowerG];
            upperbounds = [T.upperA T.upperB T.upperD T.upperE T.upperF T.upperG];
    end
end

