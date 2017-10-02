% This function selects variable x which is most frequent in both CNF and DNF
function res = Choose_SplitVar( cnf, dnf, method )

if strcmp(method , 'mostFreq')
    [~, cnfVar] = sort(sum(cnf, 1), 'descend');
    [~, dnfVar] = sort(sum(dnf, 1), 'descend');
    
    for i=1:numel(cnfVar)
        var = intersect(cnfVar(1:i), dnfVar(1:i));
        if ~isempty(var)
            res = var(1);
            return
        end
    end
end
end

