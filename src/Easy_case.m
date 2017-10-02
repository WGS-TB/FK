function CA = Easy_case(cnf,dnf)

CA = [];
% +++++++++++++++++++++++++++++++++++++++++++++++++
% Case like cnf=(1,2,3) and dnf= (1) (2) (3)
% +++++++++++++++++++++++++++++++++++++++++++++++++

if(size(cnf,1)==1 && size(dnf,1) == sum(sum(cnf)) && sum(sum(cnf))== sum(sum(dnf)))
    [~,dnf_vars] = find(dnf==1);
    cnf_vars = find(cnf==1);
    if(isequal(cnf_vars, dnf_vars'))
        CA = [];
        return
    end
end

if(size(dnf,1)==1 && size(cnf,1) == sum(sum(dnf)) && sum(sum(cnf))== sum(sum(dnf)))
    [~,cnf_vars] = find(cnf==1);
    dnf_vars = find(dnf==1);
    if(isequal(cnf_vars, dnf_vars))
        CA = [];
        return
    end
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Checking intersection between clauses and monomials
% If there is no intersection, conflict assignment would be the monomial
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% res = cnf * dnf';
% [r, c] = find(res==0);
% ind = [r c];
% % if ~isempty(ind)=0, it means that there is no zero in the res and
% % there is intersection between every monomials and clauses
% if (~isempty(ind))
%     CA = zeros(1, size(dnf,2));
%     CA(dnf(ind(1,1),:)) = 1;
% %     if (~isempty(CA))
% %         disp('Easy 1')
% %     end
%
%     return
% end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case I:  When DNF has at most two monomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(dnf,1)<=2)
    % browser()
    res = CA_CNF2(dnf, cnf);
    if(isempty(res))
        CA = [];
    else
        [~,vars]=find(dnf);
        vars = unique(vars);
        % % %         temp = setdiff(vars, res);
        
        CA = zeros(size(res,1), size(dnf,2));
        for tt=1:size(res,1)
        temp = setdiff(vars, find(res(tt,:)));
        CA(tt,temp) = 1;
        end
                if (~isempty(CA))
                    disp('Easy 2')
                end
        
    end
    return
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Case II:  When DNF has more than two monomials and CNF is the shorter one
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp =  CA_CNF2(cnf, dnf);
    if(~isempty(temp))
        %%%        CA = zeros(1, size(dnf,2));
        % % %         CA(temp) = 1;
        CA=temp;
    end
        if (~isempty(CA))
            disp('Easy 3')
        end
    
    return
end
