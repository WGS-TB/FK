% This function generates D^{m,x}_1 or C^{m,x}_0 by setting the variables in m to TRUE and calculate the remaining formula.
% Please see "how to apply sat-solving for the equivalence test of mootone normal forms", page 4
% This function inputs are:
% 1) A: D_x_1 or C_x_0
% 2) m: a vector which each element is one if that variable must be set to TRUE
% 3) type: It indicates to the form of A; It can be CNF or DNF
% Returns D^{m,x}_1 or C^{m,x}_0
function res = A_m_x(A, m, type)

if(strcmp(type,'CNF'))

    B = A;
    for j = 1 : size(A,1)
        B(j, :) = and(m, A(j,:));
    end
    chk = sum(B, 2) > 0;
    if (~isempty(chk))
        A(chk,:) =[];
    end
    if isempty(A)
        display('This is A_m_x and CNF got TRUE!!')
%         stop()
        res = true;
    else
        res = A;
    end
    
else % Type = DNF
    
    A(:,logical(m)) = 0;
    if sum(sum(A==0, 2) == size(A,2))>0
        % if one row of DNF gets empty, it makes DNF true
        display('This is A_m_x and DNF got TRUE!!')
        stop
        res = false;
    else
        res = A;
    end
    
end
end