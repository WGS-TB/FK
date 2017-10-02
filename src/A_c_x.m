% This function generates C^{c,x}_1 or D^{c,x}_0 by setting the variables in c to FALSE and calculate the remaining formula.
% Please see "how to apply sat-solving for the equivalence test of mootone normal forms", page 4
% This function inputs are:
% 1) A: C_x_1 or D_x_0 in compressed format
% 2) cc: the set of variables which must be set to FALSE
% 3) type: It indicates to the form of A; It can be CNF or DNF
% Returns C^{c,x}_1 or D^{c,x}_0
function res = A_c_x(A, cc, type)

if(strcmp(type,'CNF'))
    
    A(:,logical(cc)) = 0;
    if sum(sum(A==0, 2) == size(A,2))>0
        % if one row of CNF gets empty, it makes CNF false
        display('This is A_c_x and CNF got FALSE!!')
        res = false;
    else
        res = A;
    end
else % Type = DNF
    
    B = A;
    for j = 1 : size(A,1)
        B(j, :) = and(cc, A(j,:));
    end
    chk = sum(B, 2) > 0;
    A(chk,:) =[];

    % If a monomial is empty of any variable we remove that monomial
    if isempty(A)
        display('This is A_cc_x and DNF got FALSE!!')
%         stop()
        res = false;
    else
        res = A;
    end
end
end
