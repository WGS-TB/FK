
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% By: Nafiseh Sedaghat.
% nf_sedaghat@sfu.ca
% Copyright (C) 2018.
%
% This code was developed at Simon Fraser University.
%
% Version of February 2018.
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  The GNU General Public License is available on-line at:
%   <http://www.gnu.org/licenses/>.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
% The below code has been written based on the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5
% This script is combination of FK_M and FK_R
% For more information please see "Speeding up the Fredman-Khachiyan Dualization Algorithm B" paper.
%
% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.
% The number of columns in cnf and dnf must be the same.
% cnf_redundancy: A flag that can be 0 or 1. If it is 0, the redundancy of cnf is not checked, otherwise it is checked.
% dnf_redundancy: A flag that can be 0 or 1. If it is 0, the redundancy of dnf is not checked, otherwise it is checked.
% split_method: the method of choosing split variable; the only option is 'mostFreq' and it is for further improvements.
%
% Output:
% CA: It is the conflict assignment between the given cnf and dnf. The variables corresponding to '1' in CA are considered as True and the other variables that might be '0' or '-1' are considered as False.
% If the given cnf and dnf are equivalent, CA would be [] (NULL).
%
% Example:
% cnf = [1 1 0 0 0; 1 0 1 0 1;0 1 1 1 0];
% dnf = [1 0 0 0 0; 0 1 1 0 0];
% FK_MR(cnf, dnf, 1, 1, 'mostFreq')
%
% cnf = [1 1 0 0 1 0; 1 0 1 0 0 1;0 1 1 1 0 1];
% dnf = [1 0 0 0 1 1; 0 1 1 0 1 0];
% FK_MR(cnf, dnf, 1, 1, 'mostFreq')

function CA = FK_MR(cnf, dnf, cnf_redundancy, dnf_redundancy, split_method)

global call_counter redundancy_counter
CA = [];
call_counter = call_counter +1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                     Checking redundancy                 %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (dnf_redundancy == 1)
    if(size(dnf,1) > 1)
        dnf = Irredundant(dnf);
        redundancy_counter = redundancy_counter +1;
    end
end

if (cnf_redundancy ==1)
    if(size(cnf,1) > 1)
        cnf = Irredundant(cnf);
        redundancy_counter = redundancy_counter +1;
        
    end
end

% write2file(cnf, dnf);

nC = size(cnf,1); % Number of clauses in CNF
nD = size(dnf,1); % Number of monomials in DNF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                     Checking conditions                 %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CA = Multiple_Check_Conditions(cnf, dnf);

if (~isempty(CA))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                 Trivial Check                           %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(nD, nC) <= 2
    val = MEasy_case(cnf, dnf);
    if(~isempty(val))
        CA = val;
        return
    end
else %line 5
    
    split_var = Choose_SplitVar(cnf, dnf, 'mostFreq');
    
    D_x_1 = phi_x_1(dnf, split_var);
    D_x_0 = phi_x_0(dnf, split_var);
    C_x_1 = phi_x_1(cnf, split_var);
    C_x_0 = phi_x_0(cnf, split_var);
    
    if mu_frequent_in_A(split_var, dnf, cnf)% if split variable is at most mu-frequent in DNF (line 7)
        CA = FK_MR([C_x_0; C_x_1], D_x_1, 1, 0, split_method ); % line 8
        
        if (~isempty(CA))
            return %line 9
        end
        
        for i=1: size(C_x_0,1) %line 10
            vars = C_x_0(i,:);
            
            if(size(D_x_0,1)>0)
                D_cx_0 = A_c_x(D_x_0, vars, 'DNF');
            else
                D_cx_0 = [];
            end
            
            if(size(C_x_1,1)>0)
                C_cx_1 = A_c_x(C_x_1, vars, 'CNF');
            else
                C_cx_1 = [];
            end
            
            if(isempty(D_cx_0) && isempty(C_cx_1))
                CA = [];
            else
                CA = FK_MR(C_cx_1, D_cx_0, 1, 1, split_method); %line 11
            end
            
            if (~isempty(CA))
                CA(:, split_var) = 1; % adding split_var
                return %line 12
            end
        end %end of for (i in 1: length(C_x_0))
    else if mu_frequent_in_A(split_var, cnf, dnf) % if split variable is at most mu-frequent in CNF(line 13)
            CA = FK_MR(C_x_1, [D_x_0; D_x_1], 0, 1, split_method); %line 14
            if (~isempty(CA))
                CA(:, split_var) = 1;
                return %line 15
            end
            
            for i=1: size(D_x_0,1) %line 16
                vars = D_x_0(i,:);
                
                if(size(D_x_1,1)>0)
                    D_mx_1 = A_m_x(D_x_1, vars, 'DNF');
                else
                    D_mx_1 = [];
                end
                if(size(C_x_0,1)>0)
                    C_mx_0 = A_m_x(C_x_0, vars, 'CNF');
                else
                    C_mx_0 = [];
                end
                
                if(isempty(D_mx_1) && isempty(C_mx_0))
                    CA = [];
                else
                    CA = FK_MR(C_mx_0, D_mx_1, 1, 1, split_method); %line 17
                end
                
                if (~isempty(CA))
                    CA(:, logical(vars))= 1;
                    return
                end
            end %end of for (i in 1: length(C_x_0))
        else % line 19
            
            CA = FK_MR([C_x_0; C_x_1], D_x_1, 1, 0, split_method);  %line 20
            
            if(isempty(CA)) % line 21
                
                CA = FK_MR(C_x_1, [D_x_0; D_x_1], 0, 1, split_method); % line 22
                
                if (~isempty(CA))
                    % browser()
                    CA(:, split_var) =1; % adding split.var
                    return %line 23
                end
                
            end% end if if(length(CA)==0)
        end% end else if(mu_frequent_in_A(split.var, A=cnf, B=dnf))
    end% end else %line 5
    return
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function remove redundancy, i.e. supersets, in a given monotone Boolean function and returns an irredundant monotone Boolean function.
% MBF is a monotone Boolean function as a Binary matrix. The type of MBF, i.e. CNF or DNF, does not matter during the procedure of removing redundancy.
function MBF = Irredundant(MBF)

nrow = size(MBF, 1);
chk = false(1, nrow); % chk(i) is one if MBF(i,:) is a superset and should be removed

% Sort rows in ascending order based on the number of elements in each clause/monomial
summation = sum(MBF, 2);
[~,I] = sort(summation);
MBF = MBF(I, :);

for i=1:(nrow-1)
    a = MBF(i, :);
    for j=(i+1):nrow
        % This is the bottoleneck!!
        if ~chk(j) % if clause/monomial j is not already checked for being removed
            chk(j) = isequal(and(MBF(j,:), a), a);
        end
    end
end
MBF(chk,:) = [];

return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function check the three conditions mentioned in FK-B algorithm. The conditions are:
% 1: for each c in CNF and m in DNF => intersect(m, c) != null
% 2: DNF and CNF must have exactly the same variables
% 3: max{|m|: m \in DNF} <= |C|   &   max{|c|: c \in CNF} <= |D|

% If one of these conditions is not met in the given cnf and dnf, the function returns a cnoflict assignment.
% Inputs:
% cnf and dnf are binary matrices indicating to a monotone boolean function in the form of conjuctive normal form and disjunctive normal form.

% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.
% The number of columns in cnf and dnf must be the same.

% Output:
% CA: It is the conflict assignment(s) between the given cnf and dnf, each row corresponds to a CA. The variables corresponding to '1' in CA are considered as True and the other variables that might be '0' or '-1' are considered as False.
% If the given cnf and dnf are equivalent, CA would be [] (NULL).

function CA = Multiple_Check_Conditions(cnf, dnf)

CA = [];
%++++++++++++++++++++++++++++
%
%++++++++++++++++++++++++++++
if (isempty(cnf) || isempty(dnf))
    
    if(isempty(dnf)) % DNF is NULL, then DNF is false and we want CNF to be TRUE
        vars = 0;
        for i=1:size(cnf,1)
            vars(i) = find(cnf(i,:)==1, 1, 'first');
        end
        temp = zeros(1, size(cnf, 2));
        temp(vars) = 1;
        CA = temp;
        
    else % CNF is NULL, then CNF is true and we want DNF to be False
        [row, vars] = find(dnf==1);
        vv = zeros(1, size(dnf, 1));
        for t=1:size(dnf,1)
            ind=find(row==t, 1, 'first');
            vv(t) = vars(ind);
        end
        vars = unique(vv);
        temp = zeros(1, size(dnf, 2));
        temp(vars) = -1; % False vars
        CA = temp;
    end
    return;
end
%
%++++++++++++++++++++++++++++
%       First Condition
%++++++++++++++++++++++++++++

chk = zeros(size(cnf, 1), size(dnf, 1));
chk = (cnf * dnf') == 0;
if(sum(sum(chk))>0)
    [~, y] = find(chk==0);
    CA = unique(dnf(y, :), 'rows');
    return
end


%++++++++++++++++++++++++++++
%       Second Condition
%++++++++++++++++++++++++++++
%

[~, temp] = find(cnf);
cnf_vars = unique(temp);

[~, temp] = find(dnf);
dnf_vars = unique(temp);

check = isequal(cnf_vars, dnf_vars) || isequal(cnf_vars, dnf_vars');

if (~check)
    x = setdiff(dnf_vars, cnf_vars);
    if (numel(x)>0) % the extra variables are in dnf
        
        [~, temp] = find(dnf);
        dnf_vars = unique(temp);
        
        x = setdiff(dnf_vars, cnf_vars);
        
        [r,~] = find(dnf(:,x));
        
        temp = dnf(r,:);
        temp(:, x) = 0;
        CA = unique(temp, 'rows');
        return
        
    else % the extra variables are in cnf
        
        x = setdiff(cnf_vars, dnf_vars);
        [~, temp] = find(cnf);
        cnf_vars = unique(temp);
        
        [r,~] = find(cnf(:,x));
        tt = zeros(numel(r), size(cnf, 2));
        tt(:,cnf_vars) = 1;
        temp = tt - cnf(r,:);
        temp(:,x) = 1;
        CA = unique(temp, 'rows');
    end
    return;
end

%++++++++++++++++++++++++++++
%       Third Condition
%++++++++++++++++++++++++++++

t1 = max(sum(dnf, 2)) <= size(cnf, 1); % max{|m|: m \in DNF} <= |C|
t2 = max(sum(cnf, 2)) <= size(dnf, 1); % max{|c|: c \in CNF} <= |D|

check =  t1 && t2;

if (~check)
    % We find the conflict assignment based on waht has been said in original FK_B paper in page 2, indented part.
    if (~t1)
        [~, I] = max(sum(dnf, 2));
        longest_dnf = dnf(I,:);
        
        for i=1:((2^sum(longest_dnf))-1) % We are seeking for proper subset of longest monomial in dnf
            
            bin = dec2bin(i, sum(longest_dnf));
            S = zeros(1, size(dnf, 2));
            bb = zeros(1, length(bin));
            for u=1:length(bin)
                bb(u) = str2double(bin(u));
            end
            S(logical(longest_dnf)) = bb;
            
            chk = zeros(1, size(cnf,1));
            for j=1:size(cnf,1)
                chk(j) = sum(and(S,cnf(j,:)))>0;
            end
            if(~ismember(0, chk))
                CA = S; % we must set *proper subset* of the largest monoials to true
                return
            end
        end
    end
    
    if (~t2)
        [~, I] = max(sum(cnf, 2));
        longest_cnf = cnf(I,:);
        
        for i=1:((2^sum(longest_cnf))-1) % We are seeking for proper subset of longest monomial in dnf
            
            bin = dec2bin(i, sum(longest_cnf));
            S = zeros(1, size(cnf, 2));
            bb = zeros(1, length(bin));
            for u=1:length(bin)
                bb(u) = str2double(bin(u));
            end
            S(logical(longest_cnf)) = bb;
            
            chk = zeros(1, size(dnf,1));
            for j=1:size(dnf,1)
                chk(j) = sum(and(S,dnf(j,:)))>0;
            end
            if(~ismember(0, chk))
                CA = (sum(cnf,1)>0) - S; % we must set *proper subset* of the largest monoials to true
                return
            end
        end
    end
    
end

% If we reach here, it means that all of conditions are satisfied.
return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function is called when one of the given cnf or dnf has less than 3 clauses or monomials and it checks equivalence between cnf and dnf based on Boolean algebra.

% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.
% The number of columns in cnf and dnf must be the same.

% Output:
% CA: It is the conflict assignment(s) between the given cnf and dnf; each row corresponds to a conflict assignment. The variables corresponding to '1' in CA are considered as True and the other variables that might be '0' or '-1' are considered as False.
% If the given cnf and dnf are equivalent, CA would be [] (NULL).

function CA = MEasy_case(cnf,dnf)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case I:  When DNF has at most two monomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(size(dnf,1)<=2)
    res = CA_smallCNF(dnf, cnf);
    if(isempty(res))
        CA = [];
    else
        [~,vars]=find(dnf);
        vars = unique(vars);
        
        CA = zeros(size(res,1), size(dnf,2));
        for tt=1:size(res,1)
            temp = setdiff(vars, find(res(tt,:)));
            CA(tt,temp) = 1;
        end
    end
    return
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Case II:  When DNF has more than two monomials and CNF is the shorter one
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp =  CA_smallCNF(cnf, dnf);
    if(~isempty(temp))
        CA=temp;
    end
    return
end
end
%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

%This function returns the conflict assignment between cnf and dnf when |cnf|, i.e. the number of clauses in cnf, is <= 2. CA is directly computed based on Boolean algebra.
% It is called by Easy_case function.

% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.
% The number of columns in cnf and dnf must be the same.

% Output:
% CA: It is the conflict assignment between the given cnf and dnf. The variables corresponding to '1' in CA are considered as True and the other variables that might be '0' or '-1' are considered as False.
% If the given cnf and dnf are equivalent, CA would be [] (NULL).

function CA =  CA_smallCNF( cnf, dnf)

CA = [];

cnf_vars = find(sum(cnf, 1));

%++++++++++++++++++++
%  Case II-A: length(cnf)==1
%++++++++++++++++++++
if(size(cnf,1)==1)
    for i=1:(2^numel(cnf_vars)-1) % We are seeking for a subset
        bin = dec2bin(i, numel(cnf_vars));
        S = zeros(1, size(cnf, 2));
        bb = zeros(1, length(bin));
        for u=1:length(bin)
            bb(u) = str2double(bin(u));
        end
        S(cnf_vars) = bb;
        
        chk = zeros(1, size(dnf,1));
        for j=1:size(dnf, 1)
            chk(j) = subsetCheck(dnf(j,:), S);
        end
        
        if (sum(chk)==0) % no monomial is subset of S
            CA = S;
            return
        end
    end
end
%++++++++++++++++++++++++++++
% Case II-B:  length(cnf)==2
%++++++++++++++++++++++++++++
if(size(cnf,1)==2)
    
    %## Case 1
    S = and(cnf(1,:), cnf(2,:));
    
    S1 = 1*((cnf(1,:) - cnf(2,:))>0);
    S2 = 1*((cnf(2,:) - cnf(1,:))>0);
    
    if(sum(S1)>0 && sum(S2)>0)
        cart_mult = cartprod(find(S1), find(S2)); %Cartsian Multiplication
        temp = zeros(size(cart_mult,1), size(cnf,2));
        out = zeros(1, size(cart_mult,1));
        
        for i=1:size(cart_mult,1)
            temp(i, cart_mult(i,:))=1;
            out(i) = any(ismember(dnf,temp(i,:),'rows'));
        end
        
        % We are looknig for a member of cart_mult which is not appered in dnf
        ind = find(out==0);
        if(numel(ind)>0)
            %             CA = find(temp(ind, :));
            CA = temp(ind, :);
            return
        end
        if (sum(S)==0)
            CA = [];
            return
        end
    end
    
    %## Case 2
    for i=1:((2^sum(S))-1) % We are seeking for proper subset of A
        bin = dec2bin(i, sum(S));
        binNum = zeros(1, size(S, 2));
        xs = zeros(1,numel(bin));
        for jj=1:numel(bin)
            xs(jj)= str2double(bin(jj));
        end
        binNum(logical(S))= xs;
        
        chk = zeros(1, size(dnf, 1));
        for l=1:size(dnf,1)
            chk(l) = subsetCheck(dnf(l,:), binNum);
        end
        
        if (sum(chk)==0) % no monomial is subset of S
            CA = binNum;
            return
        end
    end
    
    %## Case 3
    
    S = (cnf(1,:) - cnf(2,:)) > 0;
    out = zeros(1, size(dnf,1));
    for i=1:size(dnf,1)
        out(i) = subsetCheck(dnf(i,:), S);
    end
    
    ind = find(out==1, 1, 'first');
    if(~isempty(ind))
        CA = S;
        return
    end
    
    %## Case 4
    
    S = (cnf(2,:) - cnf(1,:)) > 0;
    out = zeros(1, size(dnf,1));
    for i=1:size(dnf,1)
        out(i) = subsetCheck(dnf(i,:), S);
    end
    
    ind = find(out==1, 1, 'first');
    if(~isempty(ind))
        CA = S;
        return
    end
    
    CA = [];
    return
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function selects variable x which is most frequent in both CNF and DNF

% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.
% The number of columns in cnf and dnf must be the same.
% method: The method of choosing the splitting variable. Currently, the only method is 'mostFreq', and it is planted for further extension.

% Output:
% res: the chosen variable

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


%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% phi_x_0 denotes the formula that consists of terms of phi from which x is removed.
% For further information, please see the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5

% Inputs:
% MBF: Monotone Boolean Function as a binary matrix
% var: the variable that is supposed to be removed.

% Output:
% res: The monotone boolean function after removing variable 'var'.

function res = phi_x_0(MBF, var)

[rows , ~]= find(MBF(:,var)==1);
rows = unique(rows);
MBF = MBF(rows,:); % Only those rows that contain var
MBF(:,var) = 0; % removing variable var from clauses that var was included
% If a clause is empty of any variable we remove that clause

rem = sum(MBF, 2)==0;
if (sum(rem)>0)
    MBF(rem, :) = [];
end

res = MBF;

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% phi_x_1 denotes the formula that consists of all terms of phi that do not contain x.

% For further information, please see the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5

% Inputs:
% MBF: Monotone Boolean Function as a binary matrix
% var: the variable that is supposed to be removed.

% Output:
% res: The monotone boolean function consists of all terms of MBF that do not contain 'var'.

function res = phi_x_1(MBF, var)

% [rows , ~]= find(MBF(:,var)==1);
% rows = unique(rows);
% MBF = MBF(setdiff(1:size(MBF,1), rows),:); % Only those rows that does not contain var

[rows , ~]= find(MBF(:,var)==0);
rows = unique(rows);
MBF = MBF(rows,:); % Only those rows that does not contain var

% If a clause is empty of any variable we remove that clause
rem = sum(MBF, 2)==0;
if (sum(rem)>0)
    MBF(rem, :) = [];
end

res = MBF;

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function checks whether split variable satisfy the below condition:
% |{m \in A: x \in m}|/|A| <= 1/mu(|A| . |B|)
% which we say that x is at most frequent in A
% For further information, please see the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5

% Inputs:
% split_var: the variable chosen as splitting variable
% A: A monotone boolean function as a binary matrix
% B: A monotone boolean function as a binary matrix
% The number of columns in A and B must be the same.

% Output: True or False

function res = mu_frequent_in_A( split_var, A, B )

ss = sum(A, 1);
len = size(A, 1);
left_side = ss(split_var) / len;
right_side = 1/mu_function(len * size(B,1));

res = left_side <= right_side;
return

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% mu_function is needed for choosing splitting variable and is called by 'mu_frequent_in_A' function.
% For further information, please see the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5

function res = mu_function(n)

res = log(n)/log(log(n));
return

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function generates C^{c,x}_1 or D^{c,x}_0 by setting the variables in c to FALSE and calculate the remaining formula.
% For further information, please see "how to apply sat-solving for the equivalence test of mootone normal forms", page 4

% Inputs:
% A: C_x_1 or D_x_0 as a binary matrix
% cc: the set of variables which must be set to FALSE
% type: It indicates to the form of A; It can be 'CNF' or 'DNF'

% Output:
% res would be C^{c,x}_1 or D^{c,x}_0

function res = A_c_x(A, cc, type)

if(strcmp(type,'CNF'))
    
    A(:,logical(cc)) = 0;
    if sum(sum(A==0, 2) == size(A,2))>0
        % if one row of CNF gets empty, it makes CNF false
        res = [];
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
        res = [];
    else
        res = A;
    end
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function generates D^{m,x}_1 or C^{m,x}_0 by setting the variables in m to TRUE and calculate the remaining formula.
% For further information, please see "how to apply sat-solving for the equivalence test of mootone normal forms", page 4

% Inputs:
% A: D_x_1 or C_x_0
% m: a vector which each element is one if that variable must be set to TRUE
% type: It indicates to the form of A; It can be 'CNF' or 'DNF'.

% Output:
% res would be D^{m,x}_1 or C^{m,x}_0

function res = A_m_x(A, m, type)

if(strcmp(type,'CNF'))
    
    B = A;
    for j = 1 : size(A,1)
        B(j, :) = and(m, A(j,:));
    end
    chkk = sum(B, 2) > 0;
    if (~isempty(chkk))
        A(chkk,:) =[];
    end
    if isempty(A)
        res = [];
    else
        res = A;
    end
    
else % Type = DNF
    
    A(:,logical(m)) = 0;
    if sum(sum(A==0, 2) == size(A,2))>0
        % if one row of DNF gets empty, it makes DNF true
        res = [];
    else
        res = A;
    end
    
end
end
%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function has been downloaded from
% https://www.mathworks.com/matlabcentral/fileexchange/5475-cartprod--cartesian-product-of-multiple-sets
function X = cartprod(varargin)
%CARTPROD Cartesian product of multiple sets.
%
%   X = CARTPROD(A,B,C,...) returns the cartesian product of the sets 
%   A,B,C, etc, where A,B,C, are numerical vectors.  
%
%   Example: A = [-1 -3 -5];   B = [10 11];   C = [0 1];
% 
%   X = cartprod(A,B,C)
%   X =
% 
%     -5    10     0
%     -3    10     0
%     -1    10     0
%     -5    11     0
%     -3    11     0
%     -1    11     0
%     -5    10     1
%     -3    10     1
%     -1    10     1
%     -5    11     1
%     -3    11     1
%     -1    11     1
%
%   This function requires IND2SUBVECT, also available (I hope) on the MathWorks 
%   File Exchange site.


numSets = length(varargin);
for i = 1:numSets
    thisSet = sort(varargin{i});
    if ~isequal(prod(size(thisSet)),length(thisSet))
        error('All inputs must be vectors.')
    end
    if ~isnumeric(thisSet)
        error('All inputs must be numeric.')
    end
    if ~isequal(thisSet,unique(thisSet))
        error(['Input set' ' ' num2str(i) ' ' 'contains duplicated elements.'])
    end
    sizeThisSet(i) = length(thisSet);
    varargin{i} = thisSet;
end

X = zeros(prod(sizeThisSet),numSets);
for i = 1:size(X,1)
    
    % Envision imaginary n-d array with dimension "sizeThisSet" ...
    % = length(varargin{1}) x length(varargin{2}) x ...
    
    ixVect = ind2subVect(sizeThisSet,i);
    
    for j = 1:numSets
        X(i,j) = varargin{j}(ixVect(j));
    end
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function has been downloaded from:
% https://www.mathworks.com/matlabcentral/fileexchange/5476-ind2subvect--multiple-subscript-vector-from-linear-index

function X = ind2subVect(siz,ndx)
%IND2SUBVECT Multiple subscripts from linear index.
%   IND2SUBVECT is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   X = IND2SUBVECT(SIZ,IND) returns the matrix X = [I J] containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%
%   For N-D arrays, X = IND2SUBVECT(SIZ,IND) returns matrix X = [I J K ...]
%   containing the equivalent N-D array subscripts equivalent to IND for 
%   an array of size SIZ.
%
%   See also IND2SUB.  (IND2SUBVECT makes a one-line change to IND2SUB so as
%   to return a vector of N indices rather than retuning N individual
%   variables.)%IND2SUBVECT Multiple subscripts from linear index.
%   IND2SUBVECT is used to determine the equivalent subscript values
%   corresponding to a given single index into an array.
%
%   X = IND2SUBVECT(SIZ,IND) returns the matrix X = [I J] containing the
%   equivalent row and column subscripts corresponding to the index
%   matrix IND for a matrix of size SIZ.  
%
%   For N-D arrays, X = IND2SUBVECT(SIZ,IND) returns matrix X = [I J K ...]
%   containing the equivalent N-D array subscripts equivalent to IND for 
%   an array of size SIZ.
%
%   See also IND2SUB.  (IND2SUBVECT makes a one-line change to IND2SUB so as
%   to return a vector of N indices rather than returning N individual
%   variables.)
 

% All MathWorks' code from IND2SUB, except as noted:

n = length(siz);
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
for i = n:-1:1
  X(i) = floor(ndx/k(i))+1;      % replaced "varargout{i}" with "X(i)"
  ndx = rem(ndx,k(i));
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function returns true if A is subset of B where A and B are binary
% vectors.
% A and B are binary vectors
function res=subsetCheck(A,B)
res = isequal(A, and(A, B));
return
end
