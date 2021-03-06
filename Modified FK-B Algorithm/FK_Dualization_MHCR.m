
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
% This function computes dual of a given dnf using hash table and finding mutiple conflict assignments while it reduces the number of redundancy check.
% In this function, if the number of variables in CNF is less than 7, we compute its canonical form and save the corresponding permutation.
% Then, we permute dnf using the stored permutation and make the key.
% Also, prior to inserting the conflict assignment in the hash table, we permute it as well.
%
% For more information please see "Speeding up the Fredman-Khachiyan Dualization Algorithm B" paper.
% It starts the cnf by randomly finding some clauses.
%
% Inputs:
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% n_clause: The number of clauses to initialize cnf
% thresh: An integer number that is used as a threshold to use hash table or not. The hash table is used when the number of clauses and monomials are larger than 'thresh'.
% thresh = 0 means that we use hash table for every cnf and dnf without considering their size.
%
% Outputs:
% cnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% cput: It measures cputime of running the code
% cnf_len: It is an array that its i^th entry indicates to the number of clauses in iteration i.
%
% Please not that currently cput is deactivated; To use it you can uncomment 'tic' and 'cput=toc;' instructions.
%
% Example:
%
% dnf = [1 0 0 0 1 0 1 1 0; 0 1 0 1 0 1 1 0 0; 0 0 1 1 1 0 0 0 1];
% n_clause = 1;
% thresh = 0;
% [cnf, cput, cnf_len] = FK_Dualization_MHCR( dnf, n_clause, thresh );
%
% Just to see if the returned cnf is correct:
% FK(cnf, dnf, 'mostFreq') % it returns [] that means the returned cnf and dnf are dual.

function [cnf, cput, cnf_len]= FK_Dualization_MHCR( dnf, n_clause, thresh )

% tic
%+++++++++++++++++++++++++++++++
% 		Initializing CNF
%+++++++++++++++++++++++++++++++

cnf = Initialize_CNF('random', dnf, n_clause);
cnf = Irredundant(cnf);

%++++++++++++++++++++++++++++++++++++++
% 		Initializing the hash table
%++++++++++++++++++++++++++++++++++++++

StartHashing(thresh)

%+++++++++++++++++++++++++++++++
% 		Main procedure
%+++++++++++++++++++++++++++++++

cnf_len_est = 3000;  % An estimation of the number of clasues in cnf
cnf_len = zeros(1, cnf_len_est);
while_counter = 0;
iter = 0;
while(while_counter <= cnf_len_est)
    
    if mod(while_counter,100)==0
        disp(['The length of CNF is ', num2str(size(cnf, 1))]);
    end
    iter = iter + 1;
    cnf_len(iter) = size(cnf, 1);
    CA = FK_MHCR(cnf, dnf, 'mostFreq');
    CA(CA < 0) = 0;
    CA(sum(CA, 2)==0,:) = [];
    
    
    if isempty(CA)
        display('CA is NULL', '\n')
        break
    end
    
    CA = Irredundant(CA);
    new_clause = [];
    for k=1:size(CA,1)
        MFP = Maximum_False_Point(find(CA(k,:)), dnf);
        temp = zeros(1,size(dnf,2));
        temp(setdiff(1:size(dnf,2), MFP)) = 1;
        new_clause = [new_clause; temp];
    end
    new_clause = unique(new_clause, 'rows');
    cnf = [cnf; new_clause];
    while_counter = while_counter + 1;
end%end while

% cput=toc;
cput=0;
return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function randomly finds 'n_clauses' clauses included in cnf based on the given dnf.

% Inputs:
% method: currently it can be only 'randomly'; It is for future extensions.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% n_clause: The number of clauses to initialize cnf

% Output:
% cnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
%		The number of rows in cnf is equal to 'n_clause'

function cnf = Initialize_CNF( method, dnf, n_clause )

if(strcmp(method,'random'))
    
    counter = 1;
    cnf = [];
    while(size(cnf, 1) < n_clause)
        
        minimal_vars = [];
        tst = dnf;
        
        while (sum(sum(tst))>0)
            ss = sum(tst, 1);
            I = datasample(find(ss>0), 1, 'Replace',false);
            minimal_vars = union(minimal_vars, I);
            tst(logical(tst(:,I)),:) = [];
            %             tst(:,I) = 0;
        end
        minimal_vars = Minimality_Check(minimal_vars, dnf);
        new = zeros(1, size(dnf, 2));
        new(minimal_vars) = 1;
        cnf = [cnf; new];
        cnf = Irredundant(cnf);
        counter = counter + 1;
    end
    
end
return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This script checks wheather 'vars' is minimal in the given monotone boolean function (MBF). Its output is the minimal version.

function minimal_vars = Minimality_Check( vars, MBF )

mat = MBF(:, vars); % keep only columns that are in x
remList = ones(1, numel(vars));
for i=size(mat,2):-1:1
    temp = mat;
    temp(:,i)=[];
    
    chk = sum(temp, 2) > 0;
    if sum(chk) == size(temp, 1)
        mat(:,i) = 0;
        remList(i) = 0;
    end
end
minimal_vars = vars(logical(remList));
return
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

% This function starts hash tables and the related variables.
% Since they are global variables, at first we clear them to make sure everything is clear in the begining.

% Input:
% hash_threshold: An integer number that is used as a threshold to use hash table or not. The hash table is used when the number of clauses and monomials are larger than 'hash_threshold'.
% hash_threshold = 0 means that we use hash table for every cnf and dnf without considering their size.
function StartHashing(hash_threshold)

clear global cellHash cellSize repchk Hash_thresh popularKey
clear global popuInd canonic_cntr hash_cntr bigPermMat
global cellHash cellSize Hash_thresh popularKey popu_ind bigPermMat
global hash_cntr canonic_cntr canonic_fetch

cellHash = cell(10000,2);
popularKey = zeros(10000, 1); popu_ind = 1;
cellSize = 1;

canonic_fetch=0; hash_cntr=0; canonic_cntr=0;
Hash_thresh = hash_threshold;

% bigPermMat is used for finding canonical form.
for i=1:6
    if i>1
        bigPermMat{i} = Perm_BinaryMatGen(i);
    end
end

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function gets the number of variables, n, and makes a big binary
% matrix corresponding to the whole possible permutations for an MBF with n
% variables.
% It works for small 'n'.

function bigMat = Perm_BinaryMatGen(n)

% Generating all permutations for n variables
% perms is of dimension [n!, 2^n]
perms = permlist(n);

% the number of possible permutations
nrow = factorial(n);
ncol = 2^n;

bigMat = zeros(nrow*(ncol), ncol);
mat = zeros(ncol, ncol);

for i=1:nrow
    p = perms(i, :);
    mat(sub2ind(size(mat), 1:ncol, p)) = 1;
    bigMat(((i-1)*ncol)+1:(i*ncol),:) = mat;
end

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function finds minimal sets that can be zero, in this way we can have maximal set of 1

function maximal_vars = Maximum_False_Point( vars, dnf )

minimal_vars = zeros(1,size(dnf,2));
tst = dnf;
tst(:, vars) = 0;
tst(sum(tst,2)==0,:)=[];
while (sum(sum(tst))>0)
    ss = sum(tst, 1);
    [~, ind] = max(ss);
    minimal_vars(find(minimal_vars==0, 1, 'first')) = ind(1);
    tst(logical(tst(:,ind(1))),:) = [];
end
minimal_vars(minimal_vars==0)=[];
minimal_vars = Minimality_Check(minimal_vars, dnf);
maximal_vars = setdiff(1:size(dnf,2), minimal_vars);
return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% The below code has been written based on the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5
% This function is FK_MHR computing canonical form (if applicable) to store keys.
% The script returns multiple conflict assignments if there exist
% This script uses hash table to memorize cnf and dnf and their conflict assignment(s)

% In this function, if the number of variables in CNF is less than 7, we compute its canonical form and save the corresponding permutation.
% Then, we permute dnf using the stored permutation and make the key.
% Also, prior to inserting the conflict assignment in the hash table, we permute it as well.
% In this version, cell arrays are used as hash table and use the original version of CNF and DNF as key in the hash table;
% We use hash when the number of clauses/monomials are more than 'Hash_thresh'.
% This script reduces the number of redundancy checks.

% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.
% The number of columns in cnf and dnf must be the same.
% split_method: the method of choosing split variable; the only option is 'mostFreq' and it is for further improvements.

% Output:
% CA: It is the conflict assignment(s) between the given cnf and dnf; each row corresponds to a conflict assignment. The variables corresponding to '1' in CA are considered as True and the other variables that might be '0' or '-1' are considered as False.
% If the given cnf and dnf are equivalent, CA would be [] (NULL).

function CA = FK_MHCR(cnf, dnf, split_method)

global cellHash cellSize repchk Hash_thresh


thresh = Hash_thresh;
% if cellSize == 669
%     fsd=0;
% end
% Find_Duplicates(cellHash, cellSize)
CA = [];
key1 = []; key2 = []; key3 = []; key4 = []; key5 = []; key6 = [];
repchk = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                     Checking redundancy                 %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        tempCNF = [C_x_0; C_x_1];
        
%         if(size(D_x_1,1) > 1)
%             D_x_1 = Irredundant(D_x_1);
%         end
        if(size(tempCNF,1) > 1)
            tempCNF = Irredundant(tempCNF);
        end
        
        cond = (size(tempCNF,1) > thresh && size(D_x_1,1) > thresh);
        
        if(cond)
            [CA, key1, perm1] = Hashing_CanonicalForm(tempCNF, D_x_1);
        end
        if (~cond || strcmp(CA, 'None'))
            CA = FK_MHCR(tempCNF, D_x_1, split_method ); % line 8
            if(~isempty(key1))
                if (~isempty(CA))
                    CA_mapped = CA(:,perm1);
                else
                    CA_mapped = CA;
                end
                cellHash(cellSize,:) = {key1, CA_mapped};
                cellSize = cellSize + 1;
                key1 = [];
            end
        end
        if (~isempty(CA))
            return %line 9
        end
        
        for i=1: size(C_x_0,1) %line 10
            key2 = [];
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
                
                if(size(D_cx_0,1) > 1)
                    D_cx_0 = Irredundant(D_cx_0);
                end
                if(size(C_cx_1,1) > 1)
                    C_cx_1 = Irredundant(C_cx_1);
                end
                
                cond = (size(C_cx_1,1) > thresh && size(D_cx_0,1) > thresh);
                if(cond)
                    [CA, key2, perm2] = Hashing_CanonicalForm(C_cx_1, D_cx_0);
                end
                if (~cond || strcmp(CA, 'None'))
                    CA = FK_MHCR(C_cx_1, D_cx_0, split_method); %line 11
                    
                    if( ~isempty(key2))
                        if (~isempty(CA))
                            CA_mapped = CA(:,perm2);
                        else
                            CA_mapped = CA;
                        end
                        cellHash(cellSize,:) = {key2, CA_mapped};
                        cellSize = cellSize + 1;
%                         key2 = [];
                    end
                end
            end
            if (~isempty(CA))
                CA(:, split_var) = 1; % adding split_var
                return %line 12
            end
        end %end of for (i in 1: length(C_x_0))
    else if mu_frequent_in_A(split_var, cnf, dnf) % if split variable is at most mu-frequent in CNF(line 13)
            
            tempDNF = [D_x_0; D_x_1];
            
            if(size(tempDNF,1) > 1)
                tempDNF = Irredundant(tempDNF);
            end
%             if(size(C_x_1,1) > 1)
%                 C_x_1 = Irredundant(C_x_1);
%             end
            
            cond = (size(C_x_1,1) > thresh && size(tempDNF,1) > thresh);
            if(cond)
                [CA, key3, perm3] = Hashing_CanonicalForm(C_x_1, tempDNF);
            end
            if (~cond || strcmp(CA, 'None'))
                CA = FK_MHCR(C_x_1, tempDNF, split_method); %line 14
                
                if(~isempty(key3))
                    if (~isempty(CA))
                        CA_mapped = CA(:,perm3);
                    else
                        CA_mapped = CA;
                    end
                    cellHash(cellSize,:) = {key3, CA_mapped};
                    cellSize = cellSize + 1;
                    key3 = [];
                end
            end
            
            if (~isempty(CA))
                CA(:, split_var) = 1;
                return %line 15
            end
            
            for i=1: size(D_x_0,1) %line 16
                key4 = [];
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
                    
                    
                    if(size(D_mx_1,1) > 1)
                        D_mx_1 = Irredundant(D_mx_1);
                    end
                    if(size(C_mx_0,1) > 1)
                        C_mx_0 = Irredundant(C_mx_0);
                    end
                    
                    cond = (size(C_mx_0, 1) > thresh && size(D_mx_1, 1) > thresh);
                    if(cond)
                        [CA, key4, perm4] = Hashing_CanonicalForm(C_mx_0, D_mx_1);
                    end
                    if (~cond || strcmp(CA, 'None'))
                        CA = FK_MHCR(C_mx_0, D_mx_1, split_method); %line 17
                        
                        if(~isempty(key4))
                            if (~isempty(CA))
                                CA_mapped = CA(:,perm4);
                            else
                                CA_mapped = CA;
                            end
                            cellHash(cellSize,:) = {key4, CA_mapped};
                            cellSize = cellSize + 1;
%                             key4 = [];
                        end
                    end
                end
                
                if (~isempty(CA))
                    CA(:, logical(vars))= 1;
                    return
                end
            end %end of for (i in 1: length(C_x_0))
        else % line 19
            tempCNF = [C_x_0; C_x_1];
            
%             if(size(D_x_1,1) > 1)
%                 D_x_1 = Irredundant(D_x_1);
%             end
            if(size(tempCNF,1) > 1)
                tempCNF = Irredundant(tempCNF);
            end
            
            cond = (size(tempCNF,1) > thresh && size(D_x_1,1) > thresh);
            if(cond)
                [CA, key5, perm5] = Hashing_CanonicalForm(tempCNF, D_x_1);
            end
            if (~cond || strcmp(CA, 'None'))
                CA = FK_MHCR(tempCNF, D_x_1, split_method);  %line 20
                
                if(~isempty(key5))
                    if (~isempty(CA))
                        CA_mapped = CA(:,perm5);
                    else
                        CA_mapped = CA;
                    end
                    cellHash(cellSize,:) = {key5, CA_mapped};
                    cellSize = cellSize + 1;
                    key5 = [];
                end
            end
            
            if(isempty(CA)) % line 21
                tempDNF = [D_x_0; D_x_1];
                
                if(size(tempDNF,1) > 1)
                    tempDNF = Irredundant(tempDNF);
                end
%                 if(size(C_x_1,1) > 1)
%                     C_x_1 = Irredundant(C_x_1);
%                 end
                
                cond = (size(C_x_1,1) > thresh && size(tempDNF,1) > thresh);
                if(cond)
                    [CA, key6, perm6] = Hashing_CanonicalForm(C_x_1, tempDNF);
                end
                if (~cond || strcmp(CA, 'None'))
                    CA = FK_MHCR(C_x_1, tempDNF, split_method); % line 22
                    
                    if(~isempty(key6))
                        if (~isempty(CA))
                            CA_mapped = CA(:,perm6);
                        else
                            CA_mapped = CA;
                        end
                        
                        cellHash(cellSize,:) = {key6, CA_mapped};
                        cellSize = cellSize + 1;
                        key6 = [];
                    end
                end
                
                if (~isempty(CA))
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

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function gets a cnf and dnf and make key from them
% then it checks if the key is already in the hash table
% If the key is already there, it returns CA, otherwise it returns 'None'.

function [CA, key, varperm] = Hashing_CanonicalForm(cnf, dnf)
global cellHash cellSize popularKey popu_ind
global canonicalKey canonic_fetch hash_cntr canonic_cntr

canonicalKey  = 0;
[key, varperm] = MakeKey_CanonicalForm(cnf, dnf);

cmp = strcmp(key, cellHash(1:(cellSize-1),1));
if (any(cmp))
    popularKey(popu_ind) = find(cmp == 1, 1, 'first');
    popu_ind = popu_ind + 1;
    
    if (canonicalKey==1)
        canonic_fetch = canonic_fetch + 1;
    end
    r = cellHash(cmp,2);
    r = r{1,1};
    if ismatrix(r)
        CA = r;
        if(~isempty(CA))
            inverse_permutation = sortrows([varperm', (sort(varperm))'], 1);
            inverse_permutation = inverse_permutation(:,2);
            CA = CA(:, inverse_permutation);
            canonic_cntr = canonic_cntr + 1;
        end
    else
        inverse_permutation = sortrows([varperm', (sort(varperm))'], 1);
        inverse_permutation = inverse_permutation(:,2);
        
        CA = cell2mat(r{1,1});
        CA = CA(:, inverse_permutation);
        canonic_cntr = canonic_cntr + 1;
    end
    hash_cntr = hash_cntr + 1;
else
    CA = 'None';
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function makes canonical form based on "COUNTING INEQUIVALENT MONOTONE BOOLEAN FUNCTIONS" paper.

% Inputs:
% cnf: a binary matrix showing a monotone Boolean function in conjuctive normal form. Each row corresponds to a clause in CNF and each column corresponds to a variable.
% dnf: a binary matrix showing a monotone Boolean function in disjunctive normal form. Each row corresponds to a monomial in DNF and each column corresponds to a variable.

%Outputs:
% key: the key to be used in hash table
% varperm: the permutation that converts cnf to its canonical form
function [key, varperm] = MakeKey_CanonicalForm(cnf, dnf)
global canonicalKey
% n = size(cnf, 2);
cnf_perm = 1:size(cnf, 2);

delimiter = size(cnf,2);
active_vars = sum(sum(cnf,1)>0); % active variables in cnf
if(~isempty(cnf) && active_vars > 1 && active_vars < 7)
    [cnf, cnf_perm] = Find_CanonicalForm(cnf);
    %mat_temp = dnf;
    dnf(:, sort(cnf_perm)) = dnf(:, cnf_perm);
    canonicalKey = 1;
end

if(size(cnf,1)>1)
    lcnf = reshape(cnf, 1, []);
else
    lcnf = cnf;
end
% cnf_chr = strjoin(string(lcnf));
% cnf_chr = int2str(lcnf);


% reshaping DNF as a vector
if(size(dnf,1)>1)
    ldnf = reshape(dnf, 1, []);
else
    ldnf = dnf;
end

%key = strcat('CNF=[', cnf_chr, '];DNF=[', dnf_chr, ']');
comb = [lcnf, delimiter, ldnf];
key = int2str(comb);
key = char(key);

varperm = 1:size(cnf,2);
varperm(sort(cnf_perm)) = cnf_perm;
return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function gets an MBF as a binary matrix and returns its canonical form as well as the corresponding permutation
% bigPermMat is a big binary matrix including whole possible permutations

function [res, varperm]= Find_CanonicalForm(MBF)

global bigPermMat

ind = sum(MBF,1)>0; % present variables in MBF
origg_MBF = MBF;
origg_vars = find(ind);
MBF = MBF(:, ind);
nvar = sum(ind); % active variables
n = 2^nvar;
nperm = factorial(nvar);

% bigPermMat = Perm_BinaryMatGen(nvar);

% Convert MBF to a binary vector
% vec is of size [n x 1]

vec = mts_to_mbf(MBF', nvar);

%+++++++++++++++++++++++++
%+++++++++++++++++++++++++

bigRes = 1*((bigPermMat{nvar} * vec') > 0); % bigRes is of dimension [2^nvar, 1]
decvec = zeros(1, nperm);

for i=1:nperm
    binvec = bigRes(((i-1)*n)+1:(i*n));
    decvec(i) = bi2de(binvec'); %convert binary vector to a decimal number
end

[~, ind] = min(decvec);

% Canonical is a binary vector
canonical = bigRes(((ind-1)*n)+1:(ind*n));
mapping = bigPermMat{nvar}(((ind-1)*n)+1:(ind*n),:);

[perm, ~] = find(mapping);

% the rows in truth table are sorted in lexicographical descending order
% because of the defalt considered order in Tamon's code
truthTable = sortrows(getcondvects(nvar), nvar:-1:1,'descend');
ind = find(sum(truthTable, 2)==1);

permuted_truthTable = truthTable(perm,:);
permuted_ind = find(sum(permuted_truthTable, 2)==1);

perm = perm(ind);

[~, old_var] = find(truthTable(perm,:)==1);

[r, new_var] = find(permuted_truthTable(perm,:)==1);
temp = sortrows([r, new_var]);
new_var = temp(:,2);

permutation = sortrows([old_var, new_var]);
permutation = permutation(:, 2)';
permutation = origg_vars(permutation);

varperm = 1:size(origg_MBF,2);
varperm(origg_vars) = permutation;
% Convert canonical to an MBF
mapped_MBF = mts(canonical', nvar)';
res = origg_MBF;
res(:,origg_vars) = mapped_MBF;
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function result = mts_to_mbf(M, nvars)

% result = mts_to_mbf(M, nvars)
%
% Given a matrix of minimal true statements, outputs the MBF f with these minterms.
%
% Author: TJ Yusun, Simon Fraser University
% Copyright 2011
%
%  This file is part of InequivalentMBF.
%
%  InequivalentMBF is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  InequivalentMBF is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%	long with InequivalentMBF. If not, see <http://www.gnu.org/licenses/>.


comp = mat_2n(nvars);
result = zeros(1, 2^nvars);

if M == 0
    result = ones(1, 2^nvars);
    return
end

total = size(M,2);


for i = 1:total
    
    curr_col = M(:,i);
    ind = find(curr_col);
    card = size(ind);
    
    temp = ones(1, 2^nvars);
    
    for k = 1:card
        temp = bitand(temp,comp(ind(k),:));
    end
    
    result = bitor(result,temp);
    
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function [condvects] = getcondvects(i)
% GETCONDVECTS returns a matrix of binary condition vectors.
% GETCONDVECTS(I) returns a matrix of all possible binary condition
% vectors for a logical system with (I) inputs.
% INPUT: (I) shall be an integer >= 1
% OUTPUT: CONDVECTS is a binary matrix of size [2^I,I]
% METHOD: The method uses three nested FOR loops that work to target
% only those cells whose value should be true (which of course is 50%
% of the matrix).  Therefore, the algorithm scales as [0.5*i*2^i]
% which is optimal.  Furthermore, the algorithm has low memory
% footprint and overhead.
% Copyright 2011, Paul Metcalf
% Acknowledgements: James Tursa and Nico Schl�mer

g = 2;
i2 = 2^i;
condvects = false(i2,i);
for m = 1 : 1 : i
    m2 = 2^m;
    m3 = (m2/2)-1;
    i3 = i-m+1;
    for g = g : m2 : i2
        for k = 0 : 1 : m3
            condvects(g+k,i3) = true;
        end
    end
    g = m2+1;
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function result = mts(f, nvars)

% result = mts(f, nvars)
%
% Given a MBF f in truth table form, this function outputs its minimal true statements.
%
% Author: TJ Yusun, Simon Fraser University
% Copyright 2011
%
%  This file is part of InequivalentMBF.
%
%  InequivalentMBF is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  InequivalentMBF is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%	long with InequivalentMBF. If not, see <http://www.gnu.org/licenses/>.
%

comp = mat_2n(nvars);
result = zeros(nvars,0);

if f == 1
    result = zeros(nvars,1);
    return
end

while find(f) ~= 0
    
    x = ones(1, 2^nvars);
    
    curr_bit = find(f, 1, 'last');
    col = comp(:,curr_bit);
    result = [result col];
    
    ind = find(col);
    card = size(ind);
    
    for k = 1:card
        x = bitand(x,comp(ind(k),:));
    end
    
    f = max(f - x,0);
    
end
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function A = permlist(nvars)

% A = permlist(nvars)
%
% Generates all permutations of n-variable BFs
% looking at the truth table forms (so length 2^n instead of n)
%
% Author: TJ Yusun, Simon Fraser University
% Copyright 2011
%
%  This file is part of InequivalentMBF.
%
%  InequivalentMBF is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  InequivalentMBF is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%	long with InequivalentMBF. If not, see <http://www.gnu.org/licenses/>.
%

% First, generate all permutations in S^n
list = perms(1:nvars);
A = [];

% Then permute the index vector 1:2^n corresponding to
% all these permutations to get resulting
% indices for each
for i = 1:size(list,1)
    A = [A; findperm([1:2^nvars],list(i,:),nvars)];
end

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function fperm = findperm(f,sgn,nvars)

% fperm = findperm(f,sgn,nvars)
% 
% INPUTS
% nvars: number of variables
% f: function of length 2^nvars
% sgn: permutation in S^n
%
% Permutes the boolean function f (given in 0-1 form, with length 2^nvars)
% according to sgn, a permutation of nvars variables.
%
%
% Author: TJ Yusun, Simon Fraser University
% Copyright 2011
%
%  This file is part of InequivalentMBF.
%
%  InequivalentMBF is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  InequivalentMBF is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%	long with InequivalentMBF. If not, see <http://www.gnu.org/licenses/>.
%


% Check dimensions
if length(f) == 2^nvars && length(sgn) == nvars
else
    disp('One or more inputs are invalid.')
    return
end

mat = mat_2n(nvars);
mat_new = perm(mat, sgn, 1);

nperm = ones(1,2^nvars);

for k = 1:2^nvars
   col_temp = mat_new(:,k);
   
   mask = ones(1,2^nvars);
   
   for i = 1:nvars
      temp = or(and(col_temp(i),mat(i,:)),~or(col_temp(i),mat(i,:)));
      mask = and(mask,temp);
   end
    
    nperm(k) = find(mask,1);
    
end

fperm = perm(f,nperm,2);

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function mat = mat_2n(nvars)

% Outputs matrix with n rows, where the ith row contains 2^(n - i + 1) 
% blocks of length 2^(i - 1), alternating between all 1's and all 0's. 
% Each column corresponds to a truth assignment where a 1 in row i means 
% that variable x_i is set to 1.
%
% Essentially the (0,1) Rademacher matrix of order n
%
% Author: TJ Yusun, Simon Fraser University
% Copyright 2011
%
%  This file is part of InequivalentMBF.
%
%  InequivalentMBF is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  InequivalentMBF is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%	long with InequivalentMBF. If not, see <http://www.gnu.org/licenses/>.
%

if nvars == 0
    mat = [];
elseif nvars == 1
    mat = [1 0];
else
    mat_temp = mat_2n(nvars - 1);
    
    mat = [mat_temp mat_temp];
    last_row = [ones(1,2^(nvars-1)) zeros(1,2^(nvars - 1))];
    mat = [mat; last_row];
end

end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function result = perm(mat,sgn,dim)

% result = perm(mat,sgn,dim)
%
% Permutes the rows or columns of the matrix mat according to the
% permutation sgn, a row vector containing the numbers from 1 to n.
% Input dim is the dimension being permuted.
%
% Author: TJ Yusun, Simon Fraser University
% Copyright 2011
%
%  This file is part of InequivalentMBF.
%
%  InequivalentMBF is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  InequivalentMBF is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%	long with InequivalentMBF. If not, see <http://www.gnu.org/licenses/>.
%

% Check for correct size, too.

if dim == 1    
    
    if size(mat,1) == length(sgn)
    else
        disp('Incorrect dimensions')
        return
    end

    result = mat(sgn,:);

elseif dim == 2

    if size(mat,2) == length(sgn)
    else
        disp('Incorrect dimensions')
        return
    end
    result = mat(:,sgn);
else
    disp('Invalid argument for dimension. dim should be 1 or 2')
    return
end
end