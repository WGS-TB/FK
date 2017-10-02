% This function checks if there is some elements common in all clauses, 
% it exists as a single term in the other boolean function
% If it does not pass the test, it returns the conflict assignment

function [cnf, dnf, CA] = CommonVar_AllClauses( cnf, dnf )
 
CA = [];
  % Find single terms in cnf and dnf
  ind = sum(cnf, 2)==1;
  [~, cnf_singleVars] = find(cnf(ind,:)==1);
  
  ind = sum(dnf, 2)==1;
  [~, dnf_singleVars] = find(cnf(ind,:)==1);

  % Find variables that are present in all clauses/monomials
  cnf_all = find(sum(cnf, 1) == size(cnf,1));
  dnf_all = find(sum(dnf, 1) == size(dnf,1));

   comm1 = intersect(dnf_all, cnf_singleVars);
   comm2 = intersect(cnf_all, dnf_singleVars);
   
   vars = union(comm1, comm2);
   cnf(:, vars) = 0;
   dnf(:, vars) = 0;

   a = setdiff(cnf_all, dnf_singleVars);
   if ~isempty(a)
       CA = zero(1, size(cnf,2));
       CA(a) = 1;
       return
   end

   b = setdiff(dnf_all, cnf_singleVars);
   if ~isempty(b)
       false_vars = b;
       [~,dnf_vars] = find(dnf==1);
       dnf_vars = unique(dnf_vars);
       
       CA = zero(1, size(dnf,2));
       CA(setdiff(dnf_vars, false_vars)) = 1;
       return
   end
   return
end

