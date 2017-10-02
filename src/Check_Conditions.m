% Checking conditions:
% 1: for each c in CNF and m in DNF => intersect(m, c) != null
% 2: DNF and CNF must have exactly the same variables
% 3: max{|m|: m \in DNF} <= |C|   &   max{|c|: c \in CNF} <= |D|

% Inputs: CNF and DNF as binary matrices

function CA = Check_Conditions(cnf, dnf)

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
        [~, vars] = find(dnf==1);
        temp = zeros(1, size(dnf, 2));
        temp(vars) = -1; % False vars
        CA = temp;
    end
    return;
end

%++++++++++++++++++++++++++++
%       First Condition
%++++++++++++++++++++++++++++

%chk = cnf * dnf' > 0;
chk = zeros(size(cnf, 1), size(dnf, 1));
for i=1:size(cnf, 1)
    for j=1:size(dnf,1)
        chk(i,j) = sum(and(cnf(i,:), dnf(j,:)))>0;
    end
end

if(sum(chk==0)>0)
    [~, y] = find(chk==0);
    CA = dnf(y(1), :);
    return
end


%++++++++++++++++++++++++++++
%       Second Condition
%++++++++++++++++++++++++++++

[~, temp] = find(cnf);
cnf_vars = unique(temp);

[~, temp] = find(dnf);
dnf_vars = unique(temp);

check = isequal(cnf_vars, dnf_vars);

if (~check)
    x = setdiff(dnf_vars, cnf_vars);
    % if (length(setdiff(cnf.vars, dnf.vars))>0 && length(setdiff(dnf.vars, cnf.vars))>0)
    %   stop("Error!! Both CNF and DNF have a variable that does not appear in the other one. Here is Check Conditions function!")
    if (numel(x)>0) % the extra variables are in dnf
        
        [~, temp] = find(dnf);
        dnf_vars = unique(temp);
        
        
        x = setdiff(dnf_vars, cnf_vars);
        
        [r,~] = find(dnf(:,x));
        
        %         CA = setdiff(find(dnf(r(1),:)), x);
        CA = dnf(r(1),:);
        CA(x) = 0;
        return
        
    else % the extra variables are in cnf
        
        x = setdiff(cnf_vars, dnf_vars);
        [~, temp] = find(cnf);
        cnf_vars = unique(temp);
        
        [r,~] = find(cnf(:,x));
        %         CA = union(setdiff(cnf_vars, find(cnf(r(1),:))), x);
        temp = ones(1, numel(cnf_vars)) - cnf(r(1),:);
        temp(:,x) = 1;
        CA = temp;
        
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
            S(longest_dnf) = str2num(bin);
            
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
    
    if (t2==FALSE)
        [~, I] = max(sum(cnf, 2));
        longest_cnf = cnf(I,:);
        
        for i=1:((2^sum(longest_cnf))-1) % We are seeking for proper subset of longest monomial in dnf
            
            bin = dec2bin(i, sum(longest_cnf));
            S = zeros(1, size(cnf, 2));
            S(longest_cnf) = str2num(bin);
            
            chk = zeros(1, size(dnf,1));
            for j=1:size(dnf,1)
                chk(j) = sum(and(S,dnf(j,:)))>0;
            end
            if(~ismember(0, chk))
                CA = ones(1, numel(cnf_vars))- S; % we must set *proper subset* of the largest monoials to true
                % browser()
                return
            end
        end
        % stop('There is sth wrong in CHeck_Conditions' )
    end
    
end

% If we reach here, it means that all of conditions are satisfied.
if (~isempty(CA))
    CA = sort(CA);
end
return
end
