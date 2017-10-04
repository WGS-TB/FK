% The below code has been written based on the FK-B algorithm presented in "how to apply sat-solving for the equivalence test of mootone normal forms", page 5
function CA = MFK_B(cnf, dnf, DualGen, redundancychk, split_method)

CA = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                     Checking redundancy                 %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (redundancychk == 1)
    if(size(dnf,1) > 1)
        dnf = Irredundant(dnf);
    end
    if(size(cnf,1) > 1)
        cnf = Irredundant(cnf);
    end
end

nC = size(cnf,1); % Number of clauses in CNF
nD = size(dnf,1); % Number of monomials in DNF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                     Checking conditions                 %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(cnf,1) > 2 && size(dnf,1) > 2)
    [cnf, dnf, CA] = CommonVar_AllClauses( cnf, dnf );
    if (~isempty(CA))
        disp('Heuristic function works!!')
        return
%     else
%         fdsf=1;
%         cnf = cnf1;
%         dnf= dnf1;
    end
end

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
        CA = MFK_B([C_x_0; C_x_1], D_x_1, DualGen, redundancychk, split_method ); % line 8
        
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
                CA = MFK_B(C_cx_1, D_cx_0, DualGen, redundancychk, split_method); %line 11
            end
            
            if (~isempty(CA))
                CA(:, split_var) = 1; % adding split_var
                return %line 12
            end
        end %end of for (i in 1: length(C_x_0))
    else if mu_frequent_in_A(split_var, cnf, dnf) % if split variable is at most mu-frequent in CNF(line 13)
            CA = MFK_B(C_x_1, [D_x_0; D_x_1], DualGen, redundancychk, split_method); %line 14
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
                    CA = MFK_B(C_mx_0, D_mx_1, DualGen, redundancychk, split_method); %line 17
                end
                
                if (~isempty(CA))
                    CA(:, logical(vars))= 1;
                    return
                end
            end %end of for (i in 1: length(C_x_0))
        else % line 19
            
            CA = MFK_B([C_x_0; C_x_1], D_x_1, DualGen, redundancychk, split_method);  %line 20
            
            if(isempty(CA)) % line 21
                
                CA = MFK_B(C_x_1, [D_x_0; D_x_1], DualGen, redundancychk, split_method); % line 22
                
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


