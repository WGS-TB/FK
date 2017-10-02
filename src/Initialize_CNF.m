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

