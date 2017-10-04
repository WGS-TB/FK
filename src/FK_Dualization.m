function [cnf, cput, cnf_len]= FK_Dualization( dnf )

cputt=cputime;
n_clause = 1;
cnf = Initialize_CNF('random', dnf, n_clause);
cnf = Irredundant(cnf);
cnf_len = zeros(1, 1000);
while_counter = 0;
iter = 0;
while(while_counter <= 400)
    
    %     disp(['The length of CNF is ', num2str(size(cnf, 1))])
    iter = iter + 1;
    cnf_len(iter) = size(cnf, 1);
    %     if (size(cnf, 1)== 26)
    %         jio=0;
    %     end
    CA = FK_B(cnf, dnf, 1, 1, 'mostFreq');
    if (size(CA,1)>1)
        dsada=11;
    end
    if ~isempty(CA)
        CA(CA < 0) = 0;
        if (sum(CA)==0)
            stop
        end
    end
    
    if isempty(CA)
        display('Assign is NULL', '\n')
        break
    end
    
    new_clause = [];
    MFP = Maximum_False_Point(find(CA), dnf);
    temp = zeros(1,size(dnf,2));
    temp(setdiff(1:size(dnf,2), MFP)) = 1;
    new_clause = [new_clause; temp];
    
    cnf = [cnf; new_clause];
    
    if size(cnf, 1) ~= size(Irredundant(cnf),1)
        sfs=1;
    end
    while_counter = while_counter + 1;
end%end while
cput = cputime - cputt;
return
end

