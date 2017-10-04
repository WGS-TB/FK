function [cnf, cput, cnf_len]= MFK_Dualization( dnf )

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
    CA = MFK_B(cnf, dnf, 1, 1, 'mostFreq');
    
    CA(CA < 0) = 0;
    CA(sum(CA, 2)==0,:) = [];
    
    
    if isempty(CA)
        display('Assign is NULL', '\n')
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
    
%     if size(cnf, 1) ~= size(Irredundant(cnf),1)
%         sfs=1;
%     end
    while_counter = while_counter + 1;
end%end while
cput = cputime - cputt;
return
end

