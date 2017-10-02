function maximal_vars = Maximum_False_Point( vars, dnf )
% finding minimal sets that can be zero, in this way we can have maximal set of 1

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

