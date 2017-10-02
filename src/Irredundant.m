% MBF is a monotone Boolean function
function MBF = Irredundant(MBF)

nrow = size(MBF, 1);
chk = false(1, nrow); % chk(i) is one if MBF(i,:) is a superset and should be removed

% Sort rows in ascending order based on the number of elements in each clause/monomial
summation = sum(MBF, 2);
[~,I] = sort(summation);
MBF = MBF(I, :);

for i=1:nrow
    a = MBF(i, :);
    for j=(i+1):nrow
        if ~chk(j) % if clause/monomial j is not already checked for being removed
            chk(j) = isequal(and(MBF(j,:), a), a);
        end
    end
end

MBF(chk,:) = [];
return
end


