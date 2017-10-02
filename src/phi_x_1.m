% phi_x_1 denotes the formula that consists of all terms of phi that do not contain x.
% MBF: Monotone Boolean Function
function res = phi_x_1(MBF, var)

[rows , ~]= find(MBF(:,var)==1);
rows = unique(rows);
MBF = MBF(setdiff(1:size(MBF,1), rows),:); % Only those rows that does not contain var

% If a clause is empty of any variable we remove that clause
rem = sum(MBF, 2)==0;
if (sum(rem)>0)
    MBF(rem, :) = [];
end

res = MBF;

end
