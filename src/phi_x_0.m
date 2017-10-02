% phi_x_0 denotes the formula that consists of terms of phi from which x is removed
% MBF: Monotone Boolean Function
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
