%   Detailed explanation goes here
% This function checks whether split variable satisfy the below condition:
% |{m \in A: x \in m}|/|A| <= 1/mu(|A| . |B|)
% which we say that x is at most frequent in A
function res = mu_frequent_in_A( split_var, A, B )

ss = sum(A, 1);
len = size(A, 1);
left_side = ss(split_var) / len;
right_side = 1/mu_function(len * size(B,1));

res = left_side <= right_side;
return

end

