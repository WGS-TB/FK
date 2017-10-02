% This function returns true if A is subset of B where A and B are binary
% vectors.
% A and B are binary vectors
function res=subsetCheck(A,B)
res = isequal(A, and(A, B));
return
end
