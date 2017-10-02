% This function generates DNF based on Matching M(v) method

% Input: v is the number of variables
% output: DNF form

function monomials = MatchingGraph( v )

if (mod(v, 2)==1)
    stop
end
monomials = [];
counter = 1;
for i=2:2:v
    temp = zeros(1,v);
    temp(i-1:i) = 1;
    monomials(counter,:) = temp;
    counter = counter+1;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function monomials = ThresholdGraph(v)
if (mod(v, 2)==1)
    stop % v must be even
end

monomials = [];
counter = 1 ;
for i=1:v
    for j=2:2:v
        if (i < j)
            temp = zeros(1,v);
            temp(i) = 1;
            temp(j) = 1;
            monomials(counter,:) = temp;
            counter = counter+1;
        end
    end
end

res = monomials;
return
end
