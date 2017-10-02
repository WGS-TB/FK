function [Culled,length]=dup_del(A)

%% This takes a list of rows of binary-encoded columns, and deletes
%% duplicate rows by sorting the list.
%%
%% By TMS, January 2006.

int_size=32;          % Number of bits in an integer.
num_rows=size(A,1);   % Number of rows to process.
compact_cols=size(A,2);  % Number of compacted columns.

% First, find length of final column.
last_length=0;
last_col=A(:,compact_cols);
mask=uint32(ones(num_rows,1));

for i=1:int_size
  if ( ~bitand(mask,last_col) )
    % Empty column
  else
    last_length=i;
  end
  mask=bitshift(mask,1);
end
  
indic=ones(num_rows, 1);
A=sortrows(A);

for i=2:num_rows
  if( A(i-1,:)==A(i,:) )
    indic(i) = 0;  % Redundant.
  end
end
A=A(indic==1,:);
num_rows=size(A,1);

%% Sort rows.
num_vert=zeros(num_rows,1);
for i=1:compact_cols
  mask=uint32(ones(num_rows,1));
  last=int_size;
  if (i==compact_cols)
    last=last_length;
  end
  for j=1:last
    num_vert=num_vert+ne(bitand(mask,A(:,i)),0);
    mask=bitshift(mask,1);
  end
end
[blah sortorder]=sort(num_vert);
Culled=A(sortorder,:);

Culled=Culled(1:num_rows,:);
length=num_rows;

return
end