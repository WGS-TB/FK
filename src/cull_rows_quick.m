function [Culled,length]=cull_rows_quick(A)

%% This takes a list of rows of binary-encoded rows, and returns 
%% the rows from the list that are not strict supersets of any rows on 
%% the original list.
%%
%% Used for the Berge algorithm for computing a hypergraph transversal.
%%
%% By TMS, June 2006.

%% Uses the faster setup found while coding the dual generator.
%%
%% Passes nvars rather that painfully computing it.
%%

A=dup_del(A);

int_size=32;          % Number of bits in an integer.
num_rows=size(A,1);   % Number of rows to process.
ccols=size(A,2);  % Number of compacted columns.

indic=true(num_rows,1); % Rows that should be kept.

%% Look for supersets.
for i=1:num_rows
  vi=A(i,:);
  nss= (vi(1) == bitand(vi(1), A(:,1)));
  for j=2:ccols
    nss=nss & (vi(j) == bitand(vi(j), A(:,j)));
  end
  nss(i)=false;
  indic=indic & ~nss;
end

Culled=A(indic,:);
return
end