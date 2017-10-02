function [tabl, total_call_size, max_call_size, processtime] = transversal(tabl)

%% An implementation of the hypergraph traversal algorithm
%%  described in Berge's textbook on hypergraphs.
%%
%% Input:
%%  tabl contains the hypergraph H.
%%
%% Output:
%%  tabl contains the dual (transversal) hypergraph H'.

%% By: Utz-Uwe Haus, Steffen Klamt and Tamon Stephen.
%% Copyright (C) 2007.
%%
%% This code was developed at the University of Magdeburg.
%%
%% Version of December 2007.
%% 
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  The GNU General Public License is available on-line at:
%%   <http://www.gnu.org/licenses/>.

%% Minimizes insertions by "keeping" rows that intersect e_i,
%%   and only inserting for edges that miss e_i completely.
%%   Insertion and culling is done in blocks via "ss_elim.m".
%%   Also reduces printed output as compared to first version.
%%
%% Finishes with the existing postprocessing code (not timed).

int_size=32;                    % Number of bits in an integer.

disp(' ');

%%% Keep track of time used via MATLAB cputime function (flaky).
%%% Time includes preprocessing and intermediate output, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preprocessing %%%%%%%%%%%%%%%%%%%%%%
%
% We use the existing preprocessing of S. Klamt from
% FluxAnalyzer and change only the main algorithm.
%

%%% Check for empty rows.
if(sum(any(tabl,2))<size(tabl,1));
  disp('Input error: empty row.');
  tabl=[];
  return;
end

tabl=~~tabl;                      % Make tabl binary.
num_rows=size(tabl,1);            % Number of rows in original matrix.
num_cols=size(tabl,2);            % Number of columns in original matrix.

%%% Flag 0 columns & remove.
%%% Remark: these are removed and never reinserted.
notinvolved=find(~any(tabl,1));
if(~isempty(notinvolved))
  tabl(:,notinvolved)=[];
  num_cols=num_cols-length(notinvolved);
end
tabl_cols=num_cols;              % Save for later.

idxdel=zeros(1,num_cols);         % Keep track of deleted column indices.

%%% Remove all 1's columns.  
zw=find(all(tabl,1));
anzess=length(zw);
idxdel(zw)=1;
all_ones_cols=zw;               % Save for later.

disp(' ');
disp(['Number of 1 element edges in transversal: ',num2str(anzess)]);
disp(' ');

%%% Flag equivalent columns so they will be ignored until preprocessing.
nums=0;
idxsub=1:num_cols;            % List of where to put deleted columns.
for i=1:num_cols-1
  if(idxdel(i)==0)
    i_col=tabl(:,i);
    for j=i+1:num_cols
      if (idxdel(j)==0)
        if (i_col == tabl(:,j))
          idxsub(j)=-i;
          idxdel(j)=1;
          nums=nums+1;
        end
      end
    end
  end
  zw = eq(-i*ones(1,num_cols),idxsub);
  num_eq = find(zw);
end
disp(['Result: ',num2str(nums),' equivalent columns found']);
% Delete equivalent columns.
zw=find(idxdel);
tabl(:,zw)=[];
num_cols=num_cols-length(zw);
disp([num2str(num_cols),' remaining elements']);

%%% Encode tabl as bits.
tablbit=rows_to_bits(tabl);
compact_cols=size(tablbit,2);     % Number of columns of reduced matrix.
% Find length of final column.
last_length=num_cols-int_size*(compact_cols-1);

%% Purge supersets rows.
%%
%% Can speed things up if list is not pre-purged.
tablbit=cull_rows_quick(tablbit);
num_rows=size(tablbit,1);

cputt=cputime;

%% Main part of Berge's algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize 
edge_list=uint32(zeros(0,compact_cols));   % Current list of edges.

v1=tablbit(1,:);
for i=1:compact_cols
  last=int_size;
  if ( i==compact_cols )
    last=last_length;
  end
  mask=uint32(1);
  for j=1:last
    v1j=bitand(mask,v1(i));
    if (v1j) % Generate edge consisting of a single vertex of v1.
      temp=uint32(zeros(1,compact_cols));
      temp(i)=mask;
      edge_list=vertcat(edge_list,temp);
    end
    mask=bitshift(mask,1);
  end
end
ell=size(edge_list,1);
% Now edge_list is Tr(H_1).

% Keep some statistics on lengths of edge lists encountered:
max_call_size = 0; % Longest total number of edges sent to merge routine.
total_call_size=0;

% Main routine.
pointer=1;
blank=uint32(zeros(1,compact_cols));
print_update=0;
for i=2:num_rows;   % For each i, compute Tr(H_i).
  print_update=print_update+1;
  vib=tablbit(i,:);
  num_verts=0;
  for j=1:compact_cols
    mask=uint32(1);
    for l=1:int_size
      if (bitget(vib,l))
        num_verts=num_verts+1;
      end
    end
  end
  new_edge_list=uint32(zeros(ell*num_verts,compact_cols));
  pointer_nel=0;    % Current location in new_edge_list.
  keepers=zeros(ell,1);    % Rows to retain as-is.
  for j=1:ell
    vjb=edge_list(j,:);
    if ( eq(bitand(vib,vjb),blank) ) % Is edge disjoiint?
      new_edges=uint32(zeros(num_verts,compact_cols)); % Yes, generate new.
      pointer=0;     % Current location in new_edges.
      for k=1:compact_cols
        mask=uint32(1);
        for l=1:int_size
          if (bitand(vib(k),mask))
            pointer=pointer+1;
            new_edge=vjb;
            new_edge(k)=bitset(vjb(k),l);
            new_edges(pointer,:)=new_edge;
          end
          mask=bitshift(mask,1);
        end
      end
      length_new=size(new_edges,1);
      new_edge_list(pointer_nel+1:pointer_nel+length_new,:)=new_edges(1:length_new,:);
      pointer_nel=pointer_nel+length_new;
    else
      keepers(j)=1; % No, can just keep this edge.
    end
  end
  new_edge_list=new_edge_list(1:pointer_nel,:);
  nell=pointer_nel;
  call_size=size(edge_list,1)+nell;
  total_call_size = total_call_size + call_size;
  if (call_size>max_call_size)
    max_call_size=call_size;
  end
  edge_list=ss_elim(edge_list(keepers==1,:),new_edge_list,num_cols);
  ell=size(edge_list,1);
  if (print_update*20>num_rows)
    print_update=0;
    disp(['After ',num2str(i),' (of ',num2str(num_rows),') iterations: ',num2str(ell),' edges.']);
  end
end
tabl=rfb(edge_list,num_cols);
  
disp(' ');
processtime = cputime-cputt;
disp(['Ready! Computation time: ',num2str(cputime-cputt),' sec']);

disp(['Longest list merge: ',num2str(max_call_size),'.']);
disp(['Compressed transversal length: ',num2str(ell+anzess),'.']);
disp(' ');

% disp('Postprocessing ...');

%% Post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reinserts empty columns
num_sets=size(tabl,1);                  % Total number of hypergraph edges.
tabl_compact=zeros(num_sets,tabl_cols);
zw=find(~idxdel);
tabl_compact(:,zw)=tabl;
tabl=tabl_compact;


num_sets=size(tabl,1);

%% Post-processing by S. Klamt.

if(num_sets>0)
	% disp('Expanding equivalent columns ...');
 	zw=find(idxsub<0);
 	for i=1:length(zw)
 		zw2=-idxsub(zw(i));
 		zw3=find(tabl(:,zw2));
 		zw4=length(zw3);
 		tabl(num_sets+1:num_sets+zw4,:)=tabl(zw3,:);
 		tabl(num_sets+1:num_sets+zw4,zw2)=0;
 		tabl(num_sets+1:num_sets+zw4,zw(i))=1;
 		num_sets=num_sets+zw4;
 	end
end
%% Insert size one cut sets.
if(anzess>0)
         essential = zeros(anzess,tabl_cols);
         for i=1:anzess
             essential(i,all_ones_cols(i))=1;
         end
         tabl=[essential; tabl];
end
num_sets=num_sets+anzess;
disp(['Final result: ',num2str(num_sets),' edges']);
disp(' ');

 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Culled]=ss_elim(A1,A0,nvars)

%% This takes two lists of binary-encoded rows, A0 and A1.
%% It assumes that there are no duplicate or superset rows in either list,
%% but that some rows in A0 may be supersets of rows in A1 (not vice-versa).
%% It then removes the rows from A0 that are supersets of rows in A1 and
%% returns the join of the two lists.
%% The vectors rs0 and rs1 are the row sums of A0 and A1 respectively.
%%
%% An additional parameter "nvars" is included, to avoid
%% having to find where the vectors end manually inside the routine.
%% This routine ignores bits past nvars.
%%
%% By TMS, June 2006.

int_size=32;          % Number of bits in an integer.
num_rows_0=size(A0,1);   % Number of rows of possibly superset matrices.
num_rows_1=size(A1,1);   % Number of rows base matrix.
ccols=size(A0,2);  % Number of compacted columns. % Should equal size(A1,2).

last_length=nvars-(ccols-1)*int_size; % Length of last integer in clause.
last_mask=uint32(2^last_length-1); % Clear garbage past end of vectors.

indic=true(num_rows_0,1); % Rows that should be kept.

if (num_rows_1<=num_rows_0)
  %% Look for supersets.
  for i=1:num_rows_1
    vi=A1(i,:);
    vi(ccols)=bitand(vi(ccols),last_mask);
    nss= (vi(1) == bitand(vi(1), A0(:,1))); % Rows not subsets of others.
    for j=2:ccols
      nss=nss & (vi(j) == bitand(vi(j), A0(:,j)));
    end
    indic=indic & ~nss;
  end
else
  last_mask=uint32(2^int_size-double(last_mask)-1);
  for i=1:num_rows_0
    vi=A0(i,:);
    vi(ccols)=bitor(vi(ccols),last_mask);
    nss = (A1(:,1) == bitand(vi(1), A1(:,1))); % Row not subsets of others.
    for j=2:ccols
      nss=nss & (A1(:,j) == bitand(vi(j), A1(:,j)));
    end
    if (~nss)
      % No row was a subset.
    else
      indic(i)=false;
    end
  end
end

Culled=[A1; A0(indic,:)]; % Build matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bits]=rows_to_bits(bytes)

%% Compresses the rows of a 0-1 matrix into an integer matrix
%% whose bits give the old matrix.  So, the first column of the
%% original matrix becomes the least significant bits of the
%% first column of the new matrix, and so on.
%%
%% By TMS, January 2006.

int_size=32;

num_rows=size(bytes,1);
num_cols=size(bytes,2);
compact_cols=ceil(num_cols/32);
last_length=num_cols-int_size*(compact_cols-1);

bits=uint32(zeros(num_rows,compact_cols));

onesvec=ones(size(bytes,1),1);

for j=1:compact_cols
  last=int_size;          % Last column to fill.
  if (j==compact_cols)
    last=last_length;
  end
  for k=1:last
    line=bytes(:,int_size*(j-1)+k);
    bits(:,j)=bitset(bits(:,j),k*onesvec,line);
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bytes]=rfb(bits, nvars)

%% Uncompresses the rows of a matrix encoded as bits into a 
%% 0-1 integer matrix.  So, the first column of the compressed 
%% matrix encodes the first 32 columns of the the new matrix, and so on.
%%
%% There is an additional parameter that encodes the number of
%% variables.  The original "rows_from_bits" autodetects the
%% length of the rows, but this takes too much time.
%% 
%% By TMS, April 2006.

int_size=32;                % Number of bits in an integer.
num_rows=size(bits,1);      % Number of rows of bit-encoded columns.
compact_cols=size(bits,2);  % Number of input columns.

last_length=nvars-(compact_cols-1)*int_size;

% Decode.
bytes=zeros(num_rows,nvars);
for i=1:compact_cols
  curr_col=bits(:,i);
  last_index=int_size;
  if (i == compact_cols)
    last_index=last_length;
  end
  mask=uint32(zeros(num_rows,1));
  mask=bitset(mask,1);
  for j=1:last_index
    bytes(:,int_size*(i-1)+j)=bitshift(bitand(curr_col,mask),1-j);
    mask=bitshift(mask,1);
  end
end   
return

