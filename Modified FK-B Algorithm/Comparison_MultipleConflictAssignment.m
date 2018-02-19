
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% By: Nafiseh Sedaghat.
% nf_sedaghat@sfu.ca
% Copyright (C) 2018.
%
% This code was developed at Simon Fraser University.
%
% Version of February 2018.
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  The GNU General Public License is available on-line at:
%   <http://www.gnu.org/licenses/>.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
% This script compares FK_dualization and FK_dualization_M in terms of the
% number of calls to FK algorithm as well as the number of required
% iterations to compute dual of a given dnf
clear all

global call_counter
models = {'BIOMD0000000048','BIOMD0000000093','BIOMD0000000094',...
    'BIOMD0000000089', 'BIOMD0000000228', 'BIOMD0000000034', 'BIOMD0000000042',...
    'sms_33', 'ms_33','ac_200k', 'SDFP16', 'SDFP23'};

nn=1;
calls = zeros(numel(models), 2);
for nn=1:numel(models)
    model = models{nn};
    disp(strcat('nn is', num2str(nn), ' out of ', num2str(numel(models))))
    
    InputName = strcat(model, '_CNF_DNF.mat');
    
    load(InputName)
    dnf = Preprocessing( dnf );
    dnf = Irredundant(dnf);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng(6)
    call_counter = 0;
    [FK_cnf, FK_cpu, FK_len]= FK_Dualization( dnf, 1 );
    calls(nn, 1) = call_counter;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rng(6)
    call_counter = 0;
    [FKM_cnf, FKM_cpu, FKM_len]= FK_Dualization_M( dnf, 1 );
    calls(nn, 2) = call_counter;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                       Plot
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %++++++++++++++++++++++++++++++++++
    %  Comparing number of Iterations
    %++++++++++++++++++++++++++++++++++
    figure
    plot(FK_len(1:(find(FK_len==0, 1, 'first')-1)), 'r-o',...
        'MarkerSize',2,...
        'LineWidth', 2, ...
        'MarkerFaceColor','k');
    hold on;
    plot(FKM_len(1:(find(FKM_len==0, 1, 'first')-1)), 'b-o',...
        'LineWidth', 2, ...
        'MarkerSize',2,...
        'MarkerFaceColor','b');
    
    name = strrep(models(nn), '_', '-');
    title(char(name))
    
    xlabel('Iterations')
    ylabel('Length of CNF')
    xlim([1 inf])
    ylim([1 inf])
    hold off
    legend('FK','FK_M', 'Location','southeast')
    saveas(gcf, char(strcat(name, '.pdf')))
end


%++++++++++++++++++++++++++++++++++
%  Comparing number of FK calls
%++++++++++++++++++++++++++++++++++
 save('FKCalls_FKandFKM.mat', 'calls')

figure
b=bar(calls(:,[1, 2]));
b(1,1).FaceColor = 'b';
b(1,2).FaceColor = [1 171/255 0];
legend('FK', 'FK_M')

models(8) = strrep(models(8), '_', '-');
models(9) = strrep(models(9), '_', '-');
models(10) = strrep(models(10), '_', '-');

xticklabels(models)
xtickangle(55)
ylabel('Num of calls to FK algorithm')
name = strcat('FK_calls');

saveas(gcf, char(strcat(name, '.pdf')))

% close all

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function remove redundancy, i.e. supersets, in a given monotone Boolean function and returns an irredundant monotone Boolean function.
% MBF is a monotone Boolean function as a Binary matrix. The type of MBF, i.e. CNF or DNF, does not matter during the procedure of removing redundancy.
function MBF = Irredundant(MBF)

nrow = size(MBF, 1);
chk = false(1, nrow); % chk(i) is one if MBF(i,:) is a superset and should be removed

% Sort rows in ascending order based on the number of elements in each clause/monomial
summation = sum(MBF, 2);
[~,I] = sort(summation);
MBF = MBF(I, :);

for i=1:(nrow-1)
    a = MBF(i, :);
    for j=(i+1):nrow
        % This is the bottoleneck!!
        if ~chk(j) % if clause/monomial j is not already checked for being removed
            chk(j) = isequal(and(MBF(j,:), a), a);
        end
    end
end
MBF(chk,:) = [];

return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

% This function is part of the Berge algorithm for computing a hypergraph transversal for preprocessing.
% By TMS, June 2006.

function tabl = Preprocessing( tabl )

int_size=32;                    % Number of bits in an integer.

% disp(' ');


% Check for empty rows.
if(sum(any(tabl,2))<size(tabl,1));
%     disp('Input error: empty row.');
    tabl=[];
    return;
end

tabl=~~tabl;                      % Make tabl binary.
num_rows=size(tabl,1);            % Number of rows in original matrix.
num_cols=size(tabl,2);            % Number of columns in original matrix.

% Flag 0 columns & remove.
% Remark: these are removed and never reinserted.
notinvolved=find(~any(tabl,1));
if(~isempty(notinvolved))
    tabl(:,notinvolved)=[];
    num_cols=num_cols-length(notinvolved);
end
tabl_cols=num_cols;              % Save for later.

idxdel=zeros(1,num_cols);         % Keep track of deleted column indices.

% Remove all 1's columns.
zw=find(all(tabl,1));
anzess=length(zw);
idxdel(zw)=1;
all_ones_cols=zw;               % Save for later.

% disp(' ');
% disp(['Number of 1 element edges in transversal: ',num2str(anzess)]);
% disp(' ');

% Flag equivalent columns so they will be ignored until preprocessing.
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
% disp(['Result: ',num2str(nums),' equivalent columns found']);
% Delete equivalent columns.
zw=find(idxdel);
tabl(:,zw)=[];
num_cols=num_cols-length(zw);
% disp([num2str(num_cols),' remaining elements']);

% Encode tabl as bits.
tablbit=rows_to_bits(tabl);
compact_cols=size(tablbit,2);     % Number of columns of reduced matrix.
% Find length of final column.
last_length=num_cols-int_size*(compact_cols-1);

% Purge supersets rows.
%
% % Can speed things up if list is not pre-purged.
tablbit=cull_rows_quick(tablbit);
num_rows=size(tablbit,1);
tabl = rfb(tablbit, num_cols);
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function [Culled,length]=cull_rows_quick(A)

% This takes a list of rows of binary-encoded rows, and returns
% the rows from the list that are not strict supersets of any rows on
% the original list.
%
% Used for the Berge algorithm for computing a hypergraph transversal.
%
% By TMS, June 2006.

% Uses the faster setup found while coding the dual generator.
%
% Passes nvars rather that painfully computing it.
%

A=dup_del(A);

int_size=32;          % Number of bits in an integer.
num_rows=size(A,1);   % Number of rows to process.
ccols=size(A,2);  % Number of compacted columns.

indic=true(num_rows,1); % Rows that should be kept.

% Look for supersets.
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

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function [bits]=rows_to_bits(bytes)

% Compresses the rows of a 0-1 matrix into an integer matrix
% whose bits give the old matrix.  So, the first column of the
% original matrix becomes the least significant bits of the
% first column of the new matrix, and so on.
%
% Copyright Utz-Uwe Haus and Tamon Stephen, January 2006.

int_size=32;

num_rows=size(bytes,1);
num_cols=size(bytes,2);
compact_cols=ceil(num_cols/32);
last_length=num_cols-int_size*(compact_cols-1);

bits=zeros(num_rows,compact_cols,'uint32');

for j=1:compact_cols
    last=int_size;          % Last column to fill.
    if (j==compact_cols)
        last=last_length;
    end
    for k=1:last
        line=uint32(bytes(:,int_size*(j-1)+k));
        bits(:,j)=bitset(bits(:,j),k,line);
    end
end
return
end
%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function [bytes]=rfb(bits, nvars)

% Uncompresses the rows of a matrix encoded as bits into a
% 0-1 integer matrix.  So, the first column of the compressed
% matrix encodes the first 32 columns of the the new matrix, and so on.
%
% There is an additional parameter that encodes the number of
% variables.  The original "rows_from_bits" autodetects the
% length of the rows, but this takes too much time.
%
% Copyright Utz-Uwe Haus and Tamon Stephen, April 2006.

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
    mask=zeros(num_rows,1,'uint32');
    mask=bitset(mask,1);
    for j=1:last_index
        bytes(:,int_size*(i-1)+j)=bitshift(bitand(curr_col,mask),1-j);
        mask=bitshift(mask,1);
    end
end

return
end

%###########################################################################
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%###########################################################################

function [Culled,length]=dup_del(A)

% This takes a list of rows of binary-encoded columns, and deletes
% duplicate rows by sorting the list.
%
% By TMS, January 2006.

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

% Sort rows.
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