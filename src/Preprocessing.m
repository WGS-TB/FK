function tabl = Preprocessing( tabl )

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
% %% Can speed things up if list is not pre-purged.
tablbit=cull_rows_quick(tablbit);
num_rows=size(tablbit,1);
tabl = rfb(tablbit, num_cols);
end
