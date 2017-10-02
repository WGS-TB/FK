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