function list_flags = oc_list_flags(flags, bits)
%OC_LIST_FLAGS return a list of unique flag bit in matrix flags
%
%Syntax:  list_flags = oc_list_flags(flags)
%
%Inputs:
%    Required:
%        flags <NxM double> l2_flags from Ocean Color (32 bit integer)
%    Optional:
%        bits <integer of 16 for SST or 32 for Ocean Color>
%            default: 32
%
%Outputs: list_f Lx1 of int32 unique flag in the matrix
%
% Tested with Matlab R2015a, 2016a, and 2016b
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 17 August 2015
% Last update: 17 August 2015

% Check input
if nargin > 2;
   error('Too many input arguments')
elseif nargin < 1;
   error('Not enough input arguments')
end
if nargin < 2;
  bits = 32;
end;

list_flags = [];
for i=1:size(flags,1);
  for j=1:size(flags,2);
    if bits == 16;
      foo = uint8(bitget(int16(flags(i,j)),1:16));
    elseif bits == 32;
      foo = uint8(bitget(int32(flags(i,j)),1:32));
    else
      error('Unknown value for bits, should be 16 or 32');
    end;
    list_flags = unique([list_flags find(1 == foo)]);
  end;
end;

end
