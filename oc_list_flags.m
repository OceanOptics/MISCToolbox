function list_flags = oc_list_flags(flags)
%OC_LIST_FLAGS return a list of unique flag bit in matrix flags
%
%Syntax:  list_flags = oc_list_flags(flags)
%
%Inputs: 
%    Required:
%        flags NxM double containing l2_flags from Ocean Color (32 bit integer)
%
%
%Outputs: list_f Lx1 of int32 unique flag in the matrix
%
% Tested with Matlab R2015a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 17 August 2015
% Last update: 17 August 2015

% Check input
if nargin > 1;
   error('Too many input arguments')
elseif nargin < 1;
   error('Not enough input arguments')
end

list_flags = [];
for i=1:size(flags,1);
  for j=1:size(flags,2);
    foo = uint8(bitget(int32(flags(i,j)),1:32));
    list_flags = unique([list_flags find(1 == foo)]);
  end;
end;

end
    