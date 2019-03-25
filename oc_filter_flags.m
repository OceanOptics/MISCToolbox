function [val_f, nan_sel] = oc_filter_flags(val, flags, filter)
%OC_FILTER_FLAGS return the value (val_f) replacing flagged data with NaN
%
%Syntax:  val_f = oc_filter_flags( val, flags, filter )
%
%Inputs:
%    Required:
%        val   <NxM double> values (chlorophyll or par for example)
%        flags <NxM double> l2_flags from Ocean Color (32 bit integer)
%        filters <1XL unsigned integer between 1 and 32> list of flags to remove
%
%
%Outputs: val_f <Nx1 double> of val with NaN values for filtered content
%         nan_sel <NxM logical> of flags matching the filter
%
% Tested with Matlab R2015a, 2016a, and 2016b
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 17 August 2015
% Last update: 17 August 2015

% Check input
if nargin > 3;
   error('Too many input arguments')
elseif nargin < 3;
   error('Not enough input arguments')
end
if size(val) ~= size(flags);
  error('val and flags should be same size');
end;
if size(filter,1) ~= 1;
  error('filter should be 1xM');
end;


% Initialize output
%val_f = NaN(size(val));
val_f = val;

% Get values to remove
f = zeros([size(flags) size(filter,2)]);
for k=1:size(filter,2);
  f(:,:,k) = bitget(int32(flags),filter(k));
end;

% Remove values
nan_sel = logical(sum(f,3));
val_f(nan_sel) = NaN;

end