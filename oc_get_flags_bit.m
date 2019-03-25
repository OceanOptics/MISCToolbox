function flags_bit = oc_get_bit_flags(flags, disp)
%GET_OC_BIT_FLAGS return the true bit number from flags
%   can display flags from Level 2 and 3 of Ocean Color
%   http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html
%
%Syntax:  flags_bit = get_oc_bit_flags( flags )
%         flags_bit = get_oc_bit_flags( flags, true )
%
%Inputs:
%    Required:
%        flags <double> l2_flags from Ocean Color (32 bit integer)
%
%    Optional:
%        disp <boolean> display informations about flags if set to true
%           default: false
%
%Outputs: flags_bit Nx1 array of integer containing the bit number of all the flags
%
% Tested with Matlab R2015a, 2016a, and 2016b
%
% Require: oc_flag_info
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 12 August 2015
% Last update: 13 August 2015
% Acknowledgements to : hfeng and anton.korosov
%   http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?tid=1261

% Check input
if nargin > 2
   error('Too many input arguments')
elseif nargin < 1
   error('Not enough input arguments')
end
if nargin < 2
  disp = false;
end

% Get every flag with string method
% n_bit = 32;
% flags_bin = dec2bin(flags,n_bit);
% flags_bit = n_bit+1 - strfind(flags_bin, '1');

% Get every flag with binary method
foo = uint8(bitget(int32(flags),1:32));
flags_bit = find(1 == foo);

% Display flags
if disp
  % Display flag binary code
  %fprintf('%d\n', flags);
  fprintf('%s\n', mat2str(bitget(int32(flags),32:-1:1)));
  % Display flags name and description
  oc_flag_info(flags_bit);
end;

end