function oc_flag_info(flags)
%OC_FLAG_INFO display the name and description of every flag in flags
%   Information comes from flags from Level 2 and 3 of Ocean Color
%   http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.htm
%
%Syntax:  oc_flag_info(flags)
%
%Inputs: 
%    Required:
%        flags 1xN double containing l2_flags from Ocean Color (32 bit integer)
%
%
%Display: N lines with the name and description of every flag
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

% The table below shows the flags and masks that are operational
% in the Level 2 and Level 3 Ocean Color Processing.
% http://oceancolor.gsfc.nasa.gov/VALIDATION/flags.html
oc_flag(1).name='ATMFAIL'; oc_flag(1).description='Atmospheric correction failure';
oc_flag(2).name='LAND'; oc_flag(2).description='Pixel is over land';
oc_flag(3).name='PRODWARN'; oc_flag(3).description='One or more product warnings';
oc_flag(4).name='HIGLINT'; oc_flag(4).description='High sun glint';
oc_flag(5).name='HILT'; oc_flag(5).description='Observed radiance very high or saturated';
oc_flag(6).name='HISATZEN'; oc_flag(6).description='High sensor view zenith angle';
oc_flag(7).name='COASTZ'; oc_flag(7).description='Pixel is in shallow water';
oc_flag(8).name='spare'; oc_flag(8).description='spare';
oc_flag(9).name='STRAYLIGHT'; oc_flag(9).description='Straylight contamination is likely';
oc_flag(10).name='CLDICE'; oc_flag(10).description='Probable cloud or ice contamination';
oc_flag(11).name='COCCOLITH'; oc_flag(11).description='Coccolithofores detected';
oc_flag(12).name='TURBIDW'; oc_flag(12).description='Turbid water detected';
oc_flag(13).name='HISOLZEN'; oc_flag(13).description='High solar zenith';
oc_flag(14).name='spare'; oc_flag(14).description='spare';
oc_flag(15).name='LOWLW'; oc_flag(15).description='Very low water-leaving radiance (cloud shadow)';
oc_flag(16).name='CHLFAIL'; oc_flag(16).description='Derived product algorithm failure';
oc_flag(17).name='NAVWARN'; oc_flag(17).description='Navigation quality is reduced';
oc_flag(18).name='ABSAER'; oc_flag(18).description='possible absorbing aerosol (disabled)';
oc_flag(19).name='spare'; oc_flag(19).description='spare';
oc_flag(20).name='MAXAERITER'; oc_flag(20).description='Aerosol iterations exceeded max';
oc_flag(21).name='MODGLINT'; oc_flag(21).description='Moderate sun glint contamination';
oc_flag(22).name='CHLWARN'; oc_flag(22).description='Derived product quality is reduced';
oc_flag(23).name='ATMWARN'; oc_flag(23).description='Atmospheric correction is suspect';
oc_flag(24).name='spare'; oc_flag(24).description='spare';
oc_flag(25).name='SEAICE'; oc_flag(25).description='Possible sea ice contamination';
oc_flag(26).name='NAVFAIL'; oc_flag(26).description='Bad navigation';
oc_flag(27).name='FILTER'; oc_flag(27).description='Pixel rejected by user-defined filter';
oc_flag(28).name='SSTWARN'; oc_flag(28).description='SST quality is reduced';
oc_flag(29).name='SSTFAIL'; oc_flag(29).description='SST quality is bad';
oc_flag(30).name='HIPOL'; oc_flag(30).description='High degree of polarization';
oc_flag(31).name='PRODFAIL'; oc_flag(31).description='Derived product failure';
oc_flag(32).name='spare'; oc_flag(32).description='spare';

% Display flags name and description
for bit=flags;
  fprintf('%2d\t%10s\t%s\n', bit, oc_flag(bit).name, oc_flag(bit).description);
end;
  
end