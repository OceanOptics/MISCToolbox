function [ need_qc, sun_elevation ] = need_npqc( dt, lat, lon, varargin )
%NEED_NPQC Determine if the profile need a quenching correction
%   Base on a model of the sun elevation
%
% Inputs: 
%    Required:
%        dt double date and time UTC in matlab datenum format
%        lat double containing latitude
%        lon double containing longitude (wrapTo180 preferred)
%    Optional:
%        min_sun_elevation double of the minimum sun elevation to [0 90]
% Outputs:need_qc boolean
%           true need a npq correction
%           false don't need a npq correction
%         sun_elevation
%
% Examples:
% % Know if need npq correction based on sun elevation
% dt = datenum(datetime());
% lat = 44.89;
% lon = -68.66;
% if need_npqc(dt, lat, lon)
%   disp('Need NPQ correction');
% end;
%
% Required:
%   sun_position from Reda et Andreas (2003)
%
% Tested with: Matlab R2015a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 26th May 2015
% Last update: 26th May 2015
%

% Check input
if nargin > 4
   error('Too many input arguments')
elseif nargin < 1
   error('Not enough input arguments')
end

% Default arguments
min_sun_elevation = 5; % deg

% Get others arguments
for i=1:nargin-3;
  if isscalar(varargin{i})
    min_sun_elevation = varargin{i};
  else
    error('Unknown argument format');
  end;
end;

% Process
sun_elevation = sunElevation(dt, lat, lon);
if min_sun_elevation < sun_elevation
  need_qc = true;
else
  need_qc = false;
end;

end

function sun_elevation = sunElevation(dt, lat, lon)
  loc.latitude = lat;
  loc.longitude = lon;
  loc.altitude =  0;
  dt = datevec(dt); tt.UTC = 0;
  tt.year = dt(1); tt.month = dt(2); tt.day = dt(3);
  tt.hour = dt(4); tt.min = dt(5); tt.sec = dt(6);
  sun = sun_position(tt, loc);
  sun_elevation = 90-sun.zenith;
end