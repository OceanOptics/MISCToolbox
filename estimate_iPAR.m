function iPARq = estimate_iPAR(dt, lat, lon, iPAR, dtq)
% ESTIMATE_IPARNOON estimate iPAR at noon based on sun_position
%
% Input:
%   dt: double containing the date in vector [YYYY MM DD HH MM SS] in UTC 
%   lat: double latitude (in degrees, north of equator is positive)
%   lon: double longitude (in degrees, positive for east of Greenwich)
%   iPAR: double Instantaneous Photosynthetically Available Radiation
%   dtq: double of the date in vector in UTC of query iPAR
% 
% Output:
%   iPARq: double iPAR at noon in solar time
%
% Require:
%   sun_position
%
% Nils Haentjens- University of Maine
% nils.haentjens@maine.edu
% 1st July 2015

s = size(dt);
if s(1) == 6 && s(2) ~= 6;
  dt = dt';
elseif s(2) ~= 6;
  error('dt should have at list one size of 6');
end;


tt.UTC = 0;
location.altitude = 0;

for i=1:size(dt,1);
  % Actual sun elevation
  tt.year = dt(i,1);
  tt.month = dt(i,2);
  tt.day = dt(i,3);
  tt.hour = dt(i,4);
  tt.min = dt(i,5);
  tt.sec = dt(i,6);
  location.latitude = lat(i);
  location.longitude = lon(i);
  sun = sun_position(tt, location);
  elevation(i) = 90-sun.zenith*(sun.zenith<90);
  
  % Query sun elevation
  tt.year = dtq(i,1);
  tt.month = dtq(i,2);
  tt.day = dtq(i,3);
  tt.hour = dtq(i,4);
  tt.min = dtq(i,5);
  tt.sec = dtq(i,6);
  sun = sun_position(tt, location);
  elevationq(i) = 90-sun.zenith*(sun.zenith<90);
  
  % Estimate query iPAR
  iPARq(i) = iPAR * elevationq / elevation;
end;