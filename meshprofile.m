function [ X, Y, Z, V ] = meshprofile( data, resolution )
%MESHPROFILE Interpolate values between profiles
%
%Syntax:  [ X, Y, Z, V ] = scatter3m( data_struct, resolution)
%
%Inputs: 
%    Required:
%        data Nx1 struct array for N profiles
%            .lat double: latitude of the profile
%            .lon double: longitude of the profile
%            .z Mx1 double array: every depth of the profile
%            .v Mx1 double array: value at the location
%    Optional:
%        resolution struct for the resolution of the generated plot
%                  .v double: vertical resolution (depth)
%                  .h double: horrizontal resolution (lat and lon)
%           default: 10
%
%Outputs: x Px1 double array of latitude
%         y Px1 double array of longitude
%         z Px1 double array of depth
%         v Px1 double array of value
%
%Example: Interpolate between profiles to get path of the cruise
% for i=size(s,2):-1:1;
%   data(i).lat = p(s(i)).lat;
%   data(i).lon = p(s(i)).lon;
%   data(i).z = p(s(i)).z;
%   data(i).v = p(s(i)).ct;
% end;
% res.v = 20; res.h = 30;
% [xf, yf, zf, vf] = meshprofile(data, res);
%
%m-files required: none
%
% Tested with: Matlab R2015a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 20th May 2015
% Last update: 20th May 2015    

% Check input
if nargin > 2
   error('Too many input arguments')
elseif nargin < 1
   error('Not enough input arguments')
end

% Check options
if nargin < 2
  resolution.h = 10;
  resolution.v = 10;
end;

X = []; Y = []; Z = []; V = [];
for l=2:size(data,2);
  % Get data
  a=l-1; b=l;
  n_a = size(data(a).z,1); n_b = size(data(b).z,1);
  x_a = data(a).lon; x_b = data(b).lon; % azimuth
  
  % Wrap 360 if necessary to cross 180 to -180
  if abs(x_a-x_b) > abs(wrapTo360(x_a)-wrapTo360(x_b)) % wrap to 360
    x_a = wrapTo360(x_a);
    x_b = wrapTo360(x_b);
  else % wrap to 180
    x_a = wrapTo180(x_a);
    x_b = wrapTo180(x_b);
  end
  
  % Build queries
  clear('x', 'y', 'z', 'v', 'd');
  %x(1:n_a) = x_a; x(n_a+1:n_a+n_b) = x_b;
  xq = wrapTo180(linspace(x_a, x_b, resolution.h)); 
  %y(1:n_a) = y_a; y(n_a+1:n_a+n_b) = y_b;
  yq = linspace(data(a).lat, data(b).lat, resolution.h);
  z(1:n_a+n_b) = [data(a).z;data(b).z];
  v(1:n_a+n_b) = [data(a).v;data(b).v];
  dist=1; %dist = sqrt((x_a-x_b)^2+(y_a-y_b)^2);
  d(1:n_a) = 0; d(n_a+1:n_a+n_b) = dist;
  
  % Interpolate in 3D
  [dq,zq] = meshgrid(linspace(0,dist,resolution.h), linspace(min(z),max(z),resolution.v));
  vq = griddata(d,z,v,dq,zq);
  
  % Rebuild in 4D
  k=1;
  xp=zeros(resolution.h*resolution.v,1);
  yp=zeros(resolution.h*resolution.v,1);
  zp=zeros(resolution.h*resolution.v,1);
  vp=zeros(resolution.h*resolution.v,1);
  for i=1:resolution.h;
    for j=1:resolution.v;
      xp(k)=xq(i);
      yp(k)=yq(i);
      zp(k)=zq(j);
      vp(k)=vq(j,i);
      k=k+1;
    end;
  end;

  % Remove NaN
  v2rm=find(isnan(vp));
  xp(v2rm)=[]; yp(v2rm)=[]; zp(v2rm)=[]; vp(v2rm)=[];

  % Save entire dataset
  X = [X;xp];
  Y = [Y;yp];
  Z = [Z;zp];
  V = [V;vp];
  
  % Plot source
  %hold('on');
  %scatter3(x, y, -z, 20, 'ko');
end;

end

