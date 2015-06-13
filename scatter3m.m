function [ hs, hc ] = scatter3m( x, y, z, s, c, varargin )
%SCATTER3M 3D spherical scatter plot
%
%Syntax:  [ hs, hc ] = scatter3m(x, y, z, s, c, 'withCoastLine', 'withLLLines')
%
%Inputs: 
%    Required:
%        x Nx1 double array corresponding to the longitude/azimuth
%        y Nx1 double array corresponding to the latitude/elevation
%        z Nx1 double array corresponding to the depth/radius
%        s Nx1 double array or double specify the size of the circles
%        c Nx1 double array or double specify the color of the circles
%    Optional:
%        'WithCoastLine'
%        'WithLLLines'
%        zlim 1X2 double array specify the limit of the Z axis (earth rotation axis)
%        radius specify sphere radius (default:6371000 m earth radius)
%
%Outputs: hs lines handles of the scatter plot
%         hc lines handles of the coast lines
%
%Example: Plot Profiles in the Artic
% figure(3); clf(3, 'reset');
% scatter3m(xf, yf, zf*1000, 60, vf, 'withCoastLine', 'withLLLines', [60 90]);
% view(90,40);
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
if nargin > 9
   error('Too many input arguments')
elseif nargin < 5
   error('Not enough input arguments')
end

% Check options
radius = earthRadius;
withCoastLine = false;
withLLLines = false;
withLimits = false;
% Identify varargin
for i=1:nargin-5;
  if ischar(varargin{i});
    if strcmp(varargin{i}, 'withCoastLine')
      withCoastLine = true;
    elseif strcmp(varargin{i}, 'withLLLines')
      withLLLines = true;
    else
      warning('Unknown parmeter %s', varargin{i});
    end;
  elseif ismatrix(varargin{i});
    withLimits = true;
    lat_lim = varargin{i};
  elseif isscalar(varargin{i});
    radius = varargin{i};
  else
    error('Unknown argument format');
  end;
end;

hold('on');
% Longitude and latitude lines Plot
if withLLLines
  lat_lines = [-90:10:90];
  lon_lines = [-180:20:170];
  res = 100;
  for i=1:size(lat_lines,2);
    [latl, lonl, el] = sph2cart(degtorad(linspace(-180,180,res)), degtorad(linspace(lat_lines(i),lat_lines(i),res)), radius);
    plot3(latl, lonl, el, 'Color', [0.7 0.7 0.7]);
  end;
  for i=1:size(lon_lines,2);
    [latl, lonl, el] = sph2cart(degtorad(linspace(lon_lines(i),lon_lines(i),res)), degtorad(linspace(-90,90,res)), radius);
    plot3(latl, lonl, el, 'Color', [0.7 0.7 0.7]);
  end;
end;

% Coast Line Plot
if withCoastLine
  coast = load('coast');
  [latc, lonc, ec] = sph2cart(degtorad(coast.long), degtorad(coast.lat), radius);
  hc = plot3(latc, lonc, ec);
end;

% Scatter Plot
[xs, ys, zs] = sph2cart(degtorad(x), degtorad(y), radius-z);
hs = scatter3(xs, ys, zs, s, c, 'filled', 'o');
hold('off');

% Design features
axis('equal', 'off', 'tight');
colormap('jet');

% Set Limits
if withLimits;
  flag=true;
  if size(lat_lim,2) < 4
    lat_lim(3) = 0;
    lat_lim(4) = 360;
    flag=false;
  end;
  [x1, y1, z1] = sph2cart(degtorad(lat_lim(3)), degtorad(lat_lim(1)), radius-min(z));
  [x2, y2, z2] = sph2cart(degtorad(lat_lim(4)), degtorad(lat_lim(2)), radius);
  x_min = min(x1, x2); y_min = min(y1, y2); z_min = min(z1, z2);
  x_max = max(x1, x2); y_max = max(y1, y2); z_max = max(z1, z2);
  if flag
    xlim([x_min  x_max]);
    ylim([y_min  y_max]);
  end;
  zlim([z_min  z_max]);
end;

end

