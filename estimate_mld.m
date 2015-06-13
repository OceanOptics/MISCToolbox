function [ mld, mld_i ] = estimate_mld( z, varargin )
%ESTIMATE_MLD Estimates mixer layer depth (MLD)
%   The MLD is estimated directly from profiles given in input
%   Various criterion can be used to determine the MLD this function use by
%   default the fixed criterion in density of 0.03kg/m^3 difference from
%   surface. A fixed criterion in temperature (0.2 C) or a variable
%   criterion in density (which enables the definition of BLT) can also be
%   used.
%
% Inputs: 
%    Required:
%        z Nx1 double array of depth profile in m
%    Optional:
%        rho Nx1 double array of in-situ density in kg/m^3
%        t Nx1 double array of Conservative Temperature (ITS-90) in deg C
%        s Nx1 double array of Absolute Salinity in g/kg
%        method string containing the method to use:
%           fixed_temperature   threshold
%           fixed_density   threshold (default)
%           variable_density   threshold
%           fixed_density_gradient
%           average (average of value obtain by 4 previous methods)
%           robust  (median of value obtain by 4 previous methods)
%        criterion double containing the criterion
%           0.2 degre C by default if fixed_temperature method
%           0.03 kg/m^3 by default if fixed_density method
%
% Outputs:mld double containing the estimate MLD 
%         mld_i integer containing the index of the estimate of the MLD
%                 for the sorted depth with shallow depth at top
%
% Examples:
% % Get MLD from a density profile rho
% mld = estimate_mld(rho);
% % or
% mld = estimate_mld(rho, 'fixed_density');
%
% % Get MLD from a temperature profile t
% mld = estimate_mld(t, 'fixed_temperature');
%
% % Get MLD with a custom criterion custom_criterion
% mld = estimate_mld(t, 'fixed_temperature', custom_criterion);
% % or
% mld = estimate_mld(rho, custom_criterion);
%
% % Get MLD with a variable density threshold criterion
% mld = estimate_mld(z, rho, theta, sa, 'variable_density');
%
% % Get MLD with a fixed density gradient criterion
% mld = estimate_mld(z, rho, 'fixed_density_gradient');
%
% % Get MLD with the average of all the methods
% mld = estimate_mld(z, rho, theta, sa, 'average');
%
% % Get MLD with the median of all the methods
% mld = estimate_mld(z, rho, theta, sa, 'robust');
%
% Required:
%   Gibbs Seawater Oceanographic Toolbox (GSW)
%     gsw_rho(si(i),thetai(j),0);
%
% Tested with: Matlab R2015a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 26th May 2015
% Last update: 26th May 2015
%
% References:
%     Zawada, D. G., Zaneveld, J. R. V, Boss, E., Gardner, W. D.,
%       Richardson, M. J., & Mishonov, A. V. (2005). A comparison of
%       hydrographically and optically derived mixed layer depths.
%       Journal of Geophysical Research: Oceans, 110(11), 1?13.
%       http://doi.org/10.1029/2004JC002417
%     http://www.ifremer.fr/cerweb/deboyer/mld/Surface_Mixed_Layer_Depth.php

% Check input
if nargin > 6
   error('Too many input arguments')
elseif nargin < 1
   error('Not enough input arguments')
end

% Check options
% Default arguments
method = 'fixed_density';
rho_criterion = 0.03; % kg/m^3
theta_criterion = 0.2; % degre C
density_gradient_criterion =  0.01; % kg/m^3
z_ref = 10; % m % reference depth

% Get method from arguments
for i=1:nargin-1;
  if ischar(varargin{i});
    if strcmp(varargin{i}, 'fixed_temperature') ||...
        strcmp(varargin{i}, 'fixed_density') ||...
        strcmp(varargin{i}, 'variable_density') ||...
        strcmp(varargin{i}, 'fixed_density_gradient') ||...
        strcmp(varargin{i}, 'average') ||...
        strcmp(varargin{i}, 'robust');
      method = varargin{i};
    else
      warning('Unknown method %s', varargin{i});
    end;
  end;
end;

% Load other arguments
i=1;
while i <= nargin-1;
  if ismatrix(varargin{i}) && ~ischar(varargin{i}) && ~isscalar(varargin{i});
    if strcmp(method, 'fixed_density')
      rho = varargin{i};
    elseif strcmp(method, 'fixed_temperature')
      theta = varargin{i};
    elseif strcmp(method, 'variable_density')  || strcmp(method, 'average') || strcmp(method, 'robust')
      rho = varargin{i};
      theta = varargin{i+1};
      sa = varargin{i+2};
      i=i+2;
    elseif strcmp(method, 'fixed_density_gradient')
      rho = varargin{i};
    else
      error('Unknown method');
    end;
  elseif isscalar(varargin{i});
    if strcmp(method, 'fixed_density')
      rho_criterion = varargin{i};
    elseif strcmp(method, 'fixed_temperature') || strcmp(method, 'variable_density')
      theta_criterion = varargin{i};
    elseif strcmp(method, 'fixed_density_gradient')
      density_gradient_criterion = varargin{i};
    elseif  strcmp(method, 'average') || strcmp(method, 'robust')
      warning('Default criterion will be used for every method');
    else
      error('Unknown method');
    end;
  elseif ~ischar(varargin{i})
    error('Unknown argument format');
  end;
  i=i+1;
end;

%% Format Profiles
% Set orientation to Nx1 array if necessary
s = size(z);
if s(2) > s(1)
  z = z';
end;
if exist('rho', 'var');
  s = size(rho);
  if s(2) > s(1)
    rho = rho';
  end;
end;
if exist('theta', 'var');
  s = size(theta);
  if s(2) > s(1)
    theta = theta';
  end;
end;
if exist('sa', 'var');
  s = size(theta);
  if s(2) > s(1)
    sa = sa';
  end;
end;

% Set order according to depth to be sure interp1 works properly
% Shallow first and deeper last
[z, z_i] = sort(z);
if exist('rho', 'var');
  rho = rho(z_i);
end;
if exist('theta', 'var');
  theta = theta(z_i);
end;
if exist('sa', 'var');
  sa = sa(z_i);
end;

%% Find MLD according to the selected method
switch method
  case 'fixed_temperature'
    % Kara Isothermal Layer Depth (ILD)
    % Find temperature at reference depth (theta_ref)
    theta_ref = interp1(z, theta, z_ref,'linear','extrap');
    % Find depth where theta = theta_ref +- theta_criterion
    [i i]=min(abs(z-z_ref)); % Start at index of z_ref
    while i < size(theta,1) && abs(theta(i) - theta_ref) < theta_criterion;
      i = i + 1;
    end;
    mld = interp1([theta(i-1) theta(i)], [z(i-1) z(i)], theta_ref-theta_criterion); % +
    if isnan(mld)
      mld = interp1([theta(i-1) theta(i)], [z(i-1) z(i)], theta_ref+theta_criterion); % -
    end;
    mld_i = i;
  case 'fixed_density'
    % Levitus Density Difference (Lev)
    % Find temperature at reference depth (theta_ref)
    rho_ref = interp1(z, rho, z_ref,'linear','extrap');
    % Find depth by interp 
%     mld = interp1(rho, z, rho_criterion + rho_ref); % might find it way to low
%     mld_i = min(abs(z-mld));
    % Find depth where rho = rho_ref + rho_criterion by invrementation
    [i i]=min(abs(z-z_ref)); % Start at index of z_ref
    while i < size(rho,1) && rho(i) < rho_criterion + rho_ref;
      i = i + 1;
    end;
    if i > 1
      mld= interp1([rho(i-1) rho(i)], [z(i-1) z(i)], rho_ref+rho_criterion);
    else
      warning('utils:misc_lab:estimate_mld:fixed_density', 'MLD is first value available of profile');
      mld = z(1);
    end;
    mld_i = i;
  case 'variable_density'
    % Kara Density Difference (Kara)
    % Interpolate reference values
    sa_ref = interp1(z, sa, z_ref,'linear','extrap');
    theta_ref = interp1(z, theta, z_ref,'linear','extrap');
    rho_ref = interp1(z, rho, z_ref,'linear','extrap');
    [i i]=min(abs(z-z_ref)); % Start at index of z_ref
    % Compute variable density criterion threshold
    delta_rho = gsw_rho(sa_ref,theta_ref-theta_criterion,0) - gsw_rho(sa_ref,theta_ref,0);
    % Find MLD
    while i < size(rho,1) && rho(i) < delta_rho + rho_ref;
      i = i + 1;
    end;
    if i > 1
      mld= interp1([rho(i-1) rho(i)], [z(i-1) z(i)], rho_ref+delta_rho);
    else
      warning('utils:misc_lab:estimate_mld:variable_density', 'MLD is first value available of profile');
      mld = z(1);
    end;
    mld_i = i;
  case 'fixed_density_gradient'
    % Lukas-Lindstorm Fixed Density Gradient (fxGrad)
    % Interpolate rho every step (2 m)
    step = 2;
    rho2(:,1) = interp1(z, rho, z_ref:2:z(end),'linear','extrap');
    
    i=1;
    while i < size(rho2,1) && rho2(i) < density_gradient_criterion + rho2(1);
      i = i + 1;
    end;
    mld = z_ref + step * i;
    [mld_i mld_i] = min(abs(z-mld)); % Should be the shallower and not the minimum
  case 'average'
    mld(4) = estimate_mld(z, theta, 'fixed_temperature');
    mld(3) = estimate_mld(z, rho, 'fixed_density');
    mld(2) = estimate_mld(z, rho, theta, sa, 'variable_density');
    mld(1) = estimate_mld(z, rho, 'fixed_density_gradient');
    mld = mean(mld);
    [mld_i mld_i] = min(abs(z-mld));
  case 'robust'
    mld(4) = estimate_mld(z, theta, 'fixed_temperature');
    mld(3) = estimate_mld(z, rho, 'fixed_density');
    mld(2) = estimate_mld(z, rho, theta, sa, 'variable_density');
    mld(1) = estimate_mld(z, rho, 'fixed_density_gradient');
    mld = median(mld);
    [mld_i mld_i] = min(abs(z-mld));
  otherwise
    error('Unknown method');
end;

end
