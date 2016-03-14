function [ bbp, beta_p ] = estimate_bbp( beta, t, s, lambda, theta, X_p, delta )
%ESTIMATE_BBP The backscattering coefficient of particles bbp is commonly
%   estimated from measurement of scattering at a single angle in the
%   backward hemisphere beta.
%
%   $$\beta_p(\theta) &= \beta(\theta) - \beta_{sw}(\theta)\\$$
%   $$b_{bp} &= 2 \times \pi \times \chi(\theta) \times \beta_{p}(\theta)$$
%
%   beta_sw scattering at a single angle from sea water is estimated with:
%   Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
%   seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710
%   chi (X_p) is 
%
%Syntax:  [ bbp, beta_p ] = estimate_bbp( beta, t, s, lambda, theta, X_p, delta )
%
%Inputs: 
%    Required:
%        beta NxM double corresponding to the values of the VSF at one angle in m^{-1}.sr^{-1}
%        t Nx1 double corresponding to the temperature in deg C
%        s Nx1 double corresponding to the salinity in psu
%    Optional:
%        lambda 1x1 or 1xM double corresponding to the wavelength in nm
%           default: 700
%        theta 1x1 or 1xM double corresponding to the scattering angle in deg
%           default: 140 (ECO-FLBB)
%           MCOMS and ECO-FLNTU are 150
%        X_p 1x1 or 1xM double,  chi(theta) conversion coefficient at theta
%           default: interpolated from Sullivan et al. (2013)
%        delta double corresponding to the delta for Zhang et al. (2009)
%           default: 0.039 
%
%Outputs:
%   bbp NxM double corresponding to the particulate backscattering in m^{-1}
%   beta_p NxM double corresponding to particulate beta at the angle theta in m^{-1}.sr^{-1}
%
%Examples:
% [bbp_BB9, beta_p_BB9] = estimate_bbp(bb9, t, s, wl_bb9, 117);
% [bbp_VSF, beta_p_VSF] = estimate_bbp(eco_vsf, t, s, 650, scat_angle_vsf);
%
%m-files required: betasw_ZHH2009
%
% Tested with: Matlab R2015a
%              Matlab R2015b
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 11th June 2015
% Last update: 26th July 2015   

% Check Nargin
if nargin > 7
   error('Too many input arguments.')
elseif nargin < 3
   error('Not enough input arguments.')
end;
% Set default param
if ~exist('lambda','var');
  lambda = 700;
end;
if ~exist('theta','var');
  theta = 140; % 124
end;
if ~exist('delta','var');
  delta = 0.039;
end;
if ~exist('X_p','var');
  % Interpolate X_p with values from Sullivan et al. 2013
  theta_ref = 90:10:170;
  X_p_ref = [0.684 0.858 1.000 1.097 1.153 1.167 1.156 1.131 1.093];
  %sigma_ref = [0.034 0.032 0.026 0.032 0.044 0.049 0.054 0.054 0.057];
  X_p = interp1(theta_ref, X_p_ref, theta, 'spline');
end;
% Check size of input
if size(beta) ~= size(t) | size(beta) ~= size(s)
  error('beta, t and s should have same size');
end;


% Check size of lambda, theta and X_p
if size(lambda,2) == 1;
  lambda = lambda * ones(size(beta,2),1);
end;
if size(theta,2) == 1;
  theta = theta * ones(size(beta,2),1);
end;
if size(X_p,2) == 1;
  X_p = X_p * ones(size(beta,2),1);
end;

% compute beta_sw
beta_p = zeros(size(beta));
bbp = zeros(size(beta));
for i=1:size(beta,1);
  for j=1:size(beta,2);
    beta_sw = betasw_ZHH2009(lambda(j), t(i), theta(j), s(i), delta);
    beta_p(i,j) = beta(i, j) - beta_sw;
    bbp(i,j) = 2 * pi * X_p(j) * beta_p(i,j);
  end;
end;

end

