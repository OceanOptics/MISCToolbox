function [ bbp ] = compute_bbp( bb, t, s, lambda, theta, delta, X )
%COMPUTE_BBP Calculate bbp from bb with
%   Xiaodong Zhang, Lianbo Hu, and Ming-Xia He (2009), Scatteirng by pure
%   seawater: Effect of salinity, Optics Express, Vol. 17, No. 7, 5698-5710
%
%Syntax:  [ bbp ] = compute_bbp( bb, lambda, theta, delta, X )
%
%Inputs: 
%    Required:
%        bb Nx1 double array corresponding to the backscattering in m^{-1}
%        t Nx1 double array corresponding to the temperature in deg C
%        s Nx1 double array corresponding to the salinity in psu
%    Optional:
%        lambda double corresponding to the wavelength in nm
%           default: 700
%        theta double corresponding to the angle in deg
%           default: 120
%        delta double corresponding to the delta for Zhang et al. (2009)
%           default: 0.039
%        X double corresponding to X(theta) 
%           default: 1.1      
%
%Outputs: bbp Nx1 double array corresponding to the particulate backscattering
%
%Example: Plot Profiles in the Artic
% bb = [2:14]*10^(-4);
% bbp = compute_bbp(bb, 700, 120);
%
%m-files required: betasw_ZHH2009
%
% Tested with: Matlab R2015a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 11th June 2015
% Last update: 11th June 2015   

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
  theta = 120;
end;
if ~exist('delta','var');
  delta = 0.039;
end;
if ~exist('X','var');
  X = 1.1;
end;
% Check size of input
if size(bb) ~= size(t) | size(bb) ~= size(s)
  error('bb, t and shoulf have same size');
end;

% compute bb_sw
bb_sw = zeros(size(bb));
for i=1:size(bb,1);
  for j=1:size(bb,2);
    bb_sw(i, j) = betasw_ZHH2009(lambda, t(i,j), theta, s(i,j), delta);
  end;
end;

% compute bbp
bbp = 2 * pi * X * (bb - bb_sw);
end

