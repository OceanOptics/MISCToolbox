function c_phyto = estimate_cphyto(bbp, lambda, method)
%ESTIMATE_CPHYTO Phytoplankton carbon biomass Cphyto is estimated from bbp(440)
%   Two empirical relationship are available:
%     Michael J. Behrenfeld, Emmanuel Boss, David a. Siegel, and Donald M. Shea.
%   Carbon-based ocean produc- tivity and phytoplankton physiology from space. 
%   Global Biogeochemical Cycles, 19(1):1?14, 2005. ISSN 08866236. doi: 10.1029/2004GB002299.
%     Jason R. Graff, Toby K. Westberry, Allen J. Milligan, Matthew B. Brown,
%   Giorgio Dall?Olmo, Virginie van Dongen-Vogels, Kristen M. Reifel, and
%   Michael J. Behrenfeld. Analytical phytoplankton carbon measurements spanning
%   diverse ecosystems. Deep-Sea Research Part I: Oceanographic Research Papers, 
%   102:16?25, 2015. ISSN 09670637. doi: 10.1016/j.dsr.2015.04.006. 
%   URL http://dx.doi.org/10.1016/ j.dsr.2015.04.006.
% 
%   If require to change wavelenght of bbp the following method is used:
%   Emmanuel Boss, Marc Picheral, Thomas Leeuw, Alison Chase, Eric Karsenti,
%   Gabriel Gorsky, Lisa Taylor, Wayne Slade, Josephine Ras, and Herve Claustre.
%   The characteristics of particulate absorption, scattering and attenuation
%   coefficients in the surface ocean; Contribution of the Tara Oceans expedition.
%   Methods in Oceanography, 7:52?62, 2013.
%   ISSN 22111220. doi: 10.1016/j.mio.2013.11.002.
%   URL http://dx.doi. org/10.1016/j.mio.2013.11.002.
%
%Syntax:  [ c_phyto ] = estimate_cphyto( bbp, lambda, method )
%
%Inputs: 
%    Required:
%        bbp NxM double corresponding to the values of the VSF at one angle in m^{-1}
%    Optional:
%        lambda 1x1 or 1xM double corresponding to the wavelength in nm
%           default: 700
%        method string of the name of the method to use
%           default: 'Graff2015'
%           Behrenfeld2005: c_phyto = 13000 ? (bbp(440) ? 0.00035)
%           Graff2015: c_phyto = 12128 ? bbp(440) + 0.59
%
%Outputs:
%   bbp NxM double corresponding to c_phyto in mg.m^{?3}
%
%Examples:
% [poc] = estimate_cphyto(bbp);
% [poc] = estimate_cphyto(bbp, 700,'Graff2015');
%
% Tested with: Matlab R2015b
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: February 5th 2016
% Last update: February 5th 2016   

% Check Nargin
if nargin > 3
   error('Too many input arguments.')
elseif nargin < 1
   error('Not enough input arguments.')
end;
% Set default param
if ~exist('lambda','var');
  lambda = 700;
end;
if ~exist('method','var');
  method = 'Graff2015';
end;

% Check size of input/content of input
if size(beta,2) ~= size(lambda,2) || size(lambda,1) ~= 1
  error('bbp should be NxM and lambda should be 1xM');
end;

% Resize lambda
lambda = bsxfun(@times, ones(size(bbp)), lambda);

% Estimate poc
switch method
  case 'Behrenfeld2005'
    % switch to bbp(700)
    bbp_440 = bbp .* (440 / lambda) ^ (-0.78);
    % estimate poc from bbp(700)
    c_phyto = 13000 * (bbp_440 - 0.00035);
  case 'Graff2015'
    % switch to bbp(700)
    bbp_440 = bbp .* (440 / lambda) ^ (-0.78);
    % estimate poc from bbp(700)
    c_phyto = 12128 * bbp_440 + 0.59;
  otherwise
    error('Unknown method %s', method);
end;

end;