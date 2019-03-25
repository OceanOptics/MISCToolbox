function c_phyto = estimate_cphyto(bbp, lambda, method)
%ESTIMATE_CPHYTO PEstimate phytoplankton carbon biomass (Cphyto) from
%   particulate backscattering (bbp). The two methods available are based on
%   relation with bbp at 440 nm, the parameter lambda can be used to
%   shift the relation to 700 nm.
%
% /!\ The calculations used are applicable only in the top layer
%     with a maximum depth defined by max(MLD, Zeu).
%
%References:
%     Michael J. Behrenfeld, Emmanuel Boss, David a. Siegel, and Donald M. Shea.
%   Carbon-based ocean produc- tivity and phytoplankton physiology from space.
%   Global Biogeochemical Cycles, 19(1):1?14, 2005. ISSN 08866236. doi: 10.1029/2004GB002299.
%     Jason R. Graff, Toby K. Westberry, Allen J. Milligan, Matthew B. Brown,
%   Giorgio Dall'Olmo, Virginie van Dongen-Vogels, Kristen M. Reifel, and
%   Michael J. Behrenfeld. Analytical phytoplankton carbon measurements spanning
%   diverse ecosystems. Deep-Sea Research Part I: Oceanographic Research Papers,
%   102:16?25, 2015. ISSN 09670637. doi: 10.1016/j.dsr.2015.04.006.
%   URL http://dx.doi.org/10.1016/ j.dsr.2015.04.006.
%     Emmanuel Boss, Marc Picheral, Thomas Leeuw, Alison Chase, Eric Karsenti,
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
%        bbp NxM double, particulate backscattering (m^{-1})
%    Optional:
%        lambda 1x1 or 1xM double, wavelength (nm)
%           default: 700
%        method string, method to use for the estimation
%           Graff2015: based on empirical data (default)
%           Behrenfeld2005: based on a model
%
%Outputs:
%   Cphyto NxM double phytoplankton carbon (mg.m^{-3})
%
%Examples:
% [c_phyto] = estimate_cphyto(bbp);
% [c_phyto] = estimate_cphyto(bbp, 700,'Graff2015');
%
% Tested with: Matlab R2015b, R2017a
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: February 5th 2016
% Last update: Sept 8th 2017

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
if size(bbp,2) ~= size(lambda,2) || size(lambda,1) ~= 1
  error('bbp should be NxM and lambda should be 1xM');
end;

% Resize lambda
lambda = bsxfun(@times, ones(size(bbp)), lambda);

% Estimate c_phyto
switch method
  case 'Behrenfeld2005'
    % switch to bbp(440)
    bbp_440 = bbp .* (440 ./ lambda) .^ (-0.78);
    % estimate c_phyto from bbp(700)
    c_phyto = 13000 * (bbp_440 - 0.00035);
  case 'Graff2015'
    % switch to bbp(470)
    bbp_470 = bbp .* (470 ./ lambda) .^ (-0.78);
    % estimate c_phyto from bbp(700)
    c_phyto = 12128 * bbp_470 + 0.59;
  otherwise
    error('Unknown method %s', method);
end;

end