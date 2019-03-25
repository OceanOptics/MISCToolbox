function [poc, poc_lower, poc_upper] = estimate_poc(bbp, lambda, method)
%ESTIMATE_POC Particulate Organic Carbon (POC) is leanearly proportional to
%   particulate backscattering bbp, various empirical relationship exist,
%   few of them are implemented in this function
%
% /!\ The calculations used are applicable only in the top layer
%     with a maximum depth defined by max(MLD, Zeu).
%
%Syntax:  [ poc ] = estimate_poc( bbp, lambda, method )
%
%Inputs:
%    Required:
%        bbp NxM double corresponding to the values of the VSF at one angle in m^{-1}
%    Optional:
%        lambda 1x1 or 1xM double corresponding to the wavelength in nm
%           default: 700
%        method string of the name of the method to use
%           default: 'SOCCOM'
%           soccom: POC = 3.12 x 10^4 x bbp(700) + 3.04
%           an emprirical relationship built for the SOCCOM floats based on
%           the relationship between the first profile of the floats and
%           in-situ measurements taken during deployement
%           (cruises: PS89, P16S and IN2015v1)
%           NAB08_down or NAB08_up: Specific to North Atlantic in Spring
%           based on empirical relationship (n=321), with data points
%           ranging between 0-600 m, recommend downast
%
%Outputs:
%   poc NxM double corresponding to poc in mg.m^{-3}
%   poc_lower NxM double with corresponding to the lower poc estimation
%   poc_upper NxM double with corresponding to the upper poc estimation
%
%Examples:
% [poc] = estimate_poc(bbp);
% [poc] = estimate_poc(bbp, 700,'soccom');
%
%References:
%     I. Cetini?? et al., Particulate organic carbon and inherent optical
%   properties during 2008 North Atlantic bloom experiment.
%   J. Geophys. Res. Ocean. 117 (2012), doi:10.1029/2011JC007771.
%     Emmanuel Boss, Marc Picheral, Thomas Leeuw, Alison Chase, Eric Karsenti,
%   Gabriel Gorsky, Lisa Taylor, Wayne Slade, Josephine Ras, and Herve Claustre.
%   The characteristics of particulate absorption, scattering and attenuation
%   coefficients in the surface ocean; Contribution of the Tara Oceans expedition.
%   Methods in Oceanography, 7:52?62, 2013.
%   ISSN 22111220. doi: 10.1016/j.mio.2013.11.002.
%   URL http://dx.doi. org/10.1016/j.mio.2013.11.002.
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
if ~exist('method','var') || isempty(method)
  method = 'SOCCOM';
end;

% Check size of input/content of input
if size(bbp,2) ~= size(lambda,2) || size(lambda,1) ~= 1
  error('bbp should be NxM and lambda should be 1xM');
end;

% Resize lambda
lambda = bsxfun(@times, ones(size(bbp)), lambda);

% switch to bbp(700)
if lambda ~= 700
  bbp_700 = bbp .* (700 ./ lambda) .^ (-0.78);
else
  bbp_700 = bbp;
end;

% Estimate poc
switch method
  case 'SOCCOM'
    % estimate poc from bbp(700)
    poc = 3.12 * 10^4 * bbp_700 + 3.04;
    poc_lower = (3.12 * 10^4 - 2.47 * 10^3) * bbp_700 + 3.04 - 6.78;
    poc_upper = (3.12 * 10^4 + 2.47 * 10^3) * bbp_700 + 3.04 + 6.78;
  case 'NAB08_up'
    % upcast
    poc = 43317 * bbp_700 - 18.4;
    poc_lower = (43317-2092) * bbp_700 - (18.4+5.8);
    poc_upper = (43317+2092) * bbp_700 - (18.4-5.8);
  case 'NAB08_down'
    % downcast
    poc = 35422 * bbp_700 - 14.4; 
    poc_lower = (35422-1754) * bbp_700 - (14.4+5.8);
    poc_upper = (35422+1754) * bbp_700 - (14.4-5.8);
  otherwise
    error('Unknown method %s', method);
end;

end

