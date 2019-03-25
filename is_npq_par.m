function [ p_npq ] = is_npq_par( p, par, par_threshold)
%NEED_NPQC Determine if the profile need a quenching correction (NPQ)
%   Based on PAR profile
%
% Inputs: 
%    Required:
%        p <Nx1 double> pressure or depth
%        par <Nx1 double> PAR
%    Optional:
%        par_threshold <double> threshold for PAR value (default: 80)
% Outputs:
%    p_npq <double> pressure or depth at which npq correction should be started
%         0 means that the profile is not quenched
%
% Author: Nils Haentjens
% Email: nils.haentjens@maine.edu
% Created: May25, 2018
%

% Check input
if nargin > 3; error('Too many input arguments');
elseif nargin < 2; error('Not enough input arguments'); end
if nargin < 3; par_threshold = 80; end

p_npq = p(par >= par_threshold);
if isempty(p_npq)
  p_npq = 0;
else
  p_npq = max(p_npq);
end

end







% def is_npq(_p, _par, _threshold=80):
%   # Determine if profile is affected by non photochemical quenching (NPQ)
%   #   based on PAR signal
%   # OUTPUT:
%   #   0 == means that the profile is NOT quenched
%   #   0 != means that the profile is quenched and up to which depth
%   p_npq = []
%   for p, par in zip(_p, _par):
%     if par > _threshold:
%       p_npq.append(p)
%   if p_npq:
%     return max(p_npq)
%   else:
%     return 0