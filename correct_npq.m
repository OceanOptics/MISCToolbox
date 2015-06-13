function [ fchl_qc ] = correct_npq( fchl, z, mld, varargin )
%CORRECT_NPQ apply a non photochemical quenching (NPQ) correction on fl
%
% Inputs: 
%    Required:
%        fchl Nx1 array of double containing a profile of fluorescence chlorphyll
%        z Nx1 array of double containing the depth of the profile
%        mld double containing the mixed layed depth
%    Optional:
%        bb Nx1 array of double containing a profile of backscattering
%        method string specifying the method to use
%           xing (default if 3 arguments)
%             Xing et al., (2012) assumes that, in the mixed layer, chlorophyll
%             concentration is homogeneous and it proposes to extrapolate up to
%             surface the highest value of chlorophyll concentration
%             encountered in the mixed layer paper.
%           sackmann (reauire bb)
%             Sackmann et al., (2008) corrects the NPQ effect by extrapolating
%             up to the surface the fluorescence value learned (LR2) below the
%             mixed layer (MC) from the relation between fl and bb
%           average (default if more than 3 arguments)
%             an average of the previous methods
%        optimize_start_qc
%           true
%             instead of starting the quenching correction at the MLD it's
%             going to start at the maximum chlorophyll fluoresence value
%             at shallower depth
%           false
%             use MLD as starting point for the quenching correction
%
% Outputs: fchl_qc Nx1 array of double containing chlorophyll fluorescence
%            corrected for quenching 
%
% Tested with Matlab R2015a
%
% Required:
%   lr2
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 3rd February 2015
% Last update: 29th May 2015
% 
% References:
%    Xing, X., Claustre, H., & Blain, S. (2012). for in vivo chlorophyll
%      fluorescence acquired by autonomous platforms: A case study with
%      instrumented elephant seals in the Kerguelen region (Southern Ocean.
%      Limnol. Oceanogr. ?, 483?495. http://doi.org/10.4319/lom.2012.10.483
%    Sackmann, B. S., Perry, M. J., & Eriksen, C. C. (2008). Seaglider
%      observations of variability in daytime fluorescence quenching of
%      chlorophyll-a in Northeastern Pacific coastal waters. Biogeosciences
%      Discussions, 5(4), 2839?2865. http://doi.org/10.5194/bgd-5-2839-2008

% Check input
if nargin > 6
   error('Too many input arguments')
elseif nargin < 1
   error('Not enough input arguments')
end

% Default arguments
if nargin == 3
  method = 'xing';
else
  method = 'all';
end;
optimize_start_qc = false;

% Get method from arguments
for i=1:nargin-3;
  if ischar(varargin{i});
    if strcmp(varargin{i}, 'xing') ||...
        strcmp(varargin{i}, 'sackmann') ||...
        strcmp(varargin{i}, 'all');
      method = varargin{i};
    else
      warning('Unknown method %s', varargin{i});
    end;
  end;
end;

% Get others arguments
for i=1:nargin-3;
  if ismatrix(varargin{i}) && ~ischar(varargin{i}) && ~isscalar(varargin{i});
    bb = varargin{i};
  elseif isscalar(varargin{i})
    optimize_start_qc = varargin{i};
  elseif ~ischar(varargin{i})
    error('Unknown argument format');
  end;
end;

% Check missing data
if strcmp(method, 'all') || strcmp(method, 'sackmann');
  if ~exist('bb', 'var');
    error('Method %s need variable bb', method);
  end;
end;

% Format Profiles
% Set orientation to Nx1 array if necessary
s = size(z);
if s(2) > s(1)
  z = z';
end;
s = size(fchl);
if s(2) > s(1)
  fchl = fchl';
end;
if exist('bb', 'var');
  s = size(bb);
  if s(2) > s(1)
    bb = bb';
  end;
end;
% Set order according to depth to be sure interp1 works properly
% Shallow first and deeper last
[z, z_i] = sort(z);
fchl = fchl(z_i);
if exist('bb', 'var');
  bb = bb(z_i);
end;
% Build reverse index 
for i=1:size(z_i, 1);
  ri(z_i(i)) = i;
end;

% Find index of MLD
[mld_i mld_i] = min(abs(z - mld));
% Start quenching correction at MLD
start_qc = mld_i;
% Optimize start_qc with max fchl depth
if optimize_start_qc;
  % if fchl is higher shallower then take this index as starting point
  while start_qc > 2 && fchl(start_qc) <= fchl(start_qc-1) * 1.05
    start_qc = start_qc-1;
  end;
end;

% Apply correction
switch method
  case 'xing'
    fchl_qc = xing(fchl, start_qc);
    fchl_qc = fchl_qc(ri);
  case 'sackmann'
    fchl_qc = sackmann(fchl, bb, start_qc);
    fchl_qc = fchl_qc(ri);
  case 'all'
    fchl_qc1 = xing(fchl, start_qc);
    fchl_qc2 = sackmann(fchl, bb, start_qc);
    fchl_qc = avg(fchl_qc1, fchl_qc2);
    fchl_qc = fchl_qc(ri);
  otherwise
    warning('Unknown method');
end;
end

function fchl_qc = xing(fchl, start_qc)
  % XING Correct for NPQ applying "Xing et al. (2012)" method
  fchl_qc(1:start_qc,1) = fchl(start_qc);
  fchl_qc(start_qc+1:length(fchl),1) = fchl(start_qc+1:length(fchl));
end

function fchl_qc = sackmann(fchl, bb, start_qc)
  % SACKMANN Correct for NPQ applying "Sackmann et al. 2008" method
  % Learn from fchl > 0.05 && depth <= start_qc
  i = find(fchl > 0.05); % fchl > 0.05
  j = find(i >= start_qc); % from start_qc to deepest value
  if size(j,1) <= 5; % Minimum number of value to get a good profile (3 is the minimum for a robust lr2)
    b(1) = 0;
    b(2) = median(fchl(start_qc:end)) / median(bb(start_qc:end)); 
  else
    [b, flag, flag, flag] = lr2(bb(i(j)), fchl(i(j))); % if warning saying opposite sign or non correlated then use method just upper
    if flag > 0;
      b(1) = 0;
      b(2) = median(fchl(start_qc:end)) / median(bb(start_qc:end));
    end;
  end;
  % Extrapolate value up to the surface 
  fchl_qc(1:start_qc,1) = b(1) + bb(1:start_qc) .* b(2);
  fchl_qc(start_qc:length(fchl),1) = fchl(start_qc:end); 
end

function fchl_qc = avg(varargin)
  % AVG average all input array, they must be same size
  % Set varargin arrays vertically
  for i=1:nargin; 
    n = size(varargin{i});
    if n(1) < n(2);
      varargin{i} = varargin{i}(:);
    end;
  end;
  % Convert varargin to array
  var = cell2mat(varargin);
  % Compute the mean of inputs arrays
  fchl_qc = zeros(size(var,1),1);
  for i=1:size(var, 1);
    fchl_qc(i) = mean(var(i, :));
  end;
end