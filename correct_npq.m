function [ fchl_qc, qc_delta ] = correct_npq( fchl, z, start_npqc, optimize_start_qc, varargin )
%CORRECT_NPQ apply a non photochemical quenching (NPQ) correction on fl
%
% Inputs: 
%    Required:
%        fchl Nx1 array of double containing a profile of fluorescence chlorphyll (mg.m^-3)
%        z Nx1 array of double containing the depth of the profile (m or dBar)
%        start_npqc double containing depth (m or dBar) at which the
%           starts to be affected by quenching
%           can specify euphotic depth Zeu or MLD
%    Optional:
%        optimize_start_qc (optional but need to be set before any other input
%           true
%             instead of starting the quenching correction at the MLD it's
%             going to start at the maximum chlorophyll fluoresence value
%             at shallower depth
%           false
%             use start_npqc as starting point for the quenching correction
%        bbp Nx1 array of double containing a profile of backscattering (m^-1)
%        method string specifying the method to use
%           xing
%             Xing et al., (2012) assumes that, in the mixed layer, chlorophyll
%             concentration is homogeneous and it proposes to extrapolate up to
%             surface the highest value of chlorophyll concentration
%             encountered in the mixed layer paper.
%           xing2 (default if 3 arguments)
%             Same as method xing but take a three point median at mld instead
%             of a single value.
%           sackmann (require bbp)
%             Sackmann et al., (2008) corrects the NPQ effect by extrapolating
%             up to the surface the fluorescence value learned (LR2) below the
%             mixed layer (MC) from the relation between fl and bbp
%           boss (require iPar and no mld can be empty)
%             Method is in development
%           all (default if more than 3 arguments)
%             an average of results from xing2 and sackmann
%        mld double containing the mixed layed depth (m or dBar) required
%             for method sackman (depth up to which relation between chl
%             and bbp is learned by default it's the same as start_npqc
%        ipar0 double required when using method boss (micormol.photon.m^-2.s^-1)
%
% Outputs: fchl_qc Nx1 array of double containing chlorophyll fluorescence
%            corrected for quenching
%          qc_delta Nx1 array of double containing an estimate of the
%            absolute error of the quenching correction applied
%
% Tested with Matlab R2015a, R2015b and R2016a
%
% Required function:
%   lr2
%
% Author: Nils Haentjens, Ms, University of Maine
% Email: nils.haentjens@maine.edu
% Created: 3rd February 2015
% Last update: 12 Sept, 2016
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
if nargin > 7
   error('Too many input arguments')
elseif nargin < 1
   error('Not enough input arguments')
end

% Default arguments
if nargin < 4
  optimize_start_qc = false;
  method = 'xing2';
end;

% Get method from arguments
for i=1:nargin-4;
  if ischar(varargin{i});
    if strcmp(varargin{i}, 'xing') ||...
        strcmp(varargin{i}, 'xing2') ||...
        strcmp(varargin{i}, 'sackmann') ||...
        strcmp(varargin{i}, 'boss') ||...
        strcmp(varargin{i}, 'all');
      method = varargin{i};
      break;
    else
      error('Unknown method %s', varargin{i});
    end;
  end;
end;

% Get others arguments
for i=1:nargin-4;
  if ismatrix(varargin{i}) && ~ischar(varargin{i}) && ~isscalar(varargin{i});
    bbp = varargin{i};
  elseif isscalar(varargin{i}) && (strcmp(method, 'sackmann') || strcmp(method, 'all'))
    mld = varargin{i};
  elseif isscalar(varargin{i}) && strcmp(method, 'boss')
    ipar0 = varargin{i};
  elseif ~ischar(varargin{i})
    error('Unknown argument format');
  end;
end;

% Check missing data
if strcmp(method, 'all') || strcmp(method, 'sackmann');
  if ~exist('bbp', 'var');
    error('Method %s need variable bbp', method);
  end;
  if ~exist('mld', 'var');
    mld = start_npqc;
  end;
  if isempty(~isnan(bbp))
    warning('bbp contains only NaN values, switch to method xing2 only');
  end;
elseif strcmp(method, 'boss');
  if ~exist('ipar0', 'var');
    error('Method %s need variable ipar0', method);
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
if exist('bbp', 'var');
  s = size(bbp);
  if s(2) > s(1)
    bbp = bbp';
  end;
end;
% Set order according to depth to be sure interp1 works properly
% Shallow first and deeper last
[z, z_i] = sort(z);
fchl = fchl(z_i);
if exist('bbp', 'var');
  bbp = bbp(z_i);
end;
% Build reverse index 
for i=1:size(z_i, 1);
  ri(z_i(i)) = i;
end;

if ~isempty(start_npqc)
  % Find index of start_npqc
  [start_npqc_i start_npqc_i] = min(abs(z - start_npqc));
  if strcmp(method, 'all') || strcmp(method, 'sackmann');
    % Find index of MLD
    [mld_i mld_i] = min(abs(z - mld));
  end;
  % Optimize start_npqc with max fchl depth
  if optimize_start_qc;
    % if fchl is higher shallower then take this index as starting point
    while start_npqc_i > 2 && fchl(start_npqc_i) <= fchl(start_npqc_i-1) * 1.05
      start_npqc_i = start_npqc_i-1;
    end;
  end;
elseif ~strcmp(method, 'boss');
  error('Need start_npqc for all methods except boss');
end;

% Apply correction
switch method
  case 'xing'
    [fchl_qc, qc_delta] = xing(fchl, start_npqc_i);
    fchl_qc = fchl_qc(ri);
    qc_delta = qc_delta(ri);
  case 'xing2'
    [fchl_qc, qc_delta] = xing2(fchl, start_npqc_i);
    fchl_qc = fchl_qc(ri);
    qc_delta = qc_delta(ri);
  case 'sackmann'
    [fchl_qc, qc_delta] = sackmann(fchl, bbp, start_npqc_i, mld_i);
    fchl_qc = fchl_qc(ri);
    qc_delta = qc_delta(ri);
  case 'boss'
    [fchl_qc, qc_delta] = boss(fchl, z, ipar0);
    fchl_qc = fchl_qc(ri);
    qc_delta = qc_delta(ri);
  case 'all'
    fchl_qc1 = xing2(fchl, start_npqc_i);
    fchl_qc2 = sackmann(fchl, bbp, start_npqc_i, mld_i);
    [fchl_qc, qc_delta] = avg(fchl_qc1, fchl_qc2);
    fchl_qc = fchl_qc(ri);
    qc_delta = qc_delta(ri);
  otherwise
    error('Unknown method');
end;
end

function [fchl_qc, qc_delta] = xing(fchl, start_qc)
  % XING Correct for NPQ applying "Xing et al. (2012)" method
  fchl_qc(1:start_qc,1) = fchl(start_qc);
  fchl_qc(start_qc+1:length(fchl),1) = fchl(start_qc+1:length(fchl));
   % Update uncertainties
  qc_delta(1:start_qc,1) = abs(fchl_qc(1:start_qc,1) - fchl(1:start_qc,1));
  qc_delta(start_qc:length(fchl),1) = 0;
end

function [fchl_qc, qc_delta] = xing2(fchl, start_qc)
  % XING2 Correct for NPQ applying "Xing et al. (2012)" method
  %   taking a three point median at start qc instead of a single value
  if start_qc > 1 && start_qc < size(fchl,1);
    fchl_qc(1:start_qc,1) = nanmedian(fchl(start_qc-1:start_qc+1));
    qc_delta(1:start_qc,1) = nanmedian(fchl(start_qc-1:start_qc+1));
  elseif start_qc > 1;
    fchl_qc(1:start_qc,1) = nanmedian(fchl(start_qc-1:start_qc));
    qc_delta(1:start_qc,1) = nanstd(fchl(start_qc-1:start_qc));
    warning('xing2: start_qc at shallowest point');
  elseif start_qc < size(fchl,1);
    fchl_qc(1:start_qc,1) = nanmedian(fchl(start_qc:start_qc+1));
    qc_delta(1:start_qc,1) = nanstd(fchl(start_qc-1:start_qc));
    warning('xing2: start_qc at deepest point');
  else
    fchl_qc(1:start_qc,1) = fchl(start_qc:start_qc);
    qc_delta(1:start_qc,1) = abs(fchl_qc(1:start_qc,1) - fchl(1:start_qc,1));
    warning('xing2: only one value');
  end;
  fchl_qc(start_qc+1:length(fchl),1) = fchl(start_qc+1:length(fchl));
  % Update uncertainties
  qc_delta(start_qc:length(fchl),1) = 0;
end

function [fchl_qc, qc_delta] = sackmann(fchl, bbp, start_qc, mld)
  % SACKMANN Correct for NPQ applying "Sackmann et al. 2008" method
  % Learn from fchl > 0.05 && depth <= mld
  i = find(fchl > 0.05); % fchl > 0.05
  j = find(i >= mld); % from start_qc to deepest value
  if size(j,1) <= 5; % Minimum number of value to get a good profile (3 is the minimum for a robust lr2)
    b(1) = 0;
    b(2) = median(fchl(mld:end)) / median(bbp(mld:end)); 
  else
    [b, flag, flag, flag] = lr2(bbp(i(j)), fchl(i(j))); % if warning saying opposite sign or non correlated then use method just upper
    if flag > 0;
      b(1) = 0;
      b(2) = median(fchl(mld:end)) / median(bbp(mld:end));
    end;
  end;
  % Extrapolate value from start_qc up to the surface 
  fchl_qc(1:start_qc,1) = b(1) + bbp(1:start_qc) .* b(2);
  fchl_qc(start_qc:length(fchl),1) = fchl(start_qc:end);
  % Update uncertainties
  qc_delta(1:start_qc,1) = abs(fchl_qc(1:start_qc,1) - fchl(1:start_qc,1));
  qc_delta(start_qc:length(fchl),1) = 0;
end

function [fchl_qc, qc_delta] = boss(fchl, z, ipar0)
  warning('Method not yet available');
  fchl_qc = fchl;
  qc_delta = NaN(size(fchl,1),1);
end

function [fchl_qc, qc_delta] = avg(varargin)
  % AVG average all input array, they must be same size
  % Set varargin arrays vertically
  for i=1:nargin;
    n = size(varargin{i});
    if n(1) < n(2);
      varargin{i} = varargin{i}(:);
    end;
    varargin{i} = double(varargin{i});
  end;
  % Convert varargin to array
  var = cell2mat(varargin);
  % Compute the mean of inputs arrays
  fchl_qc = zeros(size(var,1),1);
  qc_delta = zeros(size(var,1),1);
  for i=1:size(var, 1);
    fchl_qc(i) = mean(var(i, :));
    qc_delta(i) = std(var(i, :));
  end;
end