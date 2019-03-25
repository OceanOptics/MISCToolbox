function bbp = estimate_bbp_from_Rrs(Rrs443, Rrs490, Rrs55x, Rrs670,  lambda)
%ESTIMATE_BBP_FROM_RRS The backscattering coefficient of particles bbp is
%   estimated from remote sensing reflectance (Rrs).
%
%   Implementation of the Quasi-Analytical Algorithm v6 is described here:
%   http://www.ioccg.org/groups/Software_OCA/QAA_v6_2014209.pdf
%   Update from November 2014 of the IOCCG Report 5 Chapter 10: 
%   Quasi-Analytical Algorithm by ZhongPing Lee, Kendall Carder
%   and Robert Arnone
%
%   If to study long-term variations from different sensors and a task 
%   requiring much higher accuracy, the slight change in wavelength matters
%   (but should use the aw corresponding to either 555 or 547 nm though).
%   If it is simply to compare satellite results with insitu (field) data,
%   that change in wavelength does not matter as uncertainties in
%   measurements and Rrs are much greater.   
%
%Syntax:  [ bbp ] = estimate_bbp( Rrs670, Rrs55x, Rrs490, Rrs443, <lambda0> )
%
%Inputs: 
%    Required:
%        Rrs443 <NxM double> remote sensing relfectance at 443 nm (sr^{-1})
%        Rrs490 <NxM double> remote sensing relfectance at 490 nm (sr^{-1})
%        Rrs55x <NxM double> remote sensing relfectance at 55x nm (sr^{-1})
%        Rrs670 <NxM double> remote sensing relfectance at 670 nm (sr^{-1})
%    Optional:
%        lambda <1x1 double> output wavelength of bbp (nm)
%           default: 700
%
%Outputs:
%   bbp <NxM double> particulate backscattering at lambda0 (m^{-1})
%
%Example:
% Rrs443 = [0.00169, 0.00119, 0.00102, 0.00346, 0.0017, 0.0008, 0.00136,...
%    0.00213, 0.00799, 0.00247, 0.00158, 0.003, 0.00402, 0.00404, 0.0041];
% Rrs490 = [0.00329, 0.00184, 0.00151, 0.00537, 0.00255, 0.00121, 0.0022,...
%    0.00381, 0.00727, 0.00269, 0.00175, 0.00341, 0.00439, 0.00411, 0.00402];
% Rrs555 = [0.00748, 0.00425, 0.0028, 0.00907, 0.00515, 0.00269, 0.00389,...
%    0.00655, 0.00208, 0.0017, 0.0015, 0.00235, 0.00265, 0.00185, 0.00169];
% Rrs670 = [0.00346, 0.00161, 0.00179, 0.00638, 0.00268, 0.00122, 0.00182,...
%    0.00282, 0.00017, 0.00013, 0.00029, 0.0004, 0.00024, 0.00014, 0.00018];
% bbp555 = estimate_bbp_from_Rrs(Rrs443, Rrs490, Rrs555, Rrs670, 555);
%
% Tested with: Matlab R2016a
%
% Author: Nils Haentjens
% Email: nils.haentjens@maine.edu
% Created: October 11, 2016

% Check input
if nargin < 5;
  lambda = 700;
elseif nargin < 4;
  error('Expect 4 Rrs as input');
elseif nargin > 5;
  error('Too many arguments');
end;
if any(size(Rrs670) ~= size(Rrs55x)) || any(size(Rrs670) ~= size(Rrs490)) ||...
    any(size(Rrs670) ~= size(Rrs443));
  error('Rrs have different sizes');
end;
if any(size(lambda) ~= [1,1]);
  error('lambda has an unexpected size');
end;


% QAA v6
rrs443 = Rrs443 ./ (0.52 + 1.7 .* Rrs443);
rrs490 = Rrs490 ./ (0.52 + 1.7 .* Rrs490);
rrs55x = Rrs55x ./ (0.52 + 1.7 .* Rrs55x);
rrs670 = Rrs670 ./ (0.52 + 1.7 .* Rrs670);

g0 = 0.089; g1 = 0.1245; % From paper
% g0 = 0.0895; g1 = 0.1247; % From Excel
% u443 = (-g0 + sqrt(g0^2 + 4 * g1 * rrs443)) / (2 * g1);
% u490 = (-g0 + sqrt(g0^2 + 4 * g1 * rrs490)) / (2 * g1);
u55x = (-g0 + sqrt(g0^2 + 4 * g1 * rrs55x)) / (2 * g1);
u670 = (-g0 + sqrt(g0^2 + 4 * g1 * rrs670)) / (2 * g1);

% Estimate bbp(lambda0) with both method (eq. 2 & 3)
% Use 55x as ref
h0 = -1.14590292783408; h1 = -1.36582826429176; h2 = -0.469266027944581;
aw55x = 0.0596; bbw55x = 0.0009;
X = log10((rrs443 + rrs490) ./ (rrs55x + 5 * (rrs670./rrs490) .* rrs670)); % Excel use Rrs instead of rrs
a55x = aw55x + 10 .^ (h0 + h1 * X + h2 * X.^2);
bbp55x = (u55x .* a55x) ./ (1 - u55x) - bbw55x;
% Use 670 as ref
aw670 = 0.439; bbw670 = 0.00034;
a670 = aw670 + 0.39 * (Rrs670 ./ (Rrs443 + Rrs490)).^1.14;
bbp670 = (u670 .* a670) ./ (1 - u670) - bbw670;

% Keep result of one method
case_QAA_v5 = (Rrs670 < 0.0015); case_QAA_v6 = ~case_QAA_v5;
lambda0 = 555 * case_QAA_v5 + 670 * case_QAA_v6;
bbp_lambda0 = bbp55x .* case_QAA_v5 + bbp670 .* case_QAA_v6;

% Shift wavelength (eq. 5 & 6)
n = 2.0 * (1 - 1.2 * exp(-0.9 * rrs443./rrs55x)); % Excel use Rrs instead of rrs
bbp = bbp_lambda0 .* (lambda0 / lambda) .^ n;

end