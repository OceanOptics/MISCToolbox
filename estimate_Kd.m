function [ Kd_490, Kd_PAR ] = estimate_Kd( chl )
%ESTIMATE_KD estimate Kd(490) and Kd(PAR) as a function of chl
%   use model from Morel et al. 2007
%   assume case 1 waters

% Estimate Kd(490) from [chl]
Kd_490 = 0.0166 + 0.07242 * chl.^0.68955;
% Estimate Kd(PAR) based on Kd(490)
if nargout >= 2
  Kd_PAR = 0.0864 + 0.884 .* Kd_490 - 0.00137 .* Kd_490.^-1;
end
end

