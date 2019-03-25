function [Lw_f0_678, Fsat] = estimate_flh(Rrs667, Rrs678, Rrs748, iPAR)
% ESTIMATE_FLH estimate fluorescence line height from remote sensing
% reflectance based on Behrenfeld, M. J., Westberry, T. K., Boss, 
% E. S., et al. (2009). Satellite-detected fluorescence reveals global
% physiology of ocean phytoplankton. Biogeosciences 6, 
% 779-795.http://dx.doi.org/10.5194/bgd-5-4235-2008
%
% NOT TESTED, nlfh available at level 2, missing Rrs748 at level 2

% Normalized water leaving radiance
F0_667 = 152.439102; F0_678 = 148.138199; F0_748 = 127.595451;
nLw667 = Rrs667 / F0_667; nLw678 = Rrs678 / F0_678; nLw748 = Rrs748 / F0_748;

% Chlorophyll fluorescence signal
Fsat = nLw678 - 70/81 * nLw667 - 11/81 * nLw748;

% Ed(0+, 670) reconstructed from iPAR using a representative spectral
% shape for the subtropics under typical sky conditions (Gregg and Carder, 1990)
Ed_0p_678 = 0.0033 * iPAR;

% Fluorescence line height
Lw_f0_678 = Fsat + Ed_0p_678 ./ F0_678; 
% with F0(678) = 148.097 mW.cm^?2.um^?1 slightly different from NASA values

end