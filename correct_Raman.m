function [Rrs412, Rrs443, Rrs490, Rrs510, Rrs555, Rrs670] = correct_Raman(Rrs412, Rrs443, Rrs490, Rrs510, Rrs555, Rrs670)
% CORRECT_RAMAN correct for Raman scattering based on Lee et al., 2013
%
%%Syntax:  [ bbp ] = estimate_bbp( Rrs670, Rrs55x, Rrs490, Rrs443, <lambda0> )
%
%Inputs: 
%    Required:
%        Rrs412 <NxM double> remote sensing relfectance at 412 nm (sr^{-1})
%        Rrs443 <NxM double> remote sensing relfectance at 443 nm (sr^{-1})
%        Rrs490 <NxM double> remote sensing relfectance at 490 nm (sr^{-1})
%        Rrs510 <NxM double> remote sensing relfectance at 510 nm (sr^{-1})
%        Rrs555 <NxM double> remote sensing relfectance at 555 nm (sr^{-1})
%        Rrs670 <NxM double> remote sensing relfectance at 670 nm (sr^{-1})
%
%Outputs:
%    Rrs412 <NxM double> remote sensing relfectance at 412 nm corrected for Raman scattering (sr^{-1})
%    Rrs443 <NxM double> remote sensing relfectance at 443 nm corrected for Raman scattering (sr^{-1})
%    Rrs490 <NxM double> remote sensing relfectance at 490 nm corrected for Raman scattering (sr^{-1})
%    Rrs510 <NxM double> remote sensing relfectance at 510 nm corrected for Raman scattering (sr^{-1})
%    Rrs555 <NxM double> remote sensing relfectance at 555 nm corrected for Raman scattering (sr^{-1})
%    Rrs670 <NxM double> remote sensing relfectance at 670 nm corrected for Raman scattering (sr^{-1})
%
%Example:
%   wl=[412,443,490,510,555,670];
%   Rrs=[0.0012	0.00169	0.00329	0.00404	0.00748	0.00346
%   0.00097	0.00119	0.00184	0.00229	0.00425	0.00161
%   0.00091	0.00102	0.00151	0.0019	0.0028	0.00179
%   0.00229	0.00346	0.00537	0.00649	0.00907	0.00638
%   0.00163	0.0017	0.00255	0.00309	0.00515	0.00268
%   0.0008	0.0008	0.00121	0.00147	0.00269	0.00122
%   0.00109	0.00136	0.0022	0.00268	0.00389	0.00182
%   0.00136	0.00213	0.00381	0.0047	0.00655	0.00282
%   0.0093	0.00799	0.00727	0.00444	0.00208	0.00017
%   0.0031	0.00247	0.00269	0.00234	0.0017	0.00013
%   0.00206	0.00158	0.00175	0.0016	0.0015	0.00029
%   0.00359	0.003	0.00341	0.00292	0.00235	0.0004
%   0.00483	0.00402	0.00439	0.00375	0.00265	0.00024
%   0.00441	0.00404	0.00411	0.00314	0.00185	0.00014
%   0.0045	0.0041	0.00402	0.00295	0.00169	0.00018];
%   [Rrs412, Rrs443, Rrs490, Rrs510, Rrs555, Rrs670] =...
%   correct_Raman(Rrs(:,1), Rrs(:,2), Rrs(:,3), Rrs(:,4), Rrs(:,5), Rrs(:,6))
%   plot(wl, Rrs, 'b'); hold('on');
%   plot(wl, [Rrs412, Rrs443, Rrs490, Rrs510, Rrs555, Rrs670], 'r');
%
% Tested with: Matlab R2016a
% 
% Author: Nils Haentjens & Emmanuel Boss
% Email: nils.haentjens@maine.edu
% Created: October 11, 2016

% Wavelength MODIS Aqua
% wl=[412,443,490,510,555,670];
% wl_aqua=[412,443,488,531,547,667];

alpha=[0.003		0.004		0.011		0.015		0.017		0.018];
beta1=[0.014		0.015		0.01		0.01		0.01		0.01];
beta2=[0.022		0.023		0.051		0.07		0.08		0.081];
% for i=1:6
%     for j=1:15
%        RF(j,i)=alpha(i)*(Rrs(j,2)/Rrs(j,5))+beta1(i)*(Rrs(j,5))^beta2(i);
%     end
% end
% Rrs=Rrs./(1+RF);

RF412 = alpha(1) * (Rrs443 ./ Rrs555) + beta1(1) * Rrs555 .^ beta2(1);
RF443 = alpha(2) * (Rrs443 ./ Rrs555) + beta1(2) * Rrs555 .^ beta2(2);
RF490 = alpha(3) * (Rrs443 ./ Rrs555) + beta1(3) * Rrs555 .^ beta2(3);
RF510 = alpha(4) * (Rrs443 ./ Rrs555) + beta1(4) * Rrs555 .^ beta2(4);
RF555 = alpha(5) * (Rrs443 ./ Rrs555) + beta1(5) * Rrs555 .^ beta2(5);
RF670 = alpha(6) * (Rrs443 ./ Rrs555) + beta1(6) * Rrs555 .^ beta2(6);

Rrs412 = Rrs412 ./ (1 + RF412);
Rrs443 = Rrs443 ./ (1 + RF443);
Rrs490 = Rrs490 ./ (1 + RF490);
Rrs510 = Rrs510 ./ (1 + RF510);
Rrs555 = Rrs555 ./ (1 + RF555);
Rrs670 = Rrs670 ./ (1 + RF670);
end