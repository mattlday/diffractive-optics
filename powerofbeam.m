function [p]=powerofbeam(I,ps)
%% Integrates over an intensity profile and returns the power of the beam
% Inputs
% I - Intensity profile, cropped such that the background is minimised
% ps - Pixel size, the spacing between points in the intensity profile. If
% the profile is from a beam profiler, it is the quoted pixel size of the
% CCD.
%
% Outputs
%
%  p - The power, in the same units inputted. This may well be arbitrary,
%  but does not matter if the final power is being used in a relative way
%  to compare to other power (i.e. for diffraction efficiency calculations)

s=size(I);
y=linspace(0,s(1)*ps,s(1));
x=linspace(0,s(2)*ps,s(2));

p=trapz(x,trapz(y,I));


end