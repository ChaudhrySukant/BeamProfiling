function [w_0, w_z_dia, sigma_b] = GaussianBeamWidth(lambda, divergence, distance_foreground)

%This function GaussianBeamWidth computes the beam width and beam shape 
%parameter at a distance of foreground [m] using standard Gaussian beam
%model equations

%Summary => 


% input:  - lambda: optical wavelength [m]
%         - divergence: angular divergence half-angle [rad]
%         - distance_foreground: foreground distance [m]


% output: - w_0 [m]: beam waist radius
%         - w_z_dia [m]: beam diameter
%         - sigma_b [m]: beam shape parameter

%           -> Remark: Please note that this function involves the Gauss 
%              beam equations used in the paper

% =========================================================================


%% computation of the beam width radius and beam shape parameter

%compute beam width radius
w_0 = lambda/(pi*divergence); % (1/e^2 [m])

%compute beam width radius at a distance 
w_z=w_0*(sqrt(1+((lambda*distance_foreground)/(pi*w_0.^2)).^2)); %beam width radius (1/e^2)[m]
w_z_dia=w_z*2; % beam width diamater (1/e^2) [m]

sigma_b=w_z_dia/4; % beam shape parameter for 1/e^2 definition. -- zf imager 5016


end

