function [geometric_distances_foreground, geometric_distances_background] = GeometricDistance(distance_foreground, distance_background, delta_elevation, delta_azimuth)

%GeometricDistance calculates the geometrical distance to the foreground
%and background. 

%Summary => The angles (in polar coordinates) of the measurement points in
%the profiles picked define the measurement beam direction. Relative 
%difference of the angles are computed by taking the difference of the 
%angles corresponding to the picked points to the angles of most orthogonal 
%point on the surfaces (step executed before calling this function). 
%Using this information, the distances to both foreground and background 
%can be calculated. These distances are in turn used to calculate the phase
%of the modulated signal to predict the measurment point (used in other
%functions).



% input:  - foreground distance [m]
%           background distance [m]
%           vertical angles difference -> delta_elevation 
%                                         [rad, (colum vector)]
%           horizontal angles difference -> delta_azimuth
%                                           [rad, (colum vector)]

% output: - distances to foreground -> geometric_distances_foreground 
%                                      [m]
%           distances to background -> geometric_distances_background 
%                                      [m]


% =========================================================================

%% the loop implemented to compute geometrical distance to both of the surfaces for a set of elevation and azimuth angle
geometric_distances_foreground = zeros(size(delta_elevation));
geometric_distances_background = zeros(size(delta_elevation));


    for i=1:length(delta_elevation)
     
        geometric_distances_foreground(i)=distance_foreground/(cos(delta_elevation(i))*cos(delta_azimuth(i)));
        geometric_distances_background(i)=distance_background/(cos(delta_elevation(i))*cos(delta_azimuth(i)));

    end

end

