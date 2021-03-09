function [azimuth_foreground, elevation_foreground] = OrthogonalPoint(normal_foreground, pc_data)

%This function OrthogonalPoint calculates the azimuth and elevation angle 
%corresponding to the orthogonal measurement point on the foreground 
%surface 

%Summary => Computing the dot product of all the points in the point cloud
%using the normal vector. This is used to check the most orthogonal
%measurmenet point by taking minimum of all incidence angles computed.


% input:  - point cloud [matrix]:
%           (n x 3): n total number of points
%           X (:,1) : x coordinates [m]
%           Y (:,2) : y coordinates [m]
%           Z (:,3) : z coordinates [m]
%         - normal vector of foreground

% output: - azimuth_foreground: azimuth (horizontal) angle [rad]
%         - elevation_foreground: elevation (vertical) angle [rad]

% =========================================================================




%% computation of the orthogonal point


normal_fg = normal_foreground / norm(normal_foreground);
if normal_fg(1) < 0
    normal_fg = -normal_fg;
end

C = dot(pc_data(:,1:3), repmat(reshape(normal_fg,[1,3]),size(pc_data,1),1), 2); %using dot product : tracing the orthogonal measurement
C = acos(C ./ sqrt(sum(pc_data(:,1:3).^2,2)));

[~,idx] = min(abs(C));

[fore_azi,fore_ele,rho_fore] = cart2sph(pc_data(:,1),pc_data(:,2),pc_data(:,3));

azimuth_foreground = fore_azi(idx); %azimuth angle of the orthogonal point
elevation_foreground = fore_ele(idx); %elevation angle of the orthogonal point


end

