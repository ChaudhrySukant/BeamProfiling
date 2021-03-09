function [FP_pos] = relativeFPcenter(elevation_angles, azimuth_angles, distances, tilt, Flag)

%This function relativeFPcenter computes the relative FP center distance to 
%the slanted edge position.

%Summary => This function quanitfies of the relative footprint on each of
%the targets. This is computed by taking the relative difference of the
%polar coordinates (index i and 1) of the points in the band of profiles.
%This shift resulting from the tilt angle between edge and scanner's
%vertical or horizontal axis is calculated. Combining all this allows one
%to calculate the relative change of the footprint center position with
%respect to the edge when moving from point 1 to point i.


% input:  - elevation [rad], azimuth angles [rad] and 3D distances [m] 
%           : n total number of points
%         - tilt : angle in [deg]
%         - Flag : [0 or 1] -> 0 for band of profiles along QUASI-VERTICAL 
%           and 1 for profiles along QUASI-HORIZONTAL edge.


% output: - FP_pos [column vector] : the relative change [m] of FP center 
%           to the edge when moving along the respective horizontal or
%           vertical edgethe relative 

%           -> Remark: Please note that this function correspond to the
%           Section 4.1

% =========================================================================


%% defining i/p params

tilt_angle_rad =tilt*(pi/180); % tilt angle in rad
delta_elevation = zeros(size(elevation_angles)); % defining the matrix for elevation angles differences (index i and 1)
delta_elevation_metric = ones(size(elevation_angles)); % defining the matrix for having the difference in [m]
delta_azimuth = zeros(size(azimuth_angles)); % azimuth angles difference (index i and 1)
delta_azimuth_metric = ones(size(azimuth_angles)); % the difference in [m]
FP_pos = zeros(size(azimuth_angles)); % relative FP center positions to the edge 
shift = ones(size(elevation_angles)); % shift resulting from tilt of the edge

%% computing the essential params for final relative position calculation

    if Flag == 0 % for vertical band of profiles -> computing horizontal FP distance to edge
        
        for i = 2:length(distances)
    
            delta_elevation(i)= elevation_angles(i) - elevation_angles(1,1);   
            delta_elevation_metric(i)= delta_elevation(i) * distances(i);
            shift(i)= delta_elevation_metric(i) * tan(tilt_angle_rad);
            delta_azimuth(i)= azimuth_angles(i) - azimuth_angles(1,1);
            delta_azimuth_metric(i)= delta_azimuth(i) * distances(i);
            FP_pos(i)= shift(i) + delta_azimuth_metric(i);
    
        end

    elseif Flag == 1 % for horizontal band of profiles -> computing vertical FP distance to edge 
    
        for i = 2:length(distances)
    
            delta_azimuth(i) = azimuth_angles(i)-azimuth_angles(1,1); 
            delta_azimuth_metric(i) = delta_azimuth(i) * distances(i);
            shift(i) = delta_azimuth_metric(i) * tan(tilt_angle_rad);
            delta_elevation(i) = elevation_angles(i) - elevation_angles(1,1);
            delta_elevation_metric(i) = delta_elevation(i) * distances(i);
            FP_pos(i) = shift(i) - delta_elevation_metric(i);
        
        end

    end

end