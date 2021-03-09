function [distance_foreground, distance_background, normal_foreground, normal_background] = estimate_fg_bg_distance(data)

%This function estimate_fg_bg_distance gives an estimate of the foregroud 
%and background distance along with the normal vectors of the respective 
%surfaces

%Summary => k-means is used to segregate points of the foreground surface
%from the background plane. Following which, an outlier run is used to pick
%up only points from foreground and fit two seperate planes to each of the 
%surfaces involved


% input:  - data[matrix]: point cloud
%           (n x 3): n total number of points
%           X (:,1) : x coordinates [m]
%           Y (:,2) : y coordinates [m]
%           Z (:,3) : z coordinates [m]

% output: - distance foreground [m]
%           distance background [m]
%           normal vector foreground
%           normal vector background


% =========================================================================



%% Estimate foreground and background distance from data using k-means

[~,elevation,rho] = cart2sph(data(:,1),data(:,2),data(:,3));

dist_hor = cos(elevation).*rho;

idx = kmeans(dist_hor, 2);
y1 = dist_hor(idx == 1);
y2 = dist_hor(idx == 2);
data1 = data(idx == 1,1:3);
data2 = data(idx == 2,1:3);


% remove outliers:
[~,idx1] = rmoutliers(y1);
[~,idx2] = rmoutliers(y2);
data1 = data1(idx1 == 0,:);
data2 = data2(idx2 == 0,:);


plane1_params = fitPlane(data1);
plane2_params = fitPlane(data2);


if plane1_params(4)<0
    plane1_params = -1*plane1_params;
end

if plane2_params(4)<0
    plane2_params = -1*plane2_params;
end


if plane1_params(4) < plane2_params(4)
    distance_foreground = plane1_params(4);
    distance_background = plane2_params(4);
    normal_foreground = plane1_params(1:3);
    normal_background = plane2_params(1:3);
else
    distance_foreground = plane2_params(4);
    distance_background = plane1_params(4);
    normal_foreground = plane2_params(1:3);
    normal_background = plane1_params(1:3);
end


end

