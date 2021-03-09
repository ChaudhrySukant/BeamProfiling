function [ini_mean] = estimate_initial_mean(distances,x_axis)

%This function estimate_initial_mean computes the intial centering of the 
%profiles which gives the centering of the points along the edge position

% input:  - distances: 3D distances from the polar coordinates of the point
%           cloud [m]
%         - x_axis: relative footprint distance to the slanted edge


% output: - ini_mean [m]: the centering position of the points in the
%           profile


% =========================================================================


%% computing centering of the points in the profile
n_clusters = 3;
idx = kmeans(reshape(distances, length(distances), 1), n_clusters);


means = zeros(n_clusters,1);

for i=1:n_clusters
    means(i) = mean(x_axis(idx==i));

end

ini_mean=median(means);

end

