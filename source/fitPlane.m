function plane_parameter = fitPlane( points )

%This function fitplane fits a plane to 3d points using PCA

% input:  - data[matrix]: point cloud
%           (n x 3): n total number of points
%           X (:,1) : x coordinates [m]
%           Y (:,2) : y coordinates [m]
%           Z (:,3) : z coordinates [m]

% output: - x [column vector]
%           plane parameter, where x = [a b c d]'

%           -> Remark: surface normal [a b c]' of the plane is normalized to
%              length 1

% =========================================================================

%% computing a plane fit


meanPoint = mean(points);


points_red = points - repmat(meanPoint, size(points,1), 1);

covMatrix = cov(points_red);

% Eigenvalue decomposition of the covariance matrix:
[eigenvectors, ~] = eig(covMatrix);

% Surface normal of the plane: 
n = eigenvectors(:,1);          
                                
n = n ./ norm(n);                                

% Plane parameter d: Euclidean distance of the plane from the origin
d = meanPoint * n;

% Plane parameters:
plane_parameter = [n; -d];

end
