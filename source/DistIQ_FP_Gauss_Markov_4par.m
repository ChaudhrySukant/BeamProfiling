function [estimated_params, residuals, rmse, degree_of_freedom, s0_squared, K_xx] = DistIQ_FP_Gauss_Markov_4par(Distances_input, position_profile,sigma_initial, mean_initial,reflectance_foreground,reflectance_background,frequency_modulation, distance_foreground, distance_background, max_iterations, delta_azimuth, delta_elevation, Flag_side)

%This function DistIQ_FP_Gauss_Markov_4par estimates the four unknown 
%parameters (beam shape parameter, beam centering, distance to the
%foreground, and to the background surface)by performing least-squares 
%estimation (LSE)of the predictions and the observed data


% input:  - Distances_input [m]: 3D distances
%         - position_profile [m]: relative FP center positon to the edge
%         - sigma_initial [m]: initial value of beam shape parameter
%         - mean_initial [m]: intial value for centering of points in the
%           profile
%         - reflectance_foreground [%]: foreground surface reflectance
%         - reflectance_background [%]: background surface reflectance
%         - frequency_modulation [m]: modulation frequencies
%         - distance_foreground [m]: distance to foreground
%         - distance_background [m]: distance to background
%         - max_iterations [count]: max number of iterations allowed
%         - delta_azimuth [rad]: azimuth angles difference
%         - delta_elevation [rad]: elevation angles difference
%         - Flag_side: location of the profiles
%           ['left'/'right'/'top'/'bottom']


% output: - esimated_params [m]: estimated parameters [beam shape 
%           parameter, beam centering, distance of foreground and distance
%           of background plane]
%         - residuals [m]: residuals of the observations
%         - rmse [m]: Root Mean Square Error
%         - degree_of_freedom [scalar]: degree of freedom
%         - s0_squared [m]: a posteriori of unit weight
%         - k_xx [m]: covariance matrix of the estimated parameters

%           -> Remark: Please note that this function correspond to the
%           Section 5.2 Beam Parameter Estimation in our paper [Table 1 and
%           Table 3 results]

% =========================================================================


                                                                                                                                                                                        
%% observations 
% ------------------------------------------------------------------------

% vector of observations (measurements)
l = reshape(Distances_input,[length(Distances_input),1]);

% number of observations (measured points)
num_obs = length(l);


%% functional model
% ------------------------------------------------------------------------

e=10^-10;

% initialization of the vector of unknowns 
x_node = [sigma_initial; mean_initial; distance_foreground; distance_background];


    l = [l; x_node];
    l_node = zeros(num_obs+length(x_node),1);
    A = zeros(num_obs + length(x_node), length(x_node));


fprintf('start adjustment...\n')

for iteration = 1:max_iterations
    
    fprintf('iteration %d\\%d\n', iteration, max_iterations)
    
    sigma = x_node(1);
    mean = x_node(2);
    distance_foreground = x_node(3);
    distance_background = x_node(4);


   [geometric_distances_foreground, geometric_distances_background] = GeometricDistance(distance_foreground, distance_background, delta_azimuth, delta_elevation);
    
   
   distance_foreground_upper3 = distance_foreground + e;
   [geometric_distances_foreground_upper3, geometric_distances_background_upper3] = GeometricDistance(distance_foreground_upper3, distance_background, delta_azimuth, delta_elevation);


   distance_foreground_lower3 = distance_foreground - e;
   [geometric_distances_foreground_lower3, geometric_distances_background_lower3] = GeometricDistance(distance_foreground_lower3, distance_background, delta_azimuth, delta_elevation);

   distance_background_upper4 = distance_background + e;
   [geometric_distances_foreground_upper4, geometric_distances_background_upper4] = GeometricDistance(distance_foreground, distance_background_upper4, delta_azimuth, delta_elevation);


   distance_background_lower4 = distance_background - e;
   [geometric_distances_foreground_lower4, geometric_distances_background_lower4] = GeometricDistance(distance_foreground, distance_background_lower4, delta_azimuth, delta_elevation);

        
   
    for j=1:length(position_profile)
        % compute approximated observation values: f(x_node)
        
 
        [~, ~, Distance_final] = DistIQ_FP(position_profile(j), sigma, mean, reflectance_foreground, reflectance_background, geometric_distances_foreground(j), geometric_distances_background(j), frequency_modulation, Flag_side);

        l_node(j)=Distance_final;

        % compute partial derivatives w.r.t. sigma (design matrix A)
        
        [~, ~, Distance_final_upper] = DistIQ_FP(position_profile(j),(sigma+e),mean,reflectance_foreground,reflectance_background, geometric_distances_foreground(j), geometric_distances_background(j), frequency_modulation, Flag_side);
        [~, ~, Distance_final_lower] = DistIQ_FP(position_profile(j),(sigma-e),mean,reflectance_foreground,reflectance_background, geometric_distances_foreground(j), geometric_distances_background(j), frequency_modulation, Flag_side);
        A(j,1) = (Distance_final_upper-Distance_final_lower)/(2*e);
        
        % compute partial derivatives w.r.t. mean (design matrix A)
        
        [~, ~, Distance_final_upper2] = DistIQ_FP(position_profile(j), sigma,(mean+e),reflectance_foreground,reflectance_background, geometric_distances_foreground(j), geometric_distances_background(j), frequency_modulation, Flag_side);
        [~, ~, Distance_final_lower2] = DistIQ_FP(position_profile(j), sigma,(mean-e),reflectance_foreground,reflectance_background, geometric_distances_foreground(j), geometric_distances_background(j), frequency_modulation, Flag_side);
        A(j,2) = (Distance_final_upper2-Distance_final_lower2)/(2*e);
        
        
        
        % compute partial derivatives w.r.t distance foreground:

        [~, ~, Distance_final_upper3] = DistIQ_FP(position_profile(j), sigma, mean, reflectance_foreground, reflectance_background, geometric_distances_foreground_upper3(j), geometric_distances_background_upper3(j), frequency_modulation, Flag_side);
        [~, ~, Distance_final_lower3] = DistIQ_FP(position_profile(j), sigma, mean, reflectance_foreground, reflectance_background, geometric_distances_foreground_lower3(j), geometric_distances_background_lower3(j), frequency_modulation, Flag_side);
        A(j,3) = (Distance_final_upper3-Distance_final_lower3)/(2*e);
        

        
        % compute partial derivatives w.r.t distance background:

        [~, ~, Distance_final_upper4] = DistIQ_FP(position_profile(j),sigma, mean, reflectance_foreground, reflectance_background, geometric_distances_foreground_upper4(j), geometric_distances_background_upper4(j), frequency_modulation, Flag_side);
        [~, ~, Distance_final_lower4] = DistIQ_FP(position_profile(j),sigma, mean, reflectance_foreground, reflectance_background, geometric_distances_foreground_lower4(j), geometric_distances_background_lower4(j), frequency_modulation, Flag_side);
        A(j,4) = (Distance_final_upper4-Distance_final_lower4)/(2*e);
        

    end
    

        l_node(j+1:end) = x_node;
        A(j+1:end,:) = eye(length(x_node));

    
    % compute reduced observations
    dl = l - l_node;

    
    % compute adjusted increments of the unknown parameters
    dx = inv(A' * A) * A' * dl;


    % update adjusted parameters
    estimated_params = x_node + dx;

    % use the adjusted unknown parameters from the previous iteration as
    % approximated values for the current iteration
    x_node = estimated_params;    
    % break condition        
    if max(abs(dx)) < 10^-8 
        fprintf('converged...')
        break               
    end


end

%% for residual computation:

residuals = A*dx - dl;
rmse = sqrt(sum(residuals.^2)/length(residuals));

%compute degree of freedom
degree_of_freedom = (num_obs - length(estimated_params));

% Variance of the unit weight
s0_squared = dot(residuals, residuals) / (degree_of_freedom);

% covariance matrix of the unknown parameters
K_xx = inv(A'*A);
 

fprintf('\n\nadjustment completed.\n\n')

end



