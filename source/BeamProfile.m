function [estimated_params, residuals, rmse, degree_of_freedom, s0_squared, K_xx, std_params, FP_pos, Total_phase, Total_Intensity, Distance_final] = BeamProfile(pc_data, pc_profiles, tilt, ReflectanceForeground, ReflectanceBackground, Flag_side, Flag_plot)

%This function BeamProfile computes the adjusted parameters by fitting our
%predicted data to the observed experimental data. The main output for
%checking the beam dimensions are the beam shape parameter which is
%contained in the estimated_params as the first entry. There are additional
%parameters outputted by the function which are listed below. 

% input:  - pc_data[matrix]:point cloud of the both surfaces together
%           (n x 3)- :
%           1st column: X-coordinates [m]
%           2nd Column: Y-coordinates [m]
%           3rd Column: Z-coordinates [m]
%         - pc_profile[matrix]: point cloud of the band of profiles
%           (n x 3)- :
%           1st column: X-coordinates [m]
%           2nd Column: Y-coordinates [m]
%           3rd Column: Z-coordinates [m]
%         - tilt angle [deg]
%         - Reflectances of foreground and background [%]
%         - Flag_side [specifying the location of the profiles picked]
%         - Flag_plot [1 for plotting the fitted parameters]

% output: - estimated_params - the estimated parameters (4 parameters)
%           after least-squares adjustment [column vector] -> beam shape 
%           parameter, beam centering, distance of foreground and distance
%           of background plane
%         - residuals [m]: residual of the observations
%         - rmse [m]: Root Mean Square Error
%         - degree_of_freedom [scalar]: degree of freedom
%         - s0_squared [m]: a posteriori of unit weight 
%         - K_xx [m]: covariance matrix of the estimated parameters
%         - std_params [m]: standard deviations of the estimated parameters
%         - FP_pos [m]: Footprint (FP) center position with respect to the edge
%         - Total_phase [rad]: estimated total phase 
%         - Total_Intensity [au]: estimated intensity 
%         - Distance_final [m]: adjusted measured distances


%           -> Remark: Please note that this function file refers to the
%           section 5.2 in our paper, entitled as Beam Parameter Estimation

% =========================================================================


%% Essential calculations 

%convert the cartesian coordinates to spherical of the profiles point cloud
[azimuth_angles,elevation_angles,distances] = cart2sph(pc_profiles(:,1),pc_profiles(:,2),pc_profiles(:,3));

% computing the relative position of the FP center to the slanted edge
switch Flag_side
    case 'bottom'
        flag = 1;
    case 'top'
        flag = 1;
    case 'left'
        flag = 0;
    case 'right'
        flag = 0;
end        

[FP_pos] = relativeFPcenter(elevation_angles, azimuth_angles, distances, tilt, flag);


% computing the foreground and background distances and
% the respective normal vectors using k-means
[distance_foreground, distance_background, normal_foreground, ~] = estimate_fg_bg_distance(pc_data);


% computing the orthogonal measurement point and the corresponding angles to it
[azimuth_foreground, elevation_foreground] = OrthogonalPoint(normal_foreground, pc_data);


% calculate the geometric distance to the foreground and to the respective background for the corresponding angles in the profiles 
%compute difference of the angles to the orthogonal pointing angles
delta_azimuth = abs(azimuth_angles - azimuth_foreground);
delta_elevation = abs(elevation_angles - elevation_foreground);


% calculate the beam parameters (beam width and beam shape parameter sigma_b) of the Z+F Imager 5016 phase-based laser scanner. 
% ------- hard coded; not to be modified --------- %

frequency_modulation = [238095238.1, 2380952.381] ; % for z+f imager 5016 [m], assumption explained in the paper
lambda = 1500e-9; % optical wavelength for ZF scanner [m]
divergence = 0.3e-3; % half-angle [rad] %from spec-sheet

% ---- calling function to compute beam shape parameter and beam width
% diameter [1/e^2 beam diameter definition]
[~, ~, sigma_b] = GaussianBeamWidth(lambda, divergence, distance_foreground);


% calculate the initial mean for estimating the centering of our profiles
[mean] = estimate_initial_mean(distances, FP_pos); % for vertical profile 



%% ----------- Parameter Estimation part ----------- %%

format long

[estimated_params, residuals, rmse, degree_of_freedom, s0_squared, K_xx] = DistIQ_FP_Gauss_Markov_4par(distances, FP_pos, sigma_b, mean, ReflectanceForeground, ReflectanceBackground, frequency_modulation, distance_foreground, distance_background, 100, delta_azimuth, delta_elevation, Flag_side);

std_params = sqrt (diag(K_xx)) * sqrt(s0_squared); % gives the standard deviation of the parameters - uncertainty sigma


%% ------- Estimating the final distances using the estimated parameters -------- %%
%%%%%%% ----------- computing the adjusted parameters ----------- %%%%%%%

[geometric_distances_foreground, geometric_distances_background] = GeometricDistance(estimated_params(3), estimated_params(4), delta_azimuth, delta_elevation);


[Total_phase, Total_Intensity, Distance_final] = DistIQ_FP(FP_pos, estimated_params(1), estimated_params(2), ReflectanceForeground, ReflectanceBackground, geometric_distances_foreground, geometric_distances_background, frequency_modulation, Flag_side);


%% ------ Optional: plotting the fit ------  %%

if Flag_plot == 1
       
    figure
    scatter(FP_pos,distances,10,'filled')
    hold on 
    scatter(FP_pos,Distance_final,10, 'filled')

    grid on
    box on
    ylabel('Distance measured [m]', 'interpreter','latex')

    switch Flag_side
        case 'bottom'
            xlabel('Relative change $\Delta\xi$ of footprint center [m]', 'interpreter','latex')
        case 'top'
            xlabel('Relative change $\Delta\xi$ of footprint center [m]', 'interpreter','latex')
        case 'left'
            xlabel('Relative change $\Delta\eta$ of footprint center [m]', 'interpreter','latex')
        case 'right'
            xlabel('Relative change $\Delta\eta$ of footprint center [m]', 'interpreter','latex')
    end        

    set(gcf,'color','white')
    % set(gca,'FontName','Times','FontSize',30);

    h=gca;                    
    h.LineWidth= 0.75;          
    h.FontSize = 15;
    h.TickLabelInterpreter = 'latex';
    hl=legend('Measured points', 'Model predictions');            
    hl.LineWidth=0.75;           
    hl.Location = 'best';
    hl.Interpreter = 'latex';


else
    return
end


end

