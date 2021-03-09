close all
clear
clc 

% change the current folder to the folder of the m.file (main script):
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
clear tmp

% add all relevant folders to the Matlab path:
addpath(genpath(fullfile(cd, 'source')));

%% -------- Input data set by the user --------- %%

%% Input the path of the point clouds

%Specify the path of the point cloud file
dir_pc = '.\Data\EHR_PQ_21m.asc';             % path of the entire point cloud 
dir_profiles = '.\Data\EHR_PQ_21m_left.asc';  % path of the selected point cloud profile

% enter the essential details for the evaluation of beam dimensions
tilt = 10;                          % tilt of the foreground plane [deg]
Flag_side = 'left';               % location of the profile; choose among ['left', 'right', 'top', 'bottom']
Flag_plot = 1;                      % 1 to plot the adjusted parameters together with the experimental data; 0 otherwise
ReflectanceForeground = 0.7331;     % reflectance of the foreground plane
ReflectanceBackground = 0.7331;     % reflectance of the background plane


%% -------- DO NOT MODIFY BELOW THIS LINE --------- %%

%% load the data
fprintf('load data...\n\n')

[pc_data] = retrieve_data(dir_pc); %the entire point cloud
[pc_profiles] = retrieve_data(dir_profiles); %the point cloud of the band of profiles to be used for analysis


%% estimate beam shape parameter

[estimated_params, residuals, rmse, degree_of_freedom, s0_squared, K_xx, std_params, FP_pos, Total_phase, Total_Intensity, Distance_final] = BeamProfile(pc_data, pc_profiles, tilt, ReflectanceForeground, ReflectanceBackground, Flag_side, Flag_plot);

beam_shape_sigma = estimated_params(1)      % [m]
beam_centering = estimated_params(2)        % [m]
distance_foreground = estimated_params(3)   % [m]
distance_background = estimated_params(4)   % [m]

beam_waist_radius = beam_shape_sigma*2      % [m]



