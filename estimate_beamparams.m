clear 
close all
clc


% Script: estimate_beamparams.m: Estimtate beam waist etc. of a 
%      Gaussian beam, based on diameter or radius observations at various
%      distances (done here for 1 dimension, e.g. either hz or vt)

% Remark -> This script refers to the section 5.2 of the paper and the
% results presented in Table 2.


%% Data

% (0) observations to be used (note: distances will be treated as fixed)

% Please uncomment the data set required for the estimation of beam parameters: For new data set please replace the entries with the new experimental results

% -------------------------------------------------------------------------------
% For horizontal beam values: please uncomment these values under left/right edge
% -------------------------------------------------------------------------------
sfac_apriori = 8;     % additional factor of obs std applied to the data below
                      % should be 1, manually changed after first results
                      % if they indicate that the std of the observations
                      % is wrong

%%%%First column - average distance of estimated foreground and background distances
%%%%Second column - estimated beam shape parameter sigma
%%%%Third column - the standard deviation of the estimated parameters

obs_dat = [...% SC, 11/2020, 0.8 mm/10 m, Premium quality, left/right edge
    6.02	0.925	0.002
    6.02	0.994	0.002
    21.44	3.511	0.009
    21.44	3.496	0.007
    25.02	4.099	0.010
    25.02	4.076	0.009
    40.02	6.674	0.023
    40.02	6.854	0.024
    51.83	8.630	0.043
    51.83	8.632	0.039
];

data_name = 'ZF5016, Nov 2020, Resolution: 0.8 mm / 10 m, premium quality, hz-beam values';

% -------------------------------------------------------------------------------
% For vertical beam values: please uncomment these values under top/bottom edge
% -------------------------------------------------------------------------------

% sfac_apriori = 10;    % should be 1, manually changed after first results
                        % if they indicate that the std of the observations
                        % is wrong

%%%%First column - average distance of estimated foreground and background distances
%%%%Second column - estimated beam shape parameter sigma
%%%%Third column - the standard deviation of the estimated parameters

% obs_dat = [...% SC, 11/2020, 0.8 mm/10 m, Premium quality, top/bottom edge
%     6.02	0.967	0.002
%     6.02	1.049	0.002
%     21.44	3.825	0.010
%     21.44	3.843	0.095
%     25.02	4.488	0.013
%     25.02	4.595	0.012
%     40.02	7.722	0.030
%     40.02	7.465	0.028
%     51.83	10.039	0.044
%     51.83	9.905	0.046
% ];
% data_name = 'ZF5016, Nov 2020, Resolution: 0.8 mm / 10 m, premium quality, vertical-beam values';


%% Definition of helper parameters

fac_obs_to_1e2diameter = 4; % factor which would  convert the observations
                            % into 1/e^2 beam diameter (e.g. 2 if obs is
                            % already 1/e^2 beam radius)
                                                 
% (1) prior knowledge to be used:
    
lambda_nm = 1500;   % carrier wavelength in nm
std_lambda_nm = 0.1;  % (assumed) standard deviation of wavelength in nm
                    % (close to 0 if fixed, large value if largely unknown)
                    
w0_mm = 1;          % 1/e^2 beam waist radius in mm
std_w0_mm = 1;      % std of prior value of beam waist radius in mm
                
z0_m  = 2;          % beam waist distance in m
std_z0_m = 0.001;       % std of prior value of beam waist distance in m


%% Parameter estimation

% (2) preparations
%

np = size(obs_dat,1);    % number of original observations
n  = np + 3;                % number of total observations (incl. direct
                            % obs of the 3 beam params)
                            
u = 3;                      % number of unknown parameters 

eps_lam = 10;               % increment to lambda for numerical differentiation (in nm)
eps_w0  = 0.1;              % increment to w0 ... (in mm)
eps_z0  = 1;                % increment to z0 ... (in m)


% extract observation distances in m
d = obs_dat(:,1);

% convert observations to corresponding 1/e^2 radius
obs_val = obs_dat(:,2) * fac_obs_to_1e2diameter * 0.5;
std_obs_val = obs_dat(:,3) * fac_obs_to_1e2diameter * 0.5 * sfac_apriori ;

% provide approximations of unknown parameters (can later be modified
% to handle situation without prior values)
w0_mm_0 = w0_mm;
z0_m_0 = z0_m;
lambda_nm_0 = lambda_nm;


max_delta = inf;
it_count = 0;
xihat = zeros(3,1);

while max_delta > 1e-6

    it_count = it_count + 1;
    
    lambda_nm_0 = lambda_nm_0 + xihat(1);
    w0_mm_0 = w0_mm_0 + xihat(2);
    z0_m_0 = z0_m_0 + xihat(3);
   
    
    % calculate observations, Jacobi matrix and vcm of observations

        % the actual observations
    yobs = [obs_val;lambda_nm;w0_mm;z0_m];


        % their standard deviations
    sigobs = [std_obs_val;std_lambda_nm;std_w0_mm;std_z0_m];    

        % the calculated observations
    yc = lfn_wz(lambda_nm_0, w0_mm_0, z0_m_0, d);
    yclc = [yc;lambda_nm_0;w0_mm_0;z0_m_0];

    % the reduced observations
    y = yobs - yclc;

    % Jacobi matrix
    A = zeros(n,u);

        % first column: lambda in nm
    A(1:np,1)      = (lfn_wz(lambda_nm_0+eps_lam, w0_mm_0, z0_m_0, d) - yc)/eps_lam;

        % second column: w0h in mm
    A(1:np,2)      = (lfn_wz(lambda_nm_0, w0_mm_0+eps_w0, z0_m_0, d) - yc)/eps_w0;

        % third column: z0h in m
    A(1:np,3)      = (lfn_wz(lambda_nm_0, w0_mm_0, z0_m_0+eps_z0, d) - yc)/eps_z0;

        % direct observations of the parameters
    A(np+(1:3),:) = eye(3);


    % vcm of obs, and weight matrix
    %Syy = diag(sigobs.*sigobs);
    h= 1./sigobs;
    P = diag(h.*h);


    % estimate
    Sxixi = inv(A'*P*A);
    xihat = Sxixi*A'*P*y;

    ehat = A*xihat - y;

    s02 = ehat'*P*ehat/(n-u);


    max_delta = max(abs(xihat));
        
    if it_count > 100
        warning 'Iterations terminated on maximum iteration count.';
        break;
    end
end



lambda_nm_hat = lambda_nm_0 + xihat(1);
w0_mm_hat     = w0_mm_0 + xihat(2);
z0_m_hat      = z0_m_0 + xihat(3);

sig_lwz_hat = sqrt(diag(Sxixi));


% calculate beam divergence half angle
%theta_fa_mrad = 2*(1e-3)*lambda_nm_hat/(pi*w0_mm_hat);   % full angle 1e-3 converts from nm to mm and from rad to mrad
theta_ha_mrad = (1e-3)*lambda_nm_hat/(pi*w0_mm_hat);     % half angle


% partial derivatives of half angle in mrad w.r.t. lambda, wo and z0
hhh = [0.001/(pi*w0_mm_hat), - 0.001*lambda_nm_hat/(pi*w0_mm_hat*w0_mm_hat), 0];
sig_theta_ha = sqrt( hhh*Sxixi*hhh' );

fprintf('\n');
fprintf('#iterat= %g\n', it_count);
fprintf('s0_post= %g\n\n', sqrt(s02));
fprintf('lambda = %7.2f nm   std = %g nm\n', lambda_nm_hat, sig_lwz_hat(1));
fprintf('w0 (beam waist radius)    =   %5.2f mm   std = %g mm\n', w0_mm_hat, sig_lwz_hat(2));
fprintf('z0 (beam focal point)    =    %4.2f m    std = %g m\n\n', z0_m_hat, sig_lwz_hat(3));
fprintf('th_ha (beam divergence half-angle)  =    %4.2f mrad std = %g mrad\n', theta_ha_mrad, sig_theta_ha);



% plot data and estimation results

figure
hpz = 0:ceil(max(obs_dat(:,1)));
hwz = lfn_wz(lambda_nm_hat,w0_mm_hat,z0_m_hat,hpz);
plot(hpz,hwz,'-','linewidth',2,'color',[0 0 0.8]);
hold on;
box on
grid on
plot(hpz,-hwz,'-','linewidth',2,'color',[0 0 0.8]);

plot(obs_dat(:,1),obs_val,'o','markersize',4,'color',[0 0 0.8]);
set(gca,'fontsize',14,'linewidth',1,'ticklength',[0.02 0.01]);
xlabel('Distance /m','interpreter','latex');
ylabel('$\frac{1}{e^{2}}$ beam limits','interpreter','latex');

title(fn2fnstr(data_name), 'interpreter','latex');

h=gca;                    
h.TickLabelInterpreter = 'latex';

% ---- local functions ---

% calculate 1/e^2 beam radius at distances z (nx1 vector, values in m) 
% from the scalar wavelength lam (in nm), the scalar 1/e^2 beam waist
% radius w0 in mm, and the scalar distance z0 in m of beam waist
function wz=lfn_wz(lam,w0,z0,z)
    inp_size = size(z);    % save size for later formatting the output
    z=z(:);         % vectorize the input for easier handling
    
                    % argument dimension checking is left out, size errors
                    % will cause the function to crash
    
                    % the 1e-3 accounts for conversion nm->mm (of lam) and
                    % m->mm (of z)
    wz = w0*sqrt(1+(lam*0.001*(z-z0)/(pi*w0*w0)).^2);
    
                    % finally, reshape the output
    wz = reshape(wz,inp_size);
end

















