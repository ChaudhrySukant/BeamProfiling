function [Total_phase, Total_Intensity, Distance_final] = DistIQ_FP(position_profile,sigma,mean,reflectance_foreground,reflectance_background, distances_foreground,distances_background, frequency_modulation, Flag_side)

%This function DistIQ_FP computes the distance to the target derived from 
%the total observed phase

% => Summary: Total phase is retrieved from the quadra- and
%             in-phase components (I and Q). These components depend on the
%             distance to the targets and on the reflected optical power 
%             within the beam FP. To compute the reflected power, the 
%             Gaussian irradiance profile is integrated at the edge of the 
%             targets (at each FP relative positions). The reflectance 
%             property of each of the target acts as a weight to compute 
%             the reflected optical power. Phase observation at each 
%             modulation wavelength/frequency is computed using the 
%             geometric distances which are scaled using the modulation 
%             wavelengths. Considering then the phasor sum of all 
%             contributions, the total observed phase is computed using 
%             IQ/- demodulation model. Similarly, the total intensity 
%             detected is computed as an outcome of the IQ/-
%             demodulation model. The measured distance is then computed
%             from the observed total phase.

% input:  - position_profile [m]: relative FP center positon to the edge
%         - sigma [m]: beam shape parameter
%         - mean [m]: centering of points in the profile
%         - reflectance_foreground [%]: foreground surface reflectance
%         - reflectance_background [%]: background surface reflectance
%         - distances_foreground [m]: geometric distances to
%           foreground
%         - distances_background [m]: geometric distances to
%           background
%         - frequency_modulation [m]: modulation frequencies
%         - Flag_side ['left' / 'top' / 'right' / 'bottom']: % a flag that
%           takes the side location of the edge along which the band of
%           profiles are picked

% output: - Total_phase [rad]: total computed phase for each measurement
%           angle as an output of IQ/- demodulation model
%         - Total_Intensity [AU]: total intensity as an output of IQ/-
%           demodulation model
%         - Distance_final [m]: estimated distance calculated from the 
%           total phase

%           -> Remark: Please note that this function correspond to the
%              Section 2 and partially to Section 3.1 in our paper. 
%              Furthermore, the ambiguity resolution algorithm is
%              implemented as explained in Section 2


% =========================================================================




%% Computation starts here for the respective powers required for IQ demodulation model:

fun = @(x) (exp(-0.5*((x-mean)/sigma).^2))./(sigma*sqrt(2*pi)); %%  Gaussian function (univariate)

q1=[];
q2=[];
power_fore=[];
power_bckgd=[];

for i=1:length(position_profile)
    
    q1(end+1,1)=integral(fun,-Inf,position_profile(i));
    q2(end+1,1)=integral(fun,position_profile(i),Inf);

    if strcmp(Flag_side, 'left') || strcmp(Flag_side, 'bottom')
       power_fore(end+1,1)=reflectance_foreground*q2(i);
       power_bckgd(end+1,1)=reflectance_background*q1(i);


    elseif strcmp(Flag_side, 'right') || strcmp(Flag_side, 'top')
        power_fore(end+1,1)=reflectance_foreground*q1(i);
        power_bckgd(end+1,1)=reflectance_background*q2(i);

        
    else
        msg = 'Error: the specified side location is not correct!';
        error(msg)
        return;
    end 

end


%% Computation of the phases:

c=3e8;

for i =1:length(distances_foreground)
    
    for w= 1:length(frequency_modulation)

    
    phase_fore(i,w)=(2*distances_foreground(i)*2*pi)/(c/frequency_modulation(w));
    phase_bckgd(i,w)=(2*distances_background(i)*2*pi)/(c/frequency_modulation(w));

    end
end



%% computation of the IQ demodulation for the total and the distances.


for i=1:length(power_fore)
    for w=1:length(frequency_modulation)
   
   NumIQ(i,w)= power_fore(i)*sin(phase_fore(i,w)) + power_bckgd(i)*sin(phase_bckgd(i,w));
   DenIQ(i,w)= power_fore(i)*cos(phase_fore(i,w)) + power_bckgd(i)*cos(phase_bckgd(i,w));
    
    
    end
end

wavelengths_modulation =c./frequency_modulation; 


    for i=1:length(position_profile)
        for w=1:length(frequency_modulation)

            IQ1 = NumIQ(i,w);
            IQ2 = DenIQ(i,w);
   
            Total_phase(i,w) = atan2(IQ1,IQ2);
            Total_Intensity(i,w) = sqrt(IQ1^2 + IQ2^2);

        end

        Total_phase= mod(Total_phase, 2*pi);
        
        %calling ambiguity resolution algorithm
        [Distance_final(i,1), ~ , ~, ~, ~] = Ambiguity_sol(Total_phase(i,:),wavelengths_modulation, max(wavelengths_modulation)/2);

    end



end




