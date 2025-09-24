clear all
close all
clc

%% Parameters

freq              = 1;      % Pacing freq, Hz
AF                = 0;        % 0 = no-AF; 1 = AF; NEED TO CHECK BEFORE USED
duration          = 3e3;    % [ms]
gender_value      = 0;      % 0 = Male; 1 = female (baseline)
ISO_param         = 0;        % ISO 0.1 or 0.02   % [uM] - SET LIGAND CONCENTRATION HERE 0-0.1
beat_analysis     = 0;        % Output APD90,CaT Amp, RMP...etc
Plot_Flag         = 1;
%% Run Population Cells


disp('Running AP Features')
frequencies = freq;
First = 1;
Second = 600;
[result] = A0_Obtain_intial_conditions(duration,frequencies,AF,ISO_param, gender_value, First, Second,Plot_Flag);



