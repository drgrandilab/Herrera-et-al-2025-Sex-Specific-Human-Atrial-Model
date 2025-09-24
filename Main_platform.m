clear all; close all; clc

%% Parameters
freq                = 1;        % Pacing freq, Hz

AF                  = 1;        % 0 = no-AF; 1 = AF;
duration            = 30e3;    % [ms]
plot_currents       = 1;        % Plot all currents
Gender_values       = 1;        % 0 = Male; 1 = female (baseline)
Plot                = 0;        % Plot Basic
Beat_analysis       = 1;        % Output APD/CaT Characterisitics
ISO_param           = 0;        % Set to 0.1uM if using
%% Run Single Cell

[t, y, result] = NH_single_cell(duration,freq,AF,Gender_values,...
plot_currents, ISO_param,Plot,Beat_analysis);
