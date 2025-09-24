close all
clear
clc

%% Parameters
% S_GNa   = par_SA(1)*(1+male_flag_p*male_change_GNa_p);
% S_GClCa = par_SA(2);
% S_GCaL  = par_SA(3);
% S_Gto   = par_SA(4)*(1+male_change_Gto_p*male_flag_p);
% S_GKur  = par_SA(5)*(1+male_change_GKur_p*male_flag_p);
% S_GKr   = par_SA(6)*(1+male_change_GKr_p*male_flag_p);
% S_GKs   = par_SA(7)*(1+male_change_GKs_p*male_flag_p);
% S_GK1   = par_SA(8);
% S_GNCX  = par_SA(9);
% S_GNaK  = par_SA(10)*(1+male_change_GNaK_p*male_flag_p);
% S_GK2P  = par_SA(11);
% S_GNaB  = par_SA(12);
% S_GCaB  = par_SA(13);
% S_GClB  = par_SA(14);
% S_GKAch = par_SA(15);
% S_GNaL = par_SA(16)*(1+male_change_GNaL_p*male_flag_p);
% S_GKp = par_SA(17);
% S_GSK = par_SA(18);
% S_GCaP  = par_SA(19)*(1+male_change_GCaP_p*male_flag_p);

parameter_names = {'GNa' 'GNaL' 'GNaB' 'VNaK' 'Gto' 'GKur'...
    'GK2P' 'GKr' 'GKs' 'GK1' 'GKp' 'GKach' 'GKCa' ...
    'GCaL' 'GCaB' 'VPMCA' 'VNCX' 'GClCa' 'GClB' 'VRYR' 'VSERCA' 'Ca Buffer-Cleft'...
    'Ca Buffer-Cyt' 'L-type k1o' 'L-type CaMKII'...
    'PLB Kmf' 'RyR EC50' 'CaSR Leak' ' CSQ Bmax' 'KiCa'...
    'KoCa' 'Kim' 'Kom'};
n_parameters = length(parameter_names);
baseline_parameters = ones(1,n_parameters);

%% Random variation
variations = 1000; % number of trials

sigmaG = 0.1*ones(1,n_parameters); % standard deviation for parameters

all_parameters = zeros(variations,n_parameters);
for ii = 1:n_parameters
    scaling = exp(sigmaG(ii)*randn(1,variations)) ;
    scaling_2 = sigmaG(ii)*randn(1,variations) ;
    newparams = baseline_parameters(ii)*scaling ;
    all_parameters(:,ii) = newparams ;
end


aaa = mean(scaling_2);
aa2 = mean(scaling);
% all_parameters % size(all_parameters)
% columns: N parameters
% rows: N trials

%%
save parameter_matrix_1000_all_parameters_1000_new_cells all_parameters
