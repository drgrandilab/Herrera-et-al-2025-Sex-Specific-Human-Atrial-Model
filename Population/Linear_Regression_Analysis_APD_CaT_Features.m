%% Linear regression (clean drop-in)
clc; clear; close all;

% -------------------------
% CONFIG
% -------------------------
AF        = 0;   % 0 = nSR, 1 = AF/cAF
Gender    = 0;   % 0 = male, 1 = female
Biomarker = 1;   % 1 = Vmax, 3 = RMP, 7 = APD50, 16 = VPLT 
%     output = [dVm_max Vm_max -RMP AP_amp APD90 APD70 APD50 APD30 Ca_max...
%         Ca_min CaT_amp CaT_rise CaT_decay_50 CaT_decay_63 Na_min VPLT freq APD20];

% -------------------------
% Parameter names (33)
% -------------------------
parameter_names = { ...
  'GNa','GNaL','GNaB','VNaK','Gto','GKur','GK2P','GKr','GKs','GK1','GKp','GKach','GKCa', ...
  'GCaL','GCaB','VPMCA','VNCX','GClCa','GClB','VRYR','VSERCA','Ca Buffer-Cleft', ...
  'Ca Buffer-Cyt','L-type k1o','L-type CaMKII','PLB Kmf','RyR EC50','CaSR Leak', ...
  'CSQ Bmax','KiCa','KoCa','Kim','Kom'};

% -------------------------
% Load parameter matrix (33 cols)
% -------------------------
load parameter_matrix_1000_all_parameters.mat   % provides: all_parameters
X_full = all_parameters;                               % (N x 33)

% -------------------------
% Load biomarker matrices
% -------------------------
if AF == 1
    if Gender == 0
        % Male AF
        load AP_CaT_Features_Male_AF.mat
        Y_all = result_AP(:, Biomarker);
    else
        % Female cAF
        load AP_CaT_Features_Female_AF.mat
        Y_all = result(:, Biomarker);
    end
else
    if Gender == 0
        % Male nSR
        load AP_CaT_Features_Male_nSR.mat
        Y_all = result(:, Biomarker);
    else
        % Female nSR
        load AP_CaT_Features_Female_nSR.mat
        Y_all = result(:, Biomarker);
    end
end

% -------------------------
% Biomarker-specific transforms (match your original logic)
% -------------------------
rows_keep = true(size(Y_all));

switch Biomarker
    case 7   % APD_50 (from your comment) – keep > 50 ms subset (as before)
        rows_keep = Y_all > 50;
        Y_all     = Y_all(rows_keep);
        X_full    = X_full(rows_keep, :);

    case 16  % V_P_L_T: flip sign
        Y_all = -Y_all;

    otherwise
        % no change
end

% -------------------------
% Build X, Y (log + zscore, safe handling)
% -------------------------
% Guard against non-positive values before log
X_pos  = all(X_full > 0, 2);
Y_pos  = (Y_all > 0);

valid  = X_pos & Y_pos & all(isfinite(X_full),2) & isfinite(Y_all);
X      = X_full(valid, :);
Y      = Y_all(valid);

Xlog   = log(X);
Ylog   = log(Y);

% z-score each column of X and Y
Xz = zscore(Xlog);
Yz = zscore(Ylog);

% -------------------------
% Linear model (fitlm adds intercept by default)
% -------------------------
mdl   = fitlm(Xz, Yz);
coefs = mdl.Coefficients;                  % table
betas = coefs.Estimate(2:end);             % exclude intercept
pvals = coefs.pValue(2:end);               % exclude intercept

% -------------------------
% Sort and split by sign
% -------------------------
[betas_sorted, idx] = sort(betas, 'descend');
names_sorted         = parameter_names(idx);
pvals_sorted         = pvals(idx);

pos_mask = betas_sorted > 0;
neg_mask = betas_sorted < 0;

pos_betas = betas_sorted(pos_mask);
pos_names = names_sorted(pos_mask);
pos_pvals = pvals_sorted(pos_mask);

neg_betas = betas_sorted(neg_mask);
neg_names = names_sorted(neg_mask);
neg_pvals = pvals_sorted(neg_mask);

% -------------------------
% Plot settings by gender
% -------------------------
if Gender == 0
    col_pos = [0 0 0];         % black
    col_neg = [0.75 0.75 0.75];% gray
else
    col_pos = [1 0 0];         % red
    col_neg = [1 0.68 0.68];   % pink
end

% -------------------------
% Figure 4: Negative (asc) + Positive (desc) coefficients as two bar plots
% -------------------------
figure(4); clf; set(gcf,'Color','w','Position',[100 100 700 520]);

% Top: positive (descending)
subplot(2,1,1);
bar(pos_betas, 'FaceColor', col_pos); box off; grid on;
set(gca,'TickDir','out','LineWidth',1.5,'FontSize',12);
ylabel('Regression Coefficient','FontWeight','bold');
title('Positive coefficients (sorted)','FontWeight','bold');
xticks(1:numel(pos_names)); xticklabels(pos_names); xtickangle(90);
ylim([min(0, min(pos_betas)*1.1), max(pos_betas)*1.1]);

% Dim labels for non-significant p>0.05
ax = gca;
xt = ax.XAxis.TickLabels;
for k = 1:numel(pos_names)
    if pos_pvals(k) > 0.05
        xt{k} = sprintf('\\color[rgb]{0.5,0.5,0.5}%s', pos_names{k});
    end
end
ax.XAxis.TickLabels = xt;

% Bottom: negative (ascending)
subplot(2,1,2);
% sort negative ascending to show largest magnitude at right
[neg_betas_sorted, ii] = sort(neg_betas, 'ascend');
neg_names_sorted       = neg_names(ii);
neg_pvals_sorted       = neg_pvals(ii);

bar(neg_betas_sorted, 'FaceColor', col_neg); box off; grid on;
set(gca,'TickDir','out','LineWidth',1.5,'FontSize',12);
ylabel('Regression Coefficient','FontWeight','bold');
title('Negative coefficients (sorted)','FontWeight','bold');
xticks(1:numel(neg_names_sorted)); xticklabels(neg_names_sorted); xtickangle(90);
ylim([min(neg_betas_sorted)*1.1, max(0, max(neg_betas_sorted)*1.1)]);

ax = gca;
xt = ax.XAxis.TickLabels;
for k = 1:numel(neg_names_sorted)
    if neg_pvals_sorted(k) > 0.05
        xt{k} = sprintf('\\color[rgb]{0.5,0.5,0.5}%s', neg_names_sorted{k});
    end
end
ax.XAxis.TickLabels = xt;

% -------------------------
% Figure 5: Top-10 by absolute magnitude (sign-colored)
% -------------------------
abs_betas = abs(betas);
[abs_sorted, Iabs] = sort(abs_betas, 'descend');
keep = min(10, numel(abs_sorted));
Iabs = Iabs(1:keep);

top_names = parameter_names(Iabs);
top_vals  = betas(Iabs);
top_pvals = pvals(Iabs);

figure(5); clf; set(gcf,'Color','w','Position',[100 100 600 450]);
hold on;
for k = 1:keep
    if top_vals(k) >= 0, fc = col_pos; else, fc = col_neg; end
    bar(k, top_vals(k), 'FaceColor', fc);
end
box off; grid on; hold off;
set(gca,'TickDir','out','LineWidth',1.5,'FontSize',12);
ylabel('Regression Coefficient','FontWeight','bold');
title('Top-10 |coefficients|','FontWeight','bold');
xticks(1:keep); xticklabels(top_names); xtickangle(90);

% Dim labels for non-significant p>0.05
ax = gca;
xt = ax.XAxis.TickLabels;
for k = 1:keep
    if top_pvals(k) > 0.05
        xt{k} = sprintf('\\color[rgb]{0.5,0.5,0.5}%s', top_names{k});
    end
end
ax.XAxis.TickLabels = xt;

% -------------------------
% Console summary
% -------------------------
fprintf('\nModel N (after filtering/log): %d samples\n', size(Xz,1));
fprintf('Num positive coeffs: %d | Num negative coeffs: %d\n', numel(pos_betas), numel(neg_betas));

ns_idx = find(pvals > 0.05);
if ~isempty(ns_idx)
    fprintf('Non-significant (p>0.05):\n');
    disp(parameter_names(ns_idx)');
else
    fprintf('All coefficients are significant at p<=0.05.\n');
end
