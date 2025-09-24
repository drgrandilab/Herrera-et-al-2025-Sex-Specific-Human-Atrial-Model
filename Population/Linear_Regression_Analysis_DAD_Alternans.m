clc; clear; close all;

%% CONFIG
N_each = 10;  
Flag   = 0;   % 0 = DAD ; 1 = Alternans

parameter_names = [{'GNa'}, {'GNaL'}, {'GNaB'}, {'VNaK'}, {'Gto'}, {'GKur'},...
    {'GK2P'}, {'GKr'}, {'GKs'},{'GK1'}, {'GKp'}, {'GKach'}, {'GKCa'}, ...
    {'GCaL'}, {'GCaB'}, {'VPMCA'}, {'VNCX'}, {'GClCa'}, {'GClB'}, {'VRYR'}, {'VSERCA'}, {'Ca Buffer-Cleft'},...
    {'Ca Buffer-Cyt'}, {'L-type k1o'}, {'L-type CaMKII'}, {'PLB Kmf'}, {'RyR EC50'}, ...
    {'CaSR Leak'}, {'CSQ Bmax'}, {'KiCa'}, {'KoCa'}, {'Kim'}, {'Kom'}];

%% RUN GENDER-SPECIFIC REGRESSIONS
dataStruct = struct();

for Gender = [0 1]  % 0 = Male, 1 = Female
    load parameter_matrix_1000_all_parameters.mat
    rm_IKACH_parameters = all_parameters;

    if Gender == 0
        load DAD_Male_AF_Threshold.mat
        pacing_threshold = result_DAD(:,1)';
        pacing_threshold(pacing_threshold==1 | pacing_threshold==0) = NaN;
        DAD_NaN_values = find(isnan(pacing_threshold));

        load Alternan_Male_AF_Threshold.mat
        pacing_alt = result_Alternans(:,1)';
        pacing_alt(pacing_alt==1 | pacing_alt==0 | pacing_alt==0.1) = NaN;
        Alternan_NaN_values = find(isnan(pacing_alt));

        Total_NaN_values_unique = unique([DAD_NaN_values Alternan_NaN_values],'first');

        load AP_CaT_Features_Male_AF.mat
        result = result_AP; %#ok<NASGU>
    else
        load DAD_Female_AF_Threshold.mat
        pacing_threshold = result_DAD(:,1)';
        pacing_threshold(pacing_threshold==1 | pacing_threshold==0) = NaN;
        DAD_NaN_values = find(isnan(pacing_threshold));

        load Alternan_Female_AF_Threshold.mat
        pacing_alt = result_Alternans(:,1)';
        pacing_alt(pacing_alt==1 | pacing_alt==0 | pacing_alt==0.1) = NaN;
        Alternan_NaN_values = find(isnan(pacing_alt));

        Total_NaN_values_unique = unique([DAD_NaN_values Alternan_NaN_values],'first');

        load AP_CaT_Features_Female_AF.mat
        % keep existing 'result' if present
    end

    % Remove rows with NaNs across aligned arrays
    for i = 1:numel(Total_NaN_values_unique)
        idx = Total_NaN_values_unique(i);
        pacing_threshold(idx) = NaN;
        pacing_alt(idx)       = NaN;
        rm_IKACH_parameters(idx,:) = NaN;
    end

    % Compact to remove NaN rows
    rm_IKACH_parameters_new = [];
    for ii = 1:33
        col = rm_IKACH_parameters(:,ii);
        rm_IKACH_parameters_new(:,ii) = col(~isnan(col));
    end
    rm_IKACH_parameters = rm_IKACH_parameters_new;

    pacing_threshold = pacing_threshold(~isnan(pacing_threshold));
    DAD = pacing_threshold;

    if Flag == 0
        out_value = DAD;
    else
        out_value = pacing_alt(~isnan(pacing_alt));
    end

    % Regression (log-transform, z-score, linear regression)
    X_clean = log(rm_IKACH_parameters);
    Y_clean = log(out_value(:));

    % Drop constant columns
    valid_cols = std(X_clean,0,1) > 0;
    X_clean = X_clean(:, valid_cols);
    param_subset = parameter_names(valid_cols);

    XX = zscore(X_clean);
    YY = zscore(Y_clean);

    [b, bint] = regress(YY, [ones(size(XX,1),1) XX]);
    mdl = fitlm(XX, YY);
    p_values = table2array(mdl.Coefficients(2:end,4));

    coeffs   = b(2:end);
    conf_int = bint(2:end,:);

    % Significant if CI doesn't cross zero AND p<0.05
    sig_mask = (conf_int(:,1).*conf_int(:,2) > 0) & (p_values < 0.05);

    coeffs   = coeffs(sig_mask);
    conf_int = conf_int(sig_mask,:);
    names    = param_subset(sig_mask);

    % Split & sort
    pos_mask = coeffs > 0;
    neg_mask = coeffs < 0;

    [~, pos_sort_idx] = sort(coeffs(pos_mask), 'descend');       % largest positive
    [~, neg_sort_idx] = sort(abs(coeffs(neg_mask)), 'descend');  % most negative by |coef|

    pos_coeffs = coeffs(pos_mask);   pos_conf = conf_int(pos_mask,:); pos_names = names(pos_mask);
    neg_coeffs = coeffs(neg_mask);   neg_conf = conf_int(neg_mask,:); neg_names = names(neg_mask);

    pos_keep = min(N_each, numel(pos_sort_idx));
    neg_keep = min(N_each, numel(neg_sort_idx));

    pos_idx_keep = pos_sort_idx(1:pos_keep);
    neg_idx_keep = neg_sort_idx(1:neg_keep);

    % Extract (handle cases with < N_each or zero)
    posVals_sorted  = pos_coeffs(pos_idx_keep);
    posCI_sorted    = pos_conf(pos_idx_keep,:);
    posNames_sorted = pos_names(pos_idx_keep);

    negVals_sorted  = neg_coeffs(neg_idx_keep);
    negCI_sorted    = neg_conf(neg_idx_keep,:);
    negNames_sorted = neg_names(neg_idx_keep);

    % --- FORCE CONSISTENT SHAPES & TYPES ---
    % Names: ensure cell column; if string, convert to cellstr first
    if isstring(posNames_sorted); posNames_sorted = cellstr(posNames_sorted); end
    if isstring(negNames_sorted); negNames_sorted = cellstr(negNames_sorted); end
    if ~iscell(posNames_sorted); posNames_sorted = num2cell(posNames_sorted); end
    if ~iscell(negNames_sorted); negNames_sorted = num2cell(negNames_sorted); end

    posNames_sorted = posNames_sorted(:);           %  Nx1 cell
    negNames_sorted = negNames_sorted(:);           %  Mx1 cell

    % Values: ensure column doubles
    posVals_sorted = posVals_sorted(:);             %  Nx1 double
    negVals_sorted = negVals_sorted(:);             %  Mx1 double

    % CIs: ensure N×2 and M×2 doubles (allow true empty 0×2)
    if isempty(posCI_sorted); posCI_sorted = zeros(0,2); end
    if isempty(negCI_sorted); negCI_sorted = zeros(0,2); end
    % (If they came as row vectors like 1×2, that’s fine; vertcat with 0×2 works.)

    % --- CONCATENATE SAFELY ---
    dataStruct(Gender+1).topParams = [negNames_sorted; posNames_sorted];
    dataStruct(Gender+1).topVals   = [negVals_sorted;  posVals_sorted];
    dataStruct(Gender+1).topCI     = [negCI_sorted;    posCI_sorted];
end

%% =======================
% PLOT MALE (unchanged style)
% =======================
male_vals  = dataStruct(1).topVals;
male_ci    = dataStruct(1).topCI;
male_names = dataStruct(1).topParams;

neg_idx = find(male_vals < 0);
pos_idx = find(male_vals > 0);

adjutsed_values_bot = male_vals(neg_idx);
adjutsed_values_top = male_vals(pos_idx);
neg_ci = male_ci(neg_idx,:);
pos_ci = male_ci(pos_idx,:);
a_n_names = male_names(neg_idx);
a_p_names = male_names(pos_idx);

figure; set(gcf,'color','w');

bar(adjutsed_values_bot,'FaceColor',[0.75 0.75 0.75]); hold on;
bar(adjutsed_values_top,'FaceColor','k');

for i = 1:length(adjutsed_values_bot)
    line([i i], [neg_ci(i,1) neg_ci(i,2)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [neg_ci(i,1) neg_ci(i,1)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [neg_ci(i,2) neg_ci(i,2)], 'Color','k','LineWidth',1.5);
    text(i, adjutsed_values_bot(i)-0.25, a_n_names{i}, ...
        'Rotation',90,'FontSize',14,'HorizontalAlignment','center');
end

for i = 1:length(adjutsed_values_top)
    line([i i], [pos_ci(i,1) pos_ci(i,2)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [pos_ci(i,1) pos_ci(i,1)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [pos_ci(i,2) pos_ci(i,2)], 'Color','k','LineWidth',1.5);
    text(i, adjutsed_values_top(i)+0.25, a_p_names{i}, ...
        'Rotation',90,'FontSize',14,'HorizontalAlignment','center');
end

yGridLines = -1:0.25:1;
xLims = xlim;
for y = yGridLines
    h = line(xLims, [y y], 'Color', [0.5 0.5 0.5 0.3], 'LineStyle', '-', 'LineWidth', 1);
    uistack(h,'bottom');
end

set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('Regression Coefficient')
set(gca,'xtick',[])
set(gca,'XLim',[0 max(length(adjutsed_values_bot),length(adjutsed_values_top))+1])
set(gca,'fontname','arial','fontsize',17)
title('Male'); ylim([-1 1]); yticks([-1 -0.5 0 0.5 1])

%% =======================
% PLOT FEMALE (unchanged style)
% =======================
female_vals  = dataStruct(2).topVals;
female_ci    = dataStruct(2).topCI;
female_names = dataStruct(2).topParams;

neg_idx = find(female_vals < 0);
pos_idx = find(female_vals > 0);

adjutsed_values_bot = female_vals(neg_idx);
adjutsed_values_top = female_vals(pos_idx);
neg_ci = female_ci(neg_idx,:);
pos_ci = female_ci(pos_idx,:);
a_n_names = female_names(neg_idx);
a_p_names = female_names(pos_idx);

figure; set(gcf,'color','w');

bar(adjutsed_values_bot,'FaceColor',[1 0.68 0.68]); hold on;
bar(adjutsed_values_top,'FaceColor','r');

for i = 1:length(adjutsed_values_bot)
    line([i i], [neg_ci(i,1) neg_ci(i,2)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [neg_ci(i,1) neg_ci(i,1)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [neg_ci(i,2) neg_ci(i,2)], 'Color','k','LineWidth',1.5);
    text(i, adjutsed_values_bot(i)-0.25, a_n_names{i}, ...
        'Rotation',90,'FontSize',14,'HorizontalAlignment','center');
end

for i = 1:length(adjutsed_values_top)
    line([i i], [pos_ci(i,1) pos_ci(i,2)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [pos_ci(i,1) pos_ci(i,1)], 'Color','k','LineWidth',1.5);
    line([i-0.2 i+0.2], [pos_ci(i,2) pos_ci(i,2)], 'Color','k','LineWidth',1.5);
    text(i, adjutsed_values_top(i)+0.25, a_p_names{i}, ...
        'Rotation',90,'FontSize',14,'HorizontalAlignment','center');
end

yGridLines = -1:0.25:1;
xLims = xlim;
for y = yGridLines
    h = line(xLims, [y y], 'Color', [0.5 0.5 0.5 0.3], 'LineStyle', '-', 'LineWidth', 1);
    uistack(h,'bottom');
end

set(gca,'box','off','tickdir','out','fontsize',12)
ylabel('Regression Coefficient')
set(gca,'xtick',[])
set(gca,'XLim',[0 max(length(adjutsed_values_bot),length(adjutsed_values_top))+1])
set(gca,'fontname','arial','fontsize',17)
title('Female'); ylim([-1 1]); yticks([-1 -0.5 0 0.5 1])
