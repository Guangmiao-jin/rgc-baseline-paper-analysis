%%% This script uses generalised linear mixed models to analyse changes in
%%% the number of light-responsive RGCs, classified as ON/OFF and
%%% transient/sustained subtypes, as well as their peak response amplitudes
%%% across age. Age is included as the main fixed effect, with retina-level
%%% variability modelled as a random effect.
%%%
%%% Rod and cone marker densities are analysed using linear
%%% mixed-effects models, with age, retinal region, and their interaction
%%% included as fixed effects, and retina-level variability included as a
%%% random effect.
%% LME VISUALIZATION: RGC SUBTYPE AMPLITUDE
%% 1. Define data directories and chronological age labels
ageFolders = {
 'D:\Data\MEA_recording\rd10_baselines\day23'
 'D:\Data\MEA_recording\rd10_baselines\day45'
 'D:\Data\MEA_recording\rd10_baselines\day60'
 'D:\Data\MEA_recording\rd10_baselines\day90'
 'D:\Data\MEA_recording\rd10_baselines\day120'
 'D:\Data\MEA_recording\rd10_baselines\day150'
 'D:\Data\MEA_recording\rd10_baselines\day200'
};
ageLabels = [23 45 60 90 120 150 200];

%% 2. Load data and aggregate into a table
T = build_retina_table(ageFolders, ageLabels);

% Define the 4 sub-types to iterate through automatically
targetGroups = {"OnTrans", "OnSus", "OffTrans", "OffSus"};

%% 3. For loop for LME modeling, ANOVA, and model diagnosis
for gi = 1:numel(targetGroups)
    currentGroup = targetGroups{gi};
    fprintf(' PROCESSING RGC SUBTYPE: %s \n', upper(currentGroup));
    
    % --- 3.1 Subsetting and Data Cleaning ---
    % Filter data for the current specific functional group
    Tg = T(T.Group == currentGroup, :);
    
    % Drop rows missing the crucial dependent variable (DeltaPeak)
    Tg = rmmissing(Tg, 'DataVariables', {'DeltaPeak'});
    
    if isempty(Tg)
        warning('No valid data found for group: %s. Skipping...', currentGroup);
        continue;
    end

    % --- 3.2 Data Type Cast for Regression and Categorical Tracking ---
    % Convert Age from categorical back to numeric double for continuous regression
    Tg.Age    = double(string(Tg.Age)); 
    Tg.Group  = categorical(Tg.Group);
    Tg.Retina = categorical(Tg.Retina);

    % Clean up unused categories to optimize categorical matrices
    Tg.Group  = removecats(Tg.Group);
    Tg.Retina = removecats(Tg.Retina);

    % --- 3.3 Linear Mixed-Effects Model (LME) Fitting ---
    % Formula: DeltaPeak ~ 1 + Age + (1|Retina)
    % Meaning: DeltaPeak is predicted by a global intercept and continuous Age (Fixed), 
    %          allowing each individual Retina to have its own random baseline intercept.
    lme = fitlme(Tg, 'DeltaPeak ~ 1 + Age + (1|Retina)');
    
    % Display fixed-effects hypothesis testing (Type III ANOVA table)
    disp('--- ANOVA Results (Fixed Effects) ---');
    disp(anova(lme))

    % --- 3.4 Model Diagnosis Plot 1: Observed vs. Predicted (Goodness of Fit) ---
    y_obs  = Tg.DeltaPeak;
    y_pred = predict(lme);   % Conditional prediction (includes estimated random effects)

    figDiag = figure('Name', sprintf('Diagnosis - %s', currentGroup), 'Color', 'w', 'Position', [100, 100, 900, 350]);
    
    % Subplot 1: Predicted vs Observed Residual Verification
    subplot(1, 3, 1); hold on;
    scatter(y_pred, y_obs, 40, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerFaceAlpha', 0.5)
    plot([min(y_obs) max(y_obs)], [min(y_obs) max(y_obs)], 'k--', 'LineWidth', 1.5)
    xlabel('Predicted \Deltapeak (spikes/s)')
    ylabel('Observed \Deltapeak (spikes/s)')
    title('Observed vs Predicted')
    axis square; box off;

    % --- 3.5 Model Diagnosis Plot 2: Q–Q Plot (Normality of Residuals) ---
    subplot(1, 3, 2);
    % Standardized normality check of model's conditional residuals
    qqplot(residuals(lme));
    title('Q–Q Plot of Residuals')
    axis square;

    % --- 3.6 Trend Visualization: Marginal Fixed Effects Prediction ---
    % Generate continuous X-axis grid based on available experimental days
    uniqueAges = sort(unique(Tg.Age));
    
    % Construct dummy table to isolate and predict population-level FIXED effects
    newT = table(uniqueAges, ...
                 categorical(repmat(currentGroup, numel(uniqueAges), 1)), ...
                 categorical(repmat("dummy", numel(uniqueAges), 1)), ...
                 'VariableNames', {'Age','Group','Retina'});

    % 'Conditional', false means we drop random effects (Z=0), plotting pure population trend
    y_trend = predict(lme, newT, 'Conditional', false);

    % Subplot 3: Raw Data Scatter vs Population Fixed Trend Curve
    subplot(1, 3, 3); hold on;
    scatter(Tg.Age, Tg.DeltaPeak, 40, 'filled', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerFaceAlpha', 0.4)
    plot(uniqueAges, y_trend, 'r-', 'LineWidth', 2.5)
    xlabel('Age (Postnatal Days)')
    ylabel('\Deltapeak (spikes/s)')
    title('Fixed-Effect Age Trend')
    axis square; box off;
    
    % Supertitle for the current loop's group figure
    sgtitle(sprintf('LME Diagnosis & Trend: %s', currentGroup), 'FontWeight', 'bold');
end


%% GLMM VISUALIZATION: RGC SUBTYPE VULNERABILITY & FRACTIONAL REGRESSION
%% 1. Configuration & Experimental Design Inputs
ageFolders = {
 'D:\Data\MEA_recording\rd10_baselines\day23'
 'D:\Data\MEA_recording\rd10_baselines\day45'
 'D:\Data\MEA_recording\rd10_baselines\day60'
 'D:\Data\MEA_recording\rd10_baselines\day90'
 'D:\Data\MEA_recording\rd10_baselines\day120'
 'D:\Data\MEA_recording\rd10_baselines\day150'
 'D:\Data\MEA_recording\rd10_baselines\day200'
};
ageLabels = [23 45 60 90 120 150 200];
% Choice 'all_units'        -> Denominator = Total recorded MEA units (nAllCells)
% Choice 'classified_only' -> Denominator = Verified responsive ON/OFF cells only
denominatorType = 'all_units'; 
%% 2. Data Gathering & Population Structure Alignment
% Collect full counts and dynamic baseline trials from directories
Tf = collect_fraction_all(ageFolders, ageLabels, denominatorType);

% Cast Age to continuous numeric double for regression fixed-effects
Tf.Age = str2double(string(Tf.Age));

% MANDATORY ENHANCEMENT: Create globally unique Retina IDs to prevent overlapping 
% random intercepts across distinct age groups (e.g., "day23_retina1")
Tf.Retina = categorical(strcat(string(Tf.Age), "_", string(Tf.Retina)));

%% 3. Generalized Linear Mixed Models (GLMM) Loop
subtypes = {'OffSus', 'OffTrans', 'OnTrans', 'OnSus'};
colors   = [0.85 0.33 0.10;  % OffSus: Burnt Orange
            0.93 0.69 0.13;  % OffTrans: Yellow/Gold
            0 0.45 0.74;     % OnTrans: Deep Blue
            0.47 0.67 0.19]; % OnSus: Apple Green

% Preallocate figure window for 2x2 multi-panel diagnostic plotting
figVulnerability = figure('Color','w','Position',[100, 100, 1050, 850]);

for i = 1:numel(subtypes)
    currentSubtype = subtypes{i};
    subplot(2, 2, i); hold on;
    
    % --- 3.1 Binomial GLMM Initialization & Fitting ---
    % Formula: SubtypeCount ~ 1 + Age + (1|Retina)
    % BinomialSize: Total units acting as N independent Bernoulli trials
    glmeFormula = sprintf('%s ~ 1 + Age + (1|Retina)', currentSubtype);
    
    glme = fitglme(Tf, glmeFormula, ...
        'Distribution', 'Binomial', ...
        'Link', 'logit', ...
        'BinomialSize', Tf.nAllCells);
    
    % --- 3.2 Raw Observation vs. Predicted Probabilities ---
    rawAge = double(Tf.Age);
    p_obs  = Tf.(currentSubtype) ./ Tf.nAllCells;
    p_pred = predict(glme); % Conditional prediction including random effects
    
    % --- 3.3 Diagnostic Quantifications & Metric Logging ---
    fixedBeta = fixedEffects(glme);             % [Intercept; Beta_Age]
    covMat    = glme.CoefficientCovariance;    % Covariance matrix of coefficients
    b_age     = fixedBeta(2);
    se_age    = sqrt(covMat(2,2));
    
    % Odds Ratio (OR) scaled per 10-day degeneration progression
    OR_10     = exp(10 * b_age);
    CI_10     = exp(10 * [b_age - 1.96*se_age, b_age + 1.96*se_age]);
    
    % Likelihood Ratio (LR) Test vs. Null Intercept Model
    nullFormula = sprintf('%s ~ 1 + (1|Retina)', currentSubtype);
    glme_null   = fitglme(Tf, nullFormula, ...
        'Distribution', 'Binomial', 'Link', 'logit', ...
        'BinomialSize', Tf.nAllCells);
    
    LR_stat = 2 * (glme.LogLikelihood - glme_null.LogLikelihood);
    p_LR    = 1 - chi2cdf(LR_stat, 1); 
    
    % Brier Score calculation (Weighted by sample size per retina)
    totalWeights = sum(Tf.nAllCells);
    brier_model  = sum(Tf.nAllCells .* (p_pred - p_obs).^2) / totalWeights;
    p_prevalence = sum(Tf.(currentSubtype)) / totalWeights;
    brier_null   = sum(Tf.nAllCells .* (p_prevalence - p_obs).^2) / totalWeights;
    
    % Output summary metrics to the MATLAB console
    fprintf('[%s] b_Age: %.4f | OR(10d): %.3f [%.3f-%.3f] | LR p: %.3e | Brier (M/N): %.3f/%.3f\n', ...
        currentSubtype, b_age, OR_10, CI_10(1), CI_10(2), p_LR, brier_model, brier_null);
    
    % --- 3.4 Plotting Stage I: Individual Sample Trajectories (Spaghetti Lines) ---
    [randEff, reInfo] = randomEffects(glme);
    b0 = randEff(strcmp(string(reInfo.Group), "Retina") & strcmp(string(reInfo.Name), "(Intercept)"));
    
    gridAge = linspace(min(rawAge), max(rawAge), 100)';
    X_matrix = [ones(numel(gridAge), 1), gridAge];
    
    % Compute the fixed population linear predictor
    eta_fixed = X_matrix * fixedBeta;
    
    % Overlay individual paths back-transformed from log-odds space
    for r = 1:numel(b0)
        y_retina_path = 1 ./ (1 + exp(-(eta_fixed + b0(r))));
        plot(gridAge, y_retina_path, 'Color', [0.88 0.88 0.88], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    
    % --- 3.5 Plotting Stage II: Population Fixed Trend + 95% Confidence Bounds ---
    se_eta = sqrt(sum((X_matrix * covMat) .* X_matrix, 2));
    y_fixed_trend = 1 ./ (1 + exp(-eta_fixed));
    ci_lower      = 1 ./ (1 + exp(-(eta_fixed - 1.96*se_eta)));
    ci_upper      = 1 ./ (1 + exp(-(eta_fixed + 1.96*se_eta)));
    
    fill([gridAge; flipud(gridAge)], [ci_lower; flipud(ci_upper)], colors(i,:), ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(gridAge, y_fixed_trend, 'Color', colors(i,:), 'LineWidth', 2.5);
    
    % --- 3.6 Plotting Stage III: Raw Retina Data Scatter Points ---
    scatter(rawAge, p_obs, 35, 'MarkerEdgeColor', colors(i,:), 'MarkerEdgeAlpha', 0.5);
    
    % --- 3.7 Plotting Stage IV: Empirical Group Medians & Interquartile Ranges (IQR) ---
    uniqueAges = unique(rawAge);
    for k = 1:numel(uniqueAges)
        empirical_p = p_obs(rawAge == uniqueAges(k));
        med_p = median(empirical_p, 'omitnan');
        iqr_p = prctile(empirical_p, [25 75]);
        
        line([uniqueAges(k) uniqueAges(k)], iqr_p, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.1, 'HandleVisibility', 'off');
        plot(uniqueAges(k), med_p, 'ko', 'MarkerFaceColor', [0.1 0.1 0.1], 'MarkerSize', 4.5, 'HandleVisibility', 'off');
    end
    
    % --- 3.8 Panel Layout Typography and Annotations ---
    title(currentSubtype, 'FontSize', 12, 'FontWeight', 'bold'); 
    box off; ylim([0 0.6]); % Fixed ceiling to prevent dynamic scaling distortions
    
    if i > 2, xlabel('Age (Postnatal Days)', 'FontSize', 10); end
    if mod(i,2) ~= 0, ylabel('Fraction of Total Units', 'FontSize', 10); end
    hudText = {
        sprintf('\\beta_{Age} = %.3f', b_age)
        sprintf('OR_{10d} = %.2f [%.2f, %.2f]', OR_10, CI_10(1), CI_10(2))
        sprintf('LR p = %.2e', p_LR)
    };
    text(ax.XLim(1) + 0.03*range(ax.XLim), ax.YLim(2) - 0.18*range(ax.YLim), hudText, ...
        'FontSize', 8.5, 'BackgroundColor', 'w', 'EdgeColor', [0.9 0.9 0.9], 'Margin', 3);
end

sgtitle('RGC Subtype Vulnerability Dynamics Across Degeneration Stages', ...
    'FontSize', 15, 'FontWeight', 'bold');



%% GLMM VISUALIZATION: CONES/RODS MARKER
%%%%%%% for loop processing staining data: INCLUDE WT IN MODEL
% 1. pre-defined excel sheet names
sheets = {'rho', 'opsinB', 'opsinRG'};
titles = {'Rhodopsin', 'OpsinB', 'OpsinRG'};

% === X axis setting: even separation ===
% rd10 ='0', 23, 45. 60, 90, 120, 150, 200
xLabels = {'0', '23', '45', '60', '90', '120', '150', '200'};
xTicks  = 0:length(xLabels)-1;  % 0,1,2,3,4,5,6,7 

% rd10 ages with WT as a psuedo zero → mapped x positions (0,1,2,3,4,5,6,7)
rd10Ages    = [0, 23, 45, 60, 90, 120, 150, 200];
rd10XPos    = 1:length(rd10Ages);  % position 0,1,2,3,4,5,6,7
ageToX = containers.Map(rd10Ages, rd10XPos);

% 2. For loop processing
for i = 1:length(sheets)
fprintf('Processing: %s...\n', titles{i});
% --- extract data and pre-processing ---
% set the directory
T = readtable("O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\whole mount raw data.xlsx", 'Sheet', sheets{i});
T.Region = categorical(T.Region);
T.Region = reordercats(T.Region, {'c','p'});
T.RetinaID = categorical(T.RetinaID);
T.logD = log(T.Density + 1e-6);
% --- LME modelling ---
lme = fitlme(T, 'logD ~ 1 + Age + Region + Age:Region + (1|RetinaID)');

% disp(anova(lme)); % print ANOVA results
% --- graph 1: Observed vs Predicted (diagnosis) ---

% fixed effects predicted values
beta    = fixedEffects(lme);
X       = designMatrix(lme, 'Fixed');
mu_fix  = X * beta;          % fixed-only predictions
var_fix = var(mu_fix);

% total fitted values (fixed + random)
mu_fit = fitted(lme); % predicted values
obs_fit = T.logD; % observed values

% Variance of random effects 
cp      = covarianceParameters(lme);
var_re  = cp{1};             % random intercept variance (RetinaID)
% Residual variance
var_res = lme.MSE;           % residual mean squared error

% Marginal R² (fixed effects only)
R2m = var_fix / (var_fix + var_re + var_res);

% Conditional R² (fixed + random effects)
R2c = (var_fix + var_re) / (var_fix + var_re + var_res);

fprintf('\n=== %s Model Fit ===\n', titles{i});
fprintf('Marginal  R² (fixed effects):           %.3f\n', R2m);
fprintf('Conditional R² (fixed + random):        %.3f\n', R2c);
fprintf('Residual variance:                       %.3f\n', var_res);
fprintf('Random effect variance (RetinaID):       %.3f\n', var_re);


figure('Name', [titles{i} ' Diagnostic']); hold on;
scatter(mu_fit, obs_fit, 60, 'filled', 'MarkerFaceAlpha', 0.5)
r = corr(mu_fit, obs_fit);
R2 = r^2;
plot([-5 5], [-5 5], 'k--', 'LineWidth', 1.5)
xlabel('Predicted log(density)'); ylabel('Observed log(density)');
title(['LME: ', titles{i}, ' Obs vs Pred']);
axis square; grid on;
x1 = xlim; y1 = ylim; 
text(x1(1)+0.05*range(x1), y1(2)-0.05*range(y1),...
    sprintf('R2 marginal = %.3f', R2m), 'FontSize',10,'BackgroundColor','w','VerticalAlignment','top')
saveas(gcf, ['O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\' sheets{i} ' fitness.png']);
close all;

% --- curve prediction and 95% CI calculation ---
ageGrid = linspace(min(T.Age), max(T.Age), 100)';

ageGrid_x = interp1([min(T.Age), max(T.Age)], ...
                        [1, length(rd10Ages)], ...
                        ageGrid);

V = lme.CoefficientCovariance;
% --- extract fixed effects coefficient ---
b0 = beta(1);
b1 = beta(2);
b2 = beta(3);
b3 = beta(4);
% --- log-scaled model ---
eq_c = sprintf('Central: D = exp(%.3f %+.3f·Age)', b0, b1);
eq_p = sprintf('Periphery: D = exp(%.3f %+.3f·Age)', b0 + b2, b1 + b3);
% construct matrix, make sure these are 'c' or 'p' not 'centre' or 'periphery'
Xc = [ones(size(ageGrid)), ageGrid, zeros(size(ageGrid)), zeros(size(ageGrid))];
Xp = [ones(size(ageGrid)), ageGrid, ones(size(ageGrid)), ageGrid];
mu_c = Xc * beta;
mu_p = Xp * beta;
se_c = sqrt(sum((Xc*V).*Xc, 2));
se_p = sqrt(sum((Xp*V).*Xp, 2));
% calculate Variance as a result of RetinaID which is a random effect
sigma_b = sqrt(cp{1});
se_re_c = sqrt(se_c.^2 + sigma_b^2);
se_re_p = sqrt(se_p.^2 + sigma_b^2);
% Exponential back-transform
Dhat_c = exp(mu_c);
Dhat_p = exp(mu_p);
CI_c_band = exp([mu_c - 1.96*se_c, mu_c + 1.96*se_c]);
CI_p_band = exp([mu_p - 1.96*se_p, mu_p + 1.96*se_p]);
RE_c_band = exp([mu_c - 1.96*se_re_c, mu_c + 1.96*se_re_c]);
RE_p_band = exp([mu_p - 1.96*se_re_p, mu_p + 1.96*se_re_p]);

T.xPos = arrayfun(@(a) ageToX(a), T.Age);

% --- graph 2: Decay Curve with Bands ---
figure('Name', [titles{i} ' Decay Curve']); hold on;
% plot Random Effect Band (background noise)
fill([ageGrid_x; flipud(ageGrid_x)], [RE_c_band(:,1); flipud(RE_c_band(:,2))], ...
[0.1 0.1 0.1], 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill([ageGrid_x; flipud(ageGrid_x)], [RE_p_band(:,1); flipud(RE_p_band(:,2))], ...
[0.1 0.1 0.1], 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
% plot Fixed Effect 95% CI (confidence interval)
fill([ageGrid_x; flipud(ageGrid_x)], [CI_c_band(:,1); flipud(CI_c_band(:,2))], ...
[1.0 0.2 0.2], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
fill([ageGrid_x; flipud(ageGrid_x)], [CI_p_band(:,1); flipud(CI_p_band(:,2))], ...
[0.2 0.2 1.0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

% plot original raw datapoints
% make sure it is labelled as 'c' or 'p'
idxC = (T.Region == 'c' | T.Region == 'central');
idxP = (T.Region == 'p' | T.Region == 'periphery');
scatter(T.xPos(idxC), T.Density(idxC), 30, 'r', 'o', 'MarkerEdgeAlpha', 0.4);
scatter(T.xPos(idxP), T.Density(idxP), 30, 'b', 'o', 'MarkerEdgeAlpha', 0.4);
% plot fit line
plot(ageGrid_x, Dhat_c, 'r-', 'LineWidth', 2.5);
plot(ageGrid_x, Dhat_p, 'b-', 'LineWidth', 2.5);
% x-/y-label
xlabel('Age (days)'); 
ylabel([titles{i}, ' density (counts/area)']);
title([titles{i}, ' Decay: Central vs Periphery']);
legend({'Central 95% CI','Periph 95% CI','Central data','Periph data','Central fit','Periph fit'}, ...
'Location','northeast');
yl = ylim;
xl = xlim;
x_text = xl(1) + 0.04 * range(xl);
y_text1 = yl(2) - 0.08 * range(yl);
y_text2 = yl(2) - 0.16 * range(yl);
text(x_text, y_text1, eq_c, 'Color', 'r', 'FontSize', 9, ...
'Interpreter', 'tex', 'FontWeight', 'bold', ...
'BackgroundColor', 'w', 'Margin', 2);
text(x_text, y_text2, eq_p, 'Color', 'b', 'FontSize', 9, ...
'Interpreter', 'tex', 'FontWeight', 'bold', ...
'BackgroundColor', 'w', 'Margin', 2);
box off; grid on;

xticks(xTicks);           % 0,1,2,3,4,5,6,7
xticklabels(xLabels);     % 0,23,45,60,90,120,150,200

hold off;
saveas(gcf, ['O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\' sheets{i} ' decay curve + 95% ci.png']);
exportgraphics(gcf, ['O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\' sheets{i} ' decay curve + 95% ci.pdf']);
close all;
end



%%%%%%for loop processing staining data: INCLUDE WT IN MODEL, but as mean ± sd
sheets = {'rho', 'opsinB', 'opsinRG'};
titles = {'Rhodopsin', 'OpsinB', 'OpsinRG'};

% === X axis setting: even separation ===
% WT=0, rd10 = 23, 45. 60, 90, 120, 150, 200
xLabels = {'23', '45', '60', '90', '120', '150', '200'};
xTicks  = 1:length(xLabels);  % 1,2,3,4,5,6,7

% rd10 actual ages → mapped x positions (WT=0，rd10 starts from 1)
rd10Ages    = [23, 45, 60, 90, 120, 150, 200];
rd10XPos    = 1:length(rd10Ages);  % position 1,2,3,4,5,6,7

% set actual age → x position projection
ageToX = containers.Map(rd10Ages, rd10XPos);

for i = 1:length(sheets)
    fprintf('Processing: %s...\n', titles{i});

    % =========================================================
    % extract data
    % =========================================================
    T = readtable("O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\whole mount raw data.xlsx", ...
                  'Sheet', sheets{i});
    T.Region   = categorical(T.Region);
    T.RetinaID = categorical(T.RetinaID);

    %  WT's Age is 0!
    isWT   = (T.Age == 0);  
    isRD10 = (T.Age > 0);

    T_wt   = T(isWT,  :);
    T_rd10 = T(isRD10,:);

    % LME(only rd10 data is included）
    T_rd10.logD = log(T_rd10.Density + 1e-6);

    lme = fitlme(T_rd10, ...
        'logD ~ 1 + Age + Region + Age:Region + (1|RetinaID)');

    % diagnosis（Observed vs Predicted）
    mu_fit  = fitted(lme);
    obs_fit = T_rd10.logD;

    figure('Name', [titles{i} ' Diagnostic']); hold on;
    scatter(mu_fit, obs_fit, 60, 'filled', 'MarkerFaceAlpha', 0.5)
    plot([-5 5], [-5 5], 'k--', 'LineWidth', 1.5)
    xlabel('Predicted log(density)');
    ylabel('Observed log(density)');
    title(['LME: ', titles{i}, ' Obs vs Pred']);
    axis square; grid on;
    close all;

    % fixed effects predicted values
    beta    = fixedEffects(lme);
    X       = designMatrix(lme, 'Fixed');
    mu_fix  = X * beta;          % fixed-only predictions
    var_fix = var(mu_fix);

    % Variance of random effects 
    cp      = covarianceParameters(lme);
    var_re  = cp{1};             % random intercept variance (RetinaID)
    % Residual variance
    var_res = lme.MSE;           % residual mean squared error
    % Marginal R² (fixed effects only)
    R2m = var_fix / (var_fix + var_re + var_res);
    % Conditional R² (fixed + random effects)
    R2c = (var_fix + var_re) / (var_fix + var_re + var_res);


    % curve prediction on the evenly separated x-axis
    % =========================================================
    % model is fitting with the real number from 23 to 200
    % predict first and then project to 1-7

    ageGrid_real = linspace(min(T_rd10.Age), max(T_rd10.Age), 200)';

    % Map real ages to x positions (linear interpolation)
    % min age → x=1, max age → x=7
    ageGrid_x = interp1([min(T_rd10.Age), max(T_rd10.Age)], ...
                        [1, length(rd10Ages)], ...
                        ageGrid_real);

    % Fixed effects
    beta = fixedEffects(lme);
    V    = lme.CoefficientCovariance;

    b0 = beta(1);  % intercept
    b1 = beta(2);  % age slope
    b2 = beta(3);  % region effect (periphery vs center)
    b3 = beta(4);  % age:region interaction

    % Design matrices（those predicted with real ages）
    Xc = [ones(size(ageGrid_real)), ageGrid_real, ...
          zeros(size(ageGrid_real)), zeros(size(ageGrid_real))];
    Xp = [ones(size(ageGrid_real)), ageGrid_real, ...
          ones(size(ageGrid_real)), ageGrid_real];

    mu_c = Xc * beta;
    mu_p = Xp * beta;

    se_c = sqrt(sum((Xc*V).*Xc, 2));
    se_p = sqrt(sum((Xp*V).*Xp, 2));

    % Random effect variance
    cp      = covarianceParameters(lme);
    sigma_b = sqrt(cp{1});

    se_re_c = sqrt(se_c.^2 + sigma_b^2);
    se_re_p = sqrt(se_p.^2 + sigma_b^2);

    % Back-transform
    Dhat_c = exp(mu_c);
    Dhat_p = exp(mu_p);

    CI_c_band = exp([mu_c - 1.96*se_c,    mu_c + 1.96*se_c]);
    CI_p_band = exp([mu_p - 1.96*se_p,    mu_p + 1.96*se_p]);
    RE_c_band = exp([mu_c - 1.96*se_re_c, mu_c + 1.96*se_re_c]);
    RE_p_band = exp([mu_p - 1.96*se_re_p, mu_p + 1.96*se_re_p]);

    % WT summary statistics（for showing the data point）
    idxC_wt = (T_wt.Region == 'c' | T_wt.Region == 'central');
    idxP_wt = (T_wt.Region == 'p' | T_wt.Region == 'periphery');

    wt_c_mean = mean(T_wt.Density(idxC_wt));
    wt_p_mean = mean(T_wt.Density(idxP_wt));
    wt_c_sem  = std(T_wt.Density(idxC_wt))  / sqrt(sum(idxC_wt));
    wt_p_sem  = std(T_wt.Density(idxP_wt))  / sqrt(sum(idxP_wt));

    % save reuslts in Excel sheets
    % 1. Raw data table
    output_raw = table();
    output_raw.RetinaID = T_rd10.RetinaID;
    output_raw.Age      = T_rd10.Age;
    output_raw.Region   = T_rd10.Region;
    output_raw.Density  = T_rd10.Density;

    % 2. WT data summary
    output_wt = table();
    output_wt.Group  = {'WT Central'; 'WT Periphery'};
    output_wt.Mean   = [wt_c_mean; wt_p_mean];
    output_wt.SEM    = [wt_c_sem;  wt_p_sem];

    % 3. Model predictions（fitted curve)
    keyAges = [23, 45, 60, 90, 120, 150, 200]';
    keyAges_x = arrayfun(@(a) ageToX(a), keyAges);

    Xc_key = [ones(size(keyAges)), keyAges, ...
           zeros(size(keyAges)), zeros(size(keyAges))];
    Xp_key = [ones(size(keyAges)), keyAges, ...
           ones(size(keyAges)), keyAges];
    mu_c_key = Xc_key * beta;
    mu_p_key = Xp_key * beta;
    se_c_key = sqrt(sum((Xc_key*V).*Xc_key, 2));
    se_p_key = sqrt(sum((Xp_key*V).*Xp_key, 2));
    output_pred = table();
    output_pred.Age              = keyAges;
    output_pred.Central_mean     = exp(mu_c_key);
    output_pred.Central_CI_lower = exp(mu_c_key - 1.96*se_c_key);
    output_pred.Central_CI_upper = exp(mu_c_key + 1.96*se_c_key);
    output_pred.Periph_mean      = exp(mu_p_key);
    output_pred.Periph_CI_lower  = exp(mu_p_key - 1.96*se_p_key);
    output_pred.Periph_CI_upper  = exp(mu_p_key + 1.96*se_p_key);

    % 4. model summary
    output_model = table();
    output_model.Parameter = {'b0 (intercept)'; 
                          'b1 (age slope)'; 
                          'b2 (region)'; 
                          'b3 (age:region)';
                          'R2_marginal';
                          'R2_conditional';
                          'AIC';
                          'BIC'};
    output_model.Value = [b0; b1; b2; b3; 
                      R2m; R2c;
                      lme.ModelCriterion.AIC;
                      lme.ModelCriterion.BIC];
    % 5. Write to Excel
    excelPath = ['O:\stemcell\MRC-Subretinal Transplantation Project\' ...
             'Rd10 paper\Figures\ephys supplementary\' ...
             sheets{i} '_prism_data.xlsx'];
    writetable(output_raw,   excelPath, 'Sheet', 'Raw Data');
    writetable(output_wt,    excelPath, 'Sheet', 'WT Reference');
    writetable(output_pred,  excelPath, 'Sheet', 'Model Predictions');
    writetable(output_model, excelPath, 'Sheet', 'Model Summary');
    
    fprintf('Excel exported: %s\n', excelPath);

    % project the x-positions of rd10 predicted datapoints to 1-7
    T_rd10.xPos = arrayfun(@(a) ageToX(a), T_rd10.Age);

    idxC_rd10 = (T_rd10.Region == 'c' | T_rd10.Region == 'central');
    idxP_rd10 = (T_rd10.Region == 'p' | T_rd10.Region == 'periphery');

    % Plotting：Decay Curve with WT reference
    figure('Name', [titles{i} ' Decay Curve']);
    hold on;

    % --- Random Effect bands（on evenly separated x-axis）---
    fill([ageGrid_x; flipud(ageGrid_x)], ...
         [RE_c_band(:,1); flipud(RE_c_band(:,2))], ...
         [0.1 0.1 0.1], 'FaceAlpha', 0.05, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    fill([ageGrid_x; flipud(ageGrid_x)], ...
         [RE_p_band(:,1); flipud(RE_p_band(:,2))], ...
         [0.1 0.1 0.1], 'FaceAlpha', 0.05, ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');

    % --- Fixed Effect 95% CI bands ---
    fill([ageGrid_x; flipud(ageGrid_x)], ...
         [CI_c_band(:,1); flipud(CI_c_band(:,2))], ...
         [1.0 0.2 0.2], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    fill([ageGrid_x; flipud(ageGrid_x)], ...
         [CI_p_band(:,1); flipud(CI_p_band(:,2))], ...
         [0.2 0.2 1.0], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    % --- rd10 raw data on evenly separately x-axis
    scatter(T_rd10.xPos(idxC_rd10), T_rd10.Density(idxC_rd10), ...
            30, 'r', 'o', 'MarkerEdgeAlpha', 0.4, ...
            'HandleVisibility', 'off');
    scatter(T_rd10.xPos(idxP_rd10), T_rd10.Density(idxP_rd10), ...
            30, 'b', 'o', 'MarkerEdgeAlpha', 0.4, ...
            'HandleVisibility', 'off');

    % --- Fitted lines on evenly separated x-axis---
    plot(ageGrid_x, Dhat_c, 'r-', 'LineWidth', 2.5);
    plot(ageGrid_x, Dhat_p, 'b-', 'LineWidth', 2.5);

    % --- WT data point at x=0 ---
    % Mean ± SEM errorbar
    % --- WT mean as horizontal dotted reference lines ---
    yline(wt_c_mean, 'r--', ...
      'LineWidth', 1.5, ...
      'Label', sprintf('WT Central (%.0f)', wt_c_mean), ...
      'LabelHorizontalAlignment', 'right', ...
      'LabelVerticalAlignment', 'top', ...
      'FontSize', 8, ...
      'HandleVisibility', 'off');
    
    yline(wt_p_mean, 'b--', ...
      'LineWidth', 1.5, ...
      'Label', sprintf('WT Periph (%.0f)', wt_p_mean), ...
      'LabelHorizontalAlignment', 'right', ...
      'LabelVerticalAlignment', 'bottom', ...
      'FontSize', 8, ...
      'HandleVisibility', 'off');

    % --- x-/y-labels ---
    xlabel('Age (postnatal days)');
    ylabel([titles{i}, ' density (counts/area)']);
    title([titles{i}, ' Decay: Central vs Periphery']);

    % x ticks
    xticks(xTicks);        % 1,2,3,4,5,6,7
    xticklabels(xLabels);  % 23,45,60,90,120,150,200

    % equations
    eq_c = sprintf('Central: D = exp(%.3f %+.3f·Age)', b0, b1);
    eq_p = sprintf('Periphery: D = exp(%.3f %+.3f·Age)', b0+b2, b1+b3);

    xl = xlim;
    yl = ylim;
    x_text  = xl(1) + 0.04 * range(xl);
    y_text1 = yl(2) - 0.08 * range(yl);
    y_text2 = yl(2) - 0.16 * range(yl);

    text(x_text, y_text1, eq_c, 'Color', 'r', 'FontSize', 9, ...
         'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 2);
    text(x_text, y_text2, eq_p, 'Color', 'b', 'FontSize', 9, ...
         'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 2);

    legend({'Central 95% CI', 'Periph 95% CI', ...
        'Central fit',    'Periph fit'}, ...
       'Location', 'northeast');

    box off; grid on;
    hold off;

   saveas(gcf, ['O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\' sheets{i} ' NEW decay curve + 95% ci.png']);
    exportgraphics(gcf, ...
        ['O:\stemcell\MRC-Subretinal Transplantation Project\Rd10 paper\Figures\ephys supplementary\' ...
         sheets{i} ' NEW decay curve + 95% ci horizontal.pdf']);
    close all;

end


