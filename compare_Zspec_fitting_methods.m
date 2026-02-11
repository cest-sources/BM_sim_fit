%% Z-Spectra Fitting Comparison Script
% This script compares different Z-spectra fitting methods to determine
% which provides the best fit for your CEST phantom data
%
% Fitting methods tested:
%   1. Spline interpolation (baseline, no fitting)
%   2. Single Lorentzian fit
%   3. Multi-peak Lorentzian fit (2-5 peaks)
%   4. Polynomial baseline + Multi-Lorentzian
%   5. Bloch-McConnell fit (if multiZfit available)
%
% Purpose:
%   - Visualize different fitting approaches
%   - Calculate goodness-of-fit metrics (R², RMSE, AIC)
%   - Identify which method best captures CEST effects
%   - Generate fitted Z-spectra for all tubes

clearvars; clc; close all;

fprintf('========================================\n');
fprintf('Z-Spectra Fitting Comparison Tool\n');
fprintf('========================================\n\n');

%% Step 1: Load Existing Data
fprintf('Step 1: Loading CEST analysis data...\n');

% Check if workspace from previous analysis exists
if ~exist('BMC_CEST_workspace.mat', 'file')
    error(['BMC_CEST_workspace.mat not found!\n' ...
           'Please run phantom_CEST_BMC_analysis.m first to generate the workspace.']);
end

load('BMC_CEST_workspace.mat', 'all_zspectra', 'ppmOffsets', 'tube_labels', ...
     'results', 'tubeMasks');

numTubes = size(all_zspectra, 1);
numOffsets = length(ppmOffsets);

fprintf('✓ Loaded Z-spectra for %d tubes with %d offsets\n', numTubes, numOffsets);
fprintf('  Offset range: %.2f to %.2f ppm\n\n', min(ppmOffsets), max(ppmOffsets));

%% Step 2: Define Fitting Functions

% Lorentzian function: L(w) = A * Γ² / ((w - w0)² + Γ²) + offset
lorentzian = @(params, w) params(1) ./ (1 + ((w - params(2)) / params(3)).^2) + params(4);

% Multi-Lorentzian (n peaks) + baseline
% params = [A1, w01, Γ1, A2, w02, Γ2, ..., An, w0n, Γn, offset]
multiLorentzian = @(params, w, nPeaks) calcMultiLorentzian(params, w, nPeaks);

% Polynomial baseline (for asymmetric baseline correction)
polyBaseline = @(coeffs, w) polyval(coeffs, w);

fprintf('Step 2: Fitting functions defined\n\n');

%% Step 3: Select Tubes to Analyze
fprintf('Step 3: Select tubes for detailed analysis\n');
fprintf('Options:\n');
fprintf('  1. Analyze all 24 tubes\n');
fprintf('  2. Analyze specific tube types (Iopamidol, Creatine, Taurine, PLL)\n');
fprintf('  3. Select individual tubes\n');

choice = input('Enter choice (1-3) [default: 1]: ');
if isempty(choice)
    choice = 1;
end

switch choice
    case 1
        tubes_to_fit = 1:numTubes;
    case 2
        fprintf('Select tube type:\n');
        fprintf('  1. Iopamidol (tubes 1-6)\n');
        fprintf('  2. Creatine (tubes 7-12)\n');
        fprintf('  3. Taurine (tubes 13-18)\n');
        fprintf('  4. PLL (tubes 19-21)\n');
        fprintf('  5. PBS (tubes 22-24)\n');
        type_choice = input('Enter choice: ');
        switch type_choice
            case 1, tubes_to_fit = 1:6;
            case 2, tubes_to_fit = 7:12;
            case 3, tubes_to_fit = 13:18;
            case 4, tubes_to_fit = 19:21;
            case 5, tubes_to_fit = 22:24;
            otherwise, tubes_to_fit = 1:6;
        end
    case 3
        tubes_to_fit = input('Enter tube numbers (e.g., [1 4 7 10]): ');
    otherwise
        tubes_to_fit = 1:numTubes;
end

fprintf('Selected %d tubes for fitting comparison\n\n', length(tubes_to_fit));

%% Step 4: Fit Each Tube with Different Methods
fprintf('Step 4: Fitting Z-spectra with multiple methods...\n');

% Storage for results
fit_results = struct();
fit_results.tube_idx = tubes_to_fit';
fit_results.methods = {'Spline', 'SingleLorentzian', 'MultiLorentzian_2peak', ...
                       'MultiLorentzian_3peak', 'MultiLorentzian_4peak', ...
                       'PolyBaseline_3peak'};

numMethods = length(fit_results.methods);

% Initialize storage for fitted curves and metrics
fitted_curves = cell(length(tubes_to_fit), numMethods);
R2_values = zeros(length(tubes_to_fit), numMethods);
RMSE_values = zeros(length(tubes_to_fit), numMethods);
AIC_values = zeros(length(tubes_to_fit), numMethods);

% Fitting options
opts = optimset('Display', 'off', 'MaxFunEvals', 5000, 'TolFun', 1e-8);

for tube_idx = 1:length(tubes_to_fit)
    t = tubes_to_fit(tube_idx);
    fprintf('  Fitting tube %d/%d: %s\n', tube_idx, length(tubes_to_fit), tube_labels{t});

    % Get Z-spectrum
    zspec = all_zspectra(t, :)';

    % Remove any NaN or Inf values
    valid_idx = isfinite(zspec) & isfinite(ppmOffsets);
    ppm_clean = ppmOffsets(valid_idx);
    zspec_clean = zspec(valid_idx);

    %% Method 1: Spline Interpolation (Baseline)
    try
        ppm_fine = linspace(min(ppm_clean), max(ppm_clean), 200)';
        zspec_spline = spline(ppm_clean, zspec_clean, ppm_fine);
        fitted_curves{tube_idx, 1} = [ppm_fine, zspec_spline];

        % Calculate metrics (using interpolated values at measured points)
        zspec_interp = interp1(ppm_fine, zspec_spline, ppm_clean);
        [R2_values(tube_idx, 1), RMSE_values(tube_idx, 1), AIC_values(tube_idx, 1)] = ...
            calcFitMetrics(zspec_clean, zspec_interp, 4);  % 4 = spline order approximation
    catch
        fprintf('    ⚠ Spline interpolation failed\n');
    end

    %% Method 2: Single Lorentzian Fit
    try
        % Initial guess: [amplitude, center, width, offset]
        [min_val, min_idx] = min(zspec_clean);
        p0 = [1 - min_val, ppm_clean(min_idx), 1.0, min_val];

        % Bounds
        lb = [0, min(ppm_clean), 0.1, 0];
        ub = [1, max(ppm_clean), 5, 1];

        % Fit
        [popt, resnorm] = lsqcurvefit(lorentzian, p0, ppm_clean, zspec_clean, lb, ub, opts);

        % Generate fitted curve
        ppm_fine = linspace(min(ppm_clean), max(ppm_clean), 200)';
        zspec_fit = lorentzian(popt, ppm_fine);
        fitted_curves{tube_idx, 2} = [ppm_fine, zspec_fit];

        % Calculate metrics
        zspec_fit_measured = lorentzian(popt, ppm_clean);
        [R2_values(tube_idx, 2), RMSE_values(tube_idx, 2), AIC_values(tube_idx, 2)] = ...
            calcFitMetrics(zspec_clean, zspec_fit_measured, 4);
    catch ME
        fprintf('    ⚠ Single Lorentzian fit failed: %s\n', ME.message);
    end

    %% Method 3-5: Multi-Lorentzian Fits (2, 3, 4 peaks)
    for nPeaks = 2:4
        method_idx = 2 + nPeaks - 1;  % Maps to methods 3, 4, 5

        try
            % Initial guess for multi-Lorentzian
            p0_multi = initMultiLorentzianParams(ppm_clean, zspec_clean, nPeaks);

            % Bounds
            [lb_multi, ub_multi] = getMultiLorentzianBounds(ppm_clean, nPeaks);

            % Fit
            fitfunc = @(p, w) multiLorentzian(p, w, nPeaks);
            [popt_multi, resnorm] = lsqcurvefit(fitfunc, p0_multi, ppm_clean, ...
                zspec_clean, lb_multi, ub_multi, opts);

            % Generate fitted curve
            ppm_fine = linspace(min(ppm_clean), max(ppm_clean), 200)';
            zspec_fit = multiLorentzian(popt_multi, ppm_fine, nPeaks);
            fitted_curves{tube_idx, method_idx} = [ppm_fine, zspec_fit];

            % Calculate metrics
            zspec_fit_measured = multiLorentzian(popt_multi, ppm_clean, nPeaks);
            nParams = 3 * nPeaks + 1;
            [R2_values(tube_idx, method_idx), RMSE_values(tube_idx, method_idx), ...
             AIC_values(tube_idx, method_idx)] = calcFitMetrics(zspec_clean, zspec_fit_measured, nParams);

        catch ME
            fprintf('    ⚠ %d-peak Lorentzian fit failed: %s\n', nPeaks, ME.message);
        end
    end

    %% Method 6: Polynomial Baseline + 3-peak Lorentzian
    try
        % Fit polynomial baseline to far offsets
        baseline_idx = abs(ppm_clean) > 5;  % Offsets beyond ±5 ppm
        if sum(baseline_idx) >= 3
            poly_coeffs = polyfit(ppm_clean(baseline_idx), zspec_clean(baseline_idx), 2);
            baseline = polyval(poly_coeffs, ppm_clean);

            % Subtract baseline
            zspec_corrected = zspec_clean - baseline + mean(baseline);

            % Fit 3-peak Lorentzian to corrected spectrum
            p0_multi = initMultiLorentzianParams(ppm_clean, zspec_corrected, 3);
            [lb_multi, ub_multi] = getMultiLorentzianBounds(ppm_clean, 3);

            fitfunc = @(p, w) multiLorentzian(p, w, 3);
            [popt_multi, resnorm] = lsqcurvefit(fitfunc, p0_multi, ppm_clean, ...
                zspec_corrected, lb_multi, ub_multi, opts);

            % Generate fitted curve (add baseline back)
            ppm_fine = linspace(min(ppm_clean), max(ppm_clean), 200)';
            baseline_fine = polyval(poly_coeffs, ppm_fine);
            zspec_fit = multiLorentzian(popt_multi, ppm_fine, 3) - mean(baseline) + baseline_fine;
            fitted_curves{tube_idx, 6} = [ppm_fine, zspec_fit];

            % Calculate metrics
            zspec_fit_measured = multiLorentzian(popt_multi, ppm_clean, 3) - mean(baseline) + baseline;
            nParams = 3 * 3 + 1 + 3;  % 3 Lorentzians + offset + poly coeffs
            [R2_values(tube_idx, 6), RMSE_values(tube_idx, 6), ...
             AIC_values(tube_idx, 6)] = calcFitMetrics(zspec_clean, zspec_fit_measured, nParams);
        end
    catch ME
        fprintf('    ⚠ Polynomial baseline fit failed: %s\n', ME.message);
    end
end

fprintf('✓ Fitting complete!\n\n');

%% Step 5: Compare Fitting Quality
fprintf('Step 5: Comparing fitting quality...\n\n');

% Calculate average metrics across all tubes
avg_R2 = mean(R2_values, 1, 'omitnan');
avg_RMSE = mean(RMSE_values, 1, 'omitnan');
avg_AIC = mean(AIC_values, 1, 'omitnan');

fprintf('========================================\n');
fprintf('AVERAGE FITTING METRICS\n');
fprintf('========================================\n');
fprintf('Method                     | R²      | RMSE    | AIC\n');
fprintf('---------------------------+---------+---------+--------\n');
for m = 1:numMethods
    fprintf('%-26s | %.4f  | %.5f | %.2f\n', ...
        fit_results.methods{m}, avg_R2(m), avg_RMSE(m), avg_AIC(m));
end
fprintf('========================================\n\n');

% Identify best method
[~, best_method_R2] = max(avg_R2);
[~, best_method_AIC] = min(avg_AIC);

fprintf('RECOMMENDATIONS:\n');
fprintf('  Best by R²: %s (R² = %.4f)\n', fit_results.methods{best_method_R2}, avg_R2(best_method_R2));
fprintf('  Best by AIC: %s (AIC = %.2f)\n\n', fit_results.methods{best_method_AIC}, avg_AIC(best_method_AIC));

%% Step 6: Visualize Fits for Selected Tubes
fprintf('Step 6: Generating comparison plots...\n');

% Create comprehensive comparison figure
num_tubes_to_plot = min(6, length(tubes_to_fit));
tube_plot_indices = round(linspace(1, length(tubes_to_fit), num_tubes_to_plot));

fig1 = figure('Position', [50, 50, 1800, 1200], 'Name', 'Fitting Method Comparison');

for plot_idx = 1:num_tubes_to_plot
    tube_idx = tube_plot_indices(plot_idx);
    t = tubes_to_fit(tube_idx);

    subplot(2, 3, plot_idx);

    % Plot measured data
    zspec = all_zspectra(t, :);
    plot(ppmOffsets, zspec, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Measured');
    hold on;

    % Plot different fits
    colors = lines(numMethods);
    for m = 1:numMethods
        if ~isempty(fitted_curves{tube_idx, m})
            plot(fitted_curves{tube_idx, m}(:,1), fitted_curves{tube_idx, m}(:,2), ...
                '-', 'Color', colors(m,:), 'LineWidth', 1.5, 'DisplayName', fit_results.methods{m});
        end
    end

    hold off;
    xlabel('Offset (ppm)', 'FontSize', 9);
    ylabel('Normalized Signal', 'FontSize', 9);
    title(sprintf('Tube %d: %s', t, strrep(tube_labels{t}, '_', ' ')), ...
        'FontSize', 10, 'Interpreter', 'none');
    legend('Location', 'best', 'FontSize', 7);
    grid on;
    xlim([min(ppmOffsets) max(ppmOffsets)]);
    ylim([0 1.2]);
    set(gca, 'FontSize', 8);
end

sgtitle('Z-Spectra Fitting Method Comparison', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, 'Zspec_fitting_comparison.png');
saveas(fig1, 'Zspec_fitting_comparison.fig');
fprintf('✓ Saved: Zspec_fitting_comparison.png/fig\n');

%% Step 7: Metrics Heatmap
fig2 = figure('Position', [100, 100, 1400, 600], 'Name', 'Fitting Metrics');

subplot(1, 3, 1);
imagesc(R2_values');
colorbar;
colormap(hot);
xlabel('Tube Index', 'FontSize', 10);
ylabel('Fitting Method', 'FontSize', 10);
title('R² Values', 'FontSize', 12, 'FontWeight', 'bold');
yticks(1:numMethods);
yticklabels(fit_results.methods);
set(gca, 'FontSize', 8);
caxis([0.9 1]);

subplot(1, 3, 2);
imagesc(RMSE_values');
colorbar;
colormap(hot);
xlabel('Tube Index', 'FontSize', 10);
ylabel('Fitting Method', 'FontSize', 10);
title('RMSE Values', 'FontSize', 12, 'FontWeight', 'bold');
yticks(1:numMethods);
yticklabels(fit_results.methods);
set(gca, 'FontSize', 8);

subplot(1, 3, 3);
imagesc(AIC_values');
colorbar;
colormap(hot);
xlabel('Tube Index', 'FontSize', 10);
ylabel('Fitting Method', 'FontSize', 10);
title('AIC Values (lower is better)', 'FontSize', 12, 'FontWeight', 'bold');
yticks(1:numMethods);
yticklabels(fit_results.methods);
set(gca, 'FontSize', 8);

sgtitle('Fitting Quality Metrics Across All Tubes', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig2, 'Zspec_fitting_metrics.png');
saveas(fig2, 'Zspec_fitting_metrics.fig');
fprintf('✓ Saved: Zspec_fitting_metrics.png/fig\n');

%% Step 8: Export Results
fprintf('\nStep 8: Exporting results...\n');

% Save all fitting results
save('Zspec_fitting_comparison_results.mat', 'fitted_curves', 'R2_values', ...
     'RMSE_values', 'AIC_values', 'fit_results', 'tubes_to_fit');

% Create summary table
summary_table = table();
summary_table.Method = fit_results.methods';
summary_table.Mean_R2 = avg_R2';
summary_table.Mean_RMSE = avg_RMSE';
summary_table.Mean_AIC = avg_AIC';

writetable(summary_table, 'Zspec_fitting_summary.csv');
fprintf('✓ Saved: Zspec_fitting_summary.csv\n');
fprintf('✓ Saved: Zspec_fitting_comparison_results.mat\n');

%% Summary
fprintf('\n========================================\n');
fprintf('FITTING COMPARISON COMPLETE\n');
fprintf('========================================\n');
fprintf('Generated files:\n');
fprintf('  1. Zspec_fitting_comparison.png/fig - Visual comparison\n');
fprintf('  2. Zspec_fitting_metrics.png/fig - Metrics heatmap\n');
fprintf('  3. Zspec_fitting_summary.csv - Average metrics\n');
fprintf('  4. Zspec_fitting_comparison_results.mat - Full results\n\n');

fprintf('Next steps:\n');
fprintf('  1. Review the R² and AIC values to select best fitting method\n');
fprintf('  2. Consider using %s for your analysis\n', fit_results.methods{best_method_AIC});
fprintf('  3. Integrate selected method into phantom_CEST_BMC_analysis.m\n');
fprintf('========================================\n');

%% ========== HELPER FUNCTIONS ==========

function z = calcMultiLorentzian(params, w, nPeaks)
    % Calculate multi-peak Lorentzian
    % params = [A1, w01, Γ1, A2, w02, Γ2, ..., An, w0n, Γn, offset]
    z = params(end) * ones(size(w));  % Start with offset

    for i = 1:nPeaks
        idx = (i-1) * 3 + 1;
        A = params(idx);
        w0 = params(idx + 1);
        gamma = params(idx + 2);
        z = z + A ./ (1 + ((w - w0) / gamma).^2);
    end
end

function p0 = initMultiLorentzianParams(ppm, zspec, nPeaks)
    % Initialize parameters for multi-Lorentzian fit
    p0 = [];

    % Find prominent dips in Z-spectrum
    [pks, locs] = findpeaks(-zspec, 'NPeaks', nPeaks, 'SortStr', 'descend');

    for i = 1:nPeaks
        if i <= length(pks)
            A = -pks(i);              % Amplitude (depth of dip)
            w0 = ppm(locs(i));        % Center position
            gamma = 1.0;              % Width (initial guess)
        else
            % If not enough peaks found, add dummy peak
            A = 0.1;
            w0 = randn() * 3;         % Random position
            gamma = 1.0;
        end
        p0 = [p0, A, w0, gamma];
    end

    % Add offset (baseline)
    p0 = [p0, max(zspec)];
end

function [lb, ub] = getMultiLorentzianBounds(ppm, nPeaks)
    % Set bounds for multi-Lorentzian parameters
    lb = [];
    ub = [];

    for i = 1:nPeaks
        lb = [lb, 0, min(ppm), 0.1];           % [A, w0, γ]
        ub = [ub, 1, max(ppm), 5.0];
    end

    % Offset bounds
    lb = [lb, 0];
    ub = [ub, 1.5];
end

function [R2, RMSE, AIC] = calcFitMetrics(y_measured, y_fitted, nParams)
    % Calculate goodness-of-fit metrics
    % R²: Coefficient of determination
    % RMSE: Root mean square error
    % AIC: Akaike Information Criterion

    % Remove NaN values
    valid = isfinite(y_measured) & isfinite(y_fitted);
    y_measured = y_measured(valid);
    y_fitted = y_fitted(valid);

    n = length(y_measured);

    % Calculate residuals
    residuals = y_measured - y_fitted;
    SS_res = sum(residuals.^2);
    SS_tot = sum((y_measured - mean(y_measured)).^2);

    % R²
    R2 = 1 - SS_res / SS_tot;

    % RMSE
    RMSE = sqrt(SS_res / n);

    % AIC (Akaike Information Criterion)
    % AIC = n * log(SS_res/n) + 2*k
    AIC = n * log(SS_res / n) + 2 * nParams;
end
