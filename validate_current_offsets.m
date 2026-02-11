%% Validate Current CEST Offsets Script
% This script validates if your current offset arrays are correct by:
% 1. Loading DICOM data with your current offsets
% 2. Checking if water saturation appears at 0 ppm (WASSR)
% 3. Checking if Z-spectrum has expected shape
% 4. Allowing manual verification of CEST peaks
%
% If offsets are CORRECT, you should see:
%   - Water saturation dip at exactly 0 ppm in WASSR data
%   - Symmetric Z-spectrum shape (before B0 correction)
%   - CEST peaks at expected chemical shifts
%
% If offsets are WRONG, you'll see:
%   - Water dip shifted from 0 ppm
%   - Asymmetric or distorted Z-spectrum
%   - CEST peaks at wrong positions

clearvars; clc; close all;

fprintf('========================================\n');
fprintf('CEST Offset Validation Tool\n');
fprintf('========================================\n\n');

%% Step 1: Load Your Current Offset Arrays
fprintf('Step 1: Using your current offset definitions...\n\n');

% YOUR CURRENT OFFSETS (from phantom_CEST_BMC_analysis.m)
B0_field = 3.0;  % Tesla
f0_MHz = 42.577 * B0_field;  % Larmor frequency (127.731 MHz)

% WASSR images: indices 4-14 (11 images)
% CORRECTED OFFSETS (validation-based correction: +48 Hz shift applied)
wassr_offsets_Hz = [288, 240, 192, 144, 96, 48, 0, -48, -96, -144, -192];
wassr_indices = 4:14;

% CEST images: indices 15-63 (49 images)
% CORRECTED OFFSETS (validation-based correction: +48 Hz shift applied)
cest_offsets_Hz = [944, 912, 880, 848, 816, 784, 752, 720, 688, 656, 624, 592, ...
                   560, 528, 496, 464, 432, 400, 368, 336, 304, 240, 176, 112, 48, ...
                   -16, -80, -144, -208, -240, -272, -304, -336, -368, -400, -432, ...
                   -464, -496, -528, -560, -592, -624, -656, -688, -720, -752, ...
                   -784, -816, -848];
cest_indices = 15:63;

% S0 reference image (index 2)
S0_index = 2;

% Convert Hz to ppm
wassr_offsets_ppm = wassr_offsets_Hz / f0_MHz;
cest_offsets_ppm = cest_offsets_Hz / f0_MHz;

fprintf('Current offset configuration:\n');
fprintf('  B0 field: %.1f T\n', B0_field);
fprintf('  Larmor frequency: %.3f MHz\n', f0_MHz);
fprintf('  WASSR: %d images (indices %d-%d)\n', length(wassr_indices), min(wassr_indices), max(wassr_indices));
fprintf('  WASSR range: %.1f to %.1f Hz (%.3f to %.3f ppm)\n', ...
    min(wassr_offsets_Hz), max(wassr_offsets_Hz), min(wassr_offsets_ppm), max(wassr_offsets_ppm));
fprintf('  CEST: %d images (indices %d-%d)\n', length(cest_indices), min(cest_indices), max(cest_indices));
fprintf('  CEST range: %.1f to %.1f Hz (%.3f to %.3f ppm)\n', ...
    min(cest_offsets_Hz), max(cest_offsets_Hz), min(cest_offsets_ppm), max(cest_offsets_ppm));
fprintf('  S0 index: %d\n\n', S0_index);

%% Step 2: Load DICOM Data
fprintf('Step 2: Loading DICOM data...\n');
DIR_CEST = uigetdir(pwd, 'Select CEST DICOM folder (all 63 files)');
if DIR_CEST == 0
    error('No folder selected. Exiting.');
end

dicomFiles = dir(fullfile(DIR_CEST, '*.dcm'));
if isempty(dicomFiles)
    dicomFiles = dir(fullfile(DIR_CEST, '*.IMA'));
end

numFiles = length(dicomFiles);
fprintf('Found %d DICOM files\n', numFiles);

if numFiles ~= 63
    warning('Expected 63 files but found %d. Proceeding anyway...', numFiles);
end

% Load all DICOM images
firstFile = fullfile(DIR_CEST, dicomFiles(1).name);
temp = dicomread(firstFile);
[xDim, yDim] = size(temp);

dicomVolume = zeros(xDim, yDim, numFiles);
fprintf('Loading images...\n');
for i = 1:numFiles
    filePath = fullfile(DIR_CEST, dicomFiles(i).name);
    dicomVolume(:,:,i) = double(dicomread(filePath));
    if mod(i, 10) == 0
        fprintf('  Loaded %d/%d\n', i, numFiles);
    end
end
fprintf('✓ All images loaded\n\n');

%% Step 3: Extract Data Using Current Offsets
fprintf('Step 3: Extracting data with current offsets...\n');

S0_image = dicomVolume(:,:,S0_index);
wasserVolume = dicomVolume(:,:,wassr_indices);
zspecVolume = dicomVolume(:,:,cest_indices);

fprintf('✓ Data extracted\n');
fprintf('  S0 image: %dx%d\n', size(S0_image));
fprintf('  WASSR volume: %dx%dx%d\n', size(wasserVolume));
fprintf('  CEST volume: %dx%dx%d\n\n', size(zspecVolume));

%% Step 4: Select ROI for Analysis
fprintf('Step 4: Select a region of interest for validation...\n');
fprintf('This will be used to calculate average Z-spectra.\n\n');

% Display S0 image
fig_roi = figure('Position', [100, 100, 800, 600]);
imagesc(S0_image);
axis image;
colormap gray;
colorbar;
title('S0 Image - Draw ROI in phantom region', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Click to draw circular ROI, then double-click inside to finish');

% Draw ROI
roi = drawcircle('Color', 'r', 'LineWidth', 2);
wait(roi);

% Create mask
roi_mask = createMask(roi);
close(fig_roi);

fprintf('✓ ROI selected (%d pixels)\n\n', sum(roi_mask(:)));

%% Step 5: Calculate Average Z-Spectra
fprintf('Step 5: Calculating average Z-spectra from ROI...\n');

% WASSR spectrum
wassr_spectrum = zeros(length(wassr_offsets_ppm), 1);
for i = 1:length(wassr_offsets_ppm)
    img_slice = wasserVolume(:,:,i);
    wassr_spectrum(i) = mean(img_slice(roi_mask), 'omitnan');
end

% Normalize WASSR by S0
S0_roi = mean(S0_image(roi_mask), 'omitnan');
wassr_spectrum_norm = wassr_spectrum / S0_roi;

% CEST Z-spectrum
cest_spectrum = zeros(length(cest_offsets_ppm), 1);
for i = 1:length(cest_offsets_ppm)
    img_slice = zspecVolume(:,:,i);
    cest_spectrum(i) = mean(img_slice(roi_mask), 'omitnan');
end

% Normalize CEST by S0
cest_spectrum_norm = cest_spectrum / S0_roi;

fprintf('✓ Spectra calculated\n\n');

%% Step 6: Validation Checks
fprintf('========================================\n');
fprintf('VALIDATION CHECKS\n');
fprintf('========================================\n\n');

%% Check 1: WASSR water peak should be at 0 ppm
fprintf('CHECK 1: Water saturation position in WASSR\n');
fprintf('----------------------------------------\n');

[min_wassr, min_wassr_idx] = min(wassr_spectrum_norm);
water_peak_ppm = wassr_offsets_ppm(min_wassr_idx);

fprintf('Water saturation dip found at: %.3f ppm\n', water_peak_ppm);

if abs(water_peak_ppm) < 0.1
    fprintf('✓ PASS: Water peak is at ~0 ppm (offset = %.3f ppm)\n', water_peak_ppm);
    fprintf('  → Your WASSR offsets appear CORRECT\n');
    wassr_correct = true;
else
    fprintf('✗ FAIL: Water peak is at %.3f ppm, should be ~0 ppm\n', water_peak_ppm);
    fprintf('  → Your WASSR offsets may be INCORRECT\n');
    fprintf('  → Offset error: %.3f ppm (%.1f Hz)\n', water_peak_ppm, water_peak_ppm * f0_MHz);
    wassr_correct = false;
end
fprintf('\n');

%% Check 2: Z-spectrum should be symmetric around 0 ppm
fprintf('CHECK 2: Z-spectrum symmetry\n');
fprintf('----------------------------------------\n');

% Calculate MTR asymmetry at far offsets (should be ~0 if symmetric)
far_offset_idx = abs(cest_offsets_ppm) > 5;  % Beyond ±5 ppm
if sum(far_offset_idx) > 0
    % Calculate asymmetry for far offsets
    asym_far = zeros(1, floor(sum(far_offset_idx)/2));
    pos_far = cest_offsets_ppm(far_offset_idx & cest_offsets_ppm > 0);

    for i = 1:length(pos_far)
        ppm_val = pos_far(i);
        idx_pos = find(abs(cest_offsets_ppm - ppm_val) == min(abs(cest_offsets_ppm - ppm_val)), 1);
        idx_neg = find(abs(cest_offsets_ppm + ppm_val) == min(abs(cest_offsets_ppm + ppm_val)), 1);

        if ~isempty(idx_pos) && ~isempty(idx_neg)
            asym_far(i) = abs(cest_spectrum_norm(idx_neg) - cest_spectrum_norm(idx_pos));
        end
    end

    mean_asym = mean(asym_far, 'omitnan');
    fprintf('Mean asymmetry at far offsets (>5 ppm): %.4f\n', mean_asym);

    if mean_asym < 0.05
        fprintf('✓ PASS: Z-spectrum is reasonably symmetric\n');
        fprintf('  → Your CEST offsets appear CORRECT\n');
        cest_correct = true;
    else
        fprintf('⚠ WARNING: Z-spectrum shows asymmetry (%.4f)\n', mean_asym);
        fprintf('  → This could indicate incorrect offsets OR B0 inhomogeneity\n');
        cest_correct = false;
    end
else
    fprintf('⚠ Not enough far offsets to check symmetry\n');
    cest_correct = true;  % Assume OK
end
fprintf('\n');

%% Check 3: Expected CEST peak positions
fprintf('CHECK 3: CEST peak positions\n');
fprintf('----------------------------------------\n');

% Calculate MTR asymmetry
MTRasym = zeros(size(cest_offsets_ppm));
for i = 1:length(cest_offsets_ppm)
    ppm_val = cest_offsets_ppm(i);
    idx_neg = find(abs(cest_offsets_ppm + ppm_val) == min(abs(cest_offsets_ppm + ppm_val)), 1);
    if ~isempty(idx_neg)
        MTRasym(i) = 100 * (cest_spectrum_norm(idx_neg) - cest_spectrum_norm(i));
    end
end

% Find peaks in MTR asymmetry (positive offsets only)
pos_idx = cest_offsets_ppm > 0;
[pks, locs] = findpeaks(MTRasym(pos_idx), 'SortStr', 'descend', 'NPeaks', 5);

if ~isempty(pks)
    fprintf('Detected CEST peaks in MTRasym:\n');
    pos_offsets = cest_offsets_ppm(pos_idx);
    for i = 1:min(3, length(pks))
        fprintf('  Peak %d: %.2f ppm (MTRasym = %.2f%%)\n', i, pos_offsets(locs(i)), pks(i));
    end

    % Check if peaks are at reasonable positions
    peak_ppms = pos_offsets(locs(1:min(3, length(locs))));
    expected_ranges = [1.5 2.5; 2.5 4.0; 4.0 6.0];  % [Creatine; Amide/Taurine; Iopamidol]

    peaks_in_range = false(length(peak_ppms), 1);
    for i = 1:length(peak_ppms)
        for j = 1:size(expected_ranges, 1)
            if peak_ppms(i) >= expected_ranges(j,1) && peak_ppms(i) <= expected_ranges(j,2)
                peaks_in_range(i) = true;
            end
        end
    end

    if sum(peaks_in_range) >= 1
        fprintf('✓ PASS: CEST peaks detected in expected ranges\n');
        fprintf('  → Your CEST offsets appear CORRECT\n');
    else
        fprintf('⚠ WARNING: CEST peaks not in expected ranges\n');
        fprintf('  → Check if peak positions match your phantom contents\n');
    end
else
    fprintf('⚠ No significant CEST peaks detected\n');
    fprintf('  → This could be normal for low-concentration phantoms\n');
end
fprintf('\n');

%% Step 7: Visualization
fprintf('========================================\n');
fprintf('VISUALIZATION\n');
fprintf('========================================\n\n');

fig = figure('Position', [50, 50, 1600, 900], 'Name', 'Offset Validation Results');

%% Plot 1: WASSR spectrum
subplot(2,3,1);
plot(wassr_offsets_ppm, wassr_spectrum_norm, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
xline(0, 'r--', 'LineWidth', 2, 'Label', 'Expected water (0 ppm)');
xline(water_peak_ppm, 'g--', 'LineWidth', 2, 'Label', sprintf('Measured (%.3f ppm)', water_peak_ppm));
hold off;
xlabel('Offset (ppm)', 'FontSize', 11);
ylabel('Normalized Signal', 'FontSize', 11);
title('WASSR: Water Saturation Check', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([min(wassr_offsets_ppm)-0.2, max(wassr_offsets_ppm)+0.2]);
legend('Location', 'best', 'FontSize', 9);

if wassr_correct
    text(0.05, 0.95, '✓ PASS', 'Units', 'normalized', 'FontSize', 14, ...
        'FontWeight', 'bold', 'Color', 'g', 'VerticalAlignment', 'top');
else
    text(0.05, 0.95, '✗ FAIL', 'Units', 'normalized', 'FontSize', 14, ...
        'FontWeight', 'bold', 'Color', 'r', 'VerticalAlignment', 'top');
end

%% Plot 2: WASSR in Hz
subplot(2,3,2);
plot(wassr_offsets_Hz, wassr_spectrum_norm, 'm-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'm');
hold on;
xline(0, 'r--', 'LineWidth', 2);
hold off;
xlabel('Offset (Hz)', 'FontSize', 11);
ylabel('Normalized Signal', 'FontSize', 11);
title('WASSR: Offsets in Hz', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

%% Plot 3: Full CEST Z-spectrum
subplot(2,3,3);
plot(cest_offsets_ppm, cest_spectrum_norm, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'b');
hold on;
xline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
% Mark expected CEST positions
xline(1.9, 'g--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Cr');
xline(2.8, 'm--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Tau');
xline(3.5, 'c--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Amide');
xline(4.3, 'r--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Iop');
xline(5.5, 'r--', 'LineWidth', 1, 'Alpha', 0.5);
hold off;
xlabel('Offset (ppm)', 'FontSize', 11);
ylabel('Normalized Signal', 'FontSize', 11);
title('CEST Z-Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([min(cest_offsets_ppm)-0.5, max(cest_offsets_ppm)+0.5]);
ylim([0 1.2]);
legend('Location', 'best', 'FontSize', 7);

%% Plot 4: MTR Asymmetry
subplot(2,3,4);
pos_idx = cest_offsets_ppm > 0;
plot(cest_offsets_ppm(pos_idx), MTRasym(pos_idx), 'r-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
hold on;
yline(0, 'k--', 'LineWidth', 1);
% Mark expected CEST positions
xline(1.9, 'g--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Cr 1.9');
xline(2.8, 'm--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Tau 2.8');
xline(3.5, 'c--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Amide 3.5');
xline(4.3, 'r--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Iop 4.3');
xline(5.5, 'r--', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Iop 5.5');
hold off;
xlabel('Offset (ppm)', 'FontSize', 11);
ylabel('MTRasym (%)', 'FontSize', 11);
title('MTR Asymmetry (Check Peak Positions)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
xlim([0, max(cest_offsets_ppm)+0.5]);
legend('Location', 'best', 'FontSize', 7);

%% Plot 5: Offset distribution
subplot(2,3,5);
hold on;
stem(wassr_indices, wassr_offsets_Hz, 'b', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'WASSR');
stem(cest_indices, cest_offsets_Hz, 'r', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'CEST');
stem(S0_index, 0, 'g', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'S0');
hold off;
xlabel('DICOM File Index', 'FontSize', 11);
ylabel('Offset (Hz)', 'FontSize', 11);
title('Your Current Offset Assignment', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;

%% Plot 6: Summary text
subplot(2,3,6);
axis off;

summary_text = {
    'VALIDATION SUMMARY',
    '==================',
    '',
    sprintf('1. WASSR water peak: %.3f ppm', water_peak_ppm),
    sprintf('   Status: %s', iif(wassr_correct, '✓ CORRECT', '✗ CHECK NEEDED')),
    '',
    sprintf('2. Z-spectrum symmetry: %.4f', mean_asym),
    sprintf('   Status: %s', iif(cest_correct, '✓ SYMMETRIC', '⚠ CHECK')),
    '',
    '3. CEST peaks detected:',
};

if ~isempty(pks)
    for i = 1:min(3, length(pks))
        summary_text{end+1} = sprintf('   %.2f ppm (%.1f%%)', pos_offsets(locs(i)), pks(i));
    end
else
    summary_text{end+1} = '   (None detected)';
end

summary_text{end+1} = '';
summary_text{end+1} = 'OVERALL VERDICT:';
if wassr_correct && cest_correct
    summary_text{end+1} = '✓✓ OFFSETS APPEAR CORRECT';
    summary_text{end+1} = '';
    summary_text{end+1} = 'You can proceed with analysis!';
    verdict_color = [0 0.6 0];
elseif wassr_correct
    summary_text{end+1} = '⚠ WASSR OK, CEST uncertain';
    summary_text{end+1} = '';
    summary_text{end+1} = 'Check CEST peak positions';
    verdict_color = [1 0.5 0];
else
    summary_text{end+1} = '✗ OFFSETS MAY BE INCORRECT';
    summary_text{end+1} = '';
    summary_text{end+1} = 'Consider adjusting offsets';
    verdict_color = [0.8 0 0];
end

text(0.1, 0.95, strjoin(summary_text, '\n'), 'FontSize', 10, 'FontName', 'Courier', ...
    'VerticalAlignment', 'top', 'Interpreter', 'none');

% Highlight verdict
text(0.1, 0.15, summary_text{end-2}, 'FontSize', 13, 'FontWeight', 'bold', ...
    'Color', verdict_color, 'VerticalAlignment', 'top', 'Interpreter', 'none');

sgtitle('CEST Offset Validation Results', 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(fig, 'Offset_validation_results.png');
saveas(fig, 'Offset_validation_results.fig');
fprintf('✓ Saved: Offset_validation_results.png/fig\n\n');

%% Step 8: Export validation report
fprintf('Step 8: Exporting validation report...\n');

report_file = 'Offset_validation_report.txt';
fid = fopen(report_file, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'CEST OFFSET VALIDATION REPORT\n');
fprintf(fid, '========================================\n\n');

fprintf(fid, 'Date: %s\n\n', datestr(now));

fprintf(fid, 'CURRENT CONFIGURATION:\n');
fprintf(fid, '  B0 field: %.1f T\n', B0_field);
fprintf(fid, '  Larmor frequency: %.3f MHz\n', f0_MHz);
fprintf(fid, '  WASSR indices: %s\n', mat2str(wassr_indices));
fprintf(fid, '  CEST indices: %s\n', mat2str(cest_indices));
fprintf(fid, '  S0 index: %d\n\n', S0_index);

fprintf(fid, 'VALIDATION RESULTS:\n\n');

fprintf(fid, '1. WASSR Water Peak:\n');
fprintf(fid, '   Expected: 0.000 ppm\n');
fprintf(fid, '   Measured: %.3f ppm\n', water_peak_ppm);
fprintf(fid, '   Offset error: %.3f ppm (%.1f Hz)\n', water_peak_ppm, water_peak_ppm * f0_MHz);
fprintf(fid, '   Status: %s\n\n', iif(wassr_correct, 'PASS', 'FAIL'));

fprintf(fid, '2. Z-Spectrum Symmetry:\n');
fprintf(fid, '   Mean asymmetry (>5ppm): %.4f\n', mean_asym);
fprintf(fid, '   Status: %s\n\n', iif(cest_correct, 'PASS', 'WARNING'));

fprintf(fid, '3. CEST Peaks Detected:\n');
if ~isempty(pks)
    for i = 1:min(3, length(pks))
        fprintf(fid, '   Peak %d: %.2f ppm (MTRasym = %.2f%%)\n', i, pos_offsets(locs(i)), pks(i));
    end
else
    fprintf(fid, '   (No significant peaks detected)\n');
end
fprintf(fid, '\n');

fprintf(fid, 'OVERALL VERDICT:\n');
if wassr_correct && cest_correct
    fprintf(fid, '   ✓✓ OFFSETS APPEAR CORRECT\n');
    fprintf(fid, '   You can proceed with your CEST analysis.\n');
elseif wassr_correct
    fprintf(fid, '   ⚠ WASSR offsets OK, CEST uncertain\n');
    fprintf(fid, '   Verify CEST peak positions match your phantom contents.\n');
else
    fprintf(fid, '   ✗ OFFSETS MAY BE INCORRECT\n');
    fprintf(fid, '   Water peak offset by %.1f Hz. Consider adjusting offsets.\n', water_peak_ppm * f0_MHz);
end

fprintf(fid, '\n========================================\n');

fclose(fid);

fprintf('✓ Saved: %s\n\n', report_file);

%% Final message
fprintf('========================================\n');
fprintf('VALIDATION COMPLETE\n');
fprintf('========================================\n\n');

if wassr_correct && cest_correct
    fprintf('✓✓ Your current offsets appear CORRECT!\n');
    fprintf('   Water saturation is at %.3f ppm (expected 0 ppm)\n', water_peak_ppm);
    fprintf('   Z-spectrum shows good symmetry\n\n');
    fprintf('   → You can proceed with phantom_CEST_BMC_analysis.m\n\n');
else
    fprintf('⚠ Your current offsets may need adjustment:\n\n');

    if ~wassr_correct
        fprintf('   Issue 1: Water peak is at %.3f ppm (should be 0 ppm)\n', water_peak_ppm);
        fprintf('   Suggested fix: Add %.1f Hz to all WASSR offsets\n', -water_peak_ppm * f0_MHz);
        fprintf('   Or: Subtract %.3f ppm from all WASSR offsets\n\n', water_peak_ppm);

        fprintf('   Updated WASSR offsets would be:\n');
        fprintf('   wassr_offsets_Hz = %s;\n\n', ...
            mat2str(round(wassr_offsets_Hz - water_peak_ppm * f0_MHz)));
    end

    if ~cest_correct
        fprintf('   Issue 2: Z-spectrum asymmetry detected\n');
        fprintf('   This could be due to:\n');
        fprintf('     - Incorrect CEST offsets\n');
        fprintf('     - B0 field inhomogeneity (normal, corrected later)\n');
        fprintf('     - Strong CEST effects (normal for high concentration)\n\n');
    end
end

fprintf('Files generated:\n');
fprintf('  1. Offset_validation_results.png/fig\n');
fprintf('  2. Offset_validation_report.txt\n');
fprintf('========================================\n');

%% Helper function
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
