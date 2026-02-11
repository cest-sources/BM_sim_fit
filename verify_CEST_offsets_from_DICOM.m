%% DICOM Offset Verification Script for GE 3T Scanner
% This script extracts and verifies frequency offsets from DICOM headers
% Specifically designed for GE clinical scanner CEST acquisitions
%
% Purpose:
%   1. Extract frequency offsets from DICOM headers
%   2. Verify the offsets match your acquisition protocol
%   3. Identify WASSR vs CEST images
%   4. Display offset information in both Hz and ppm
%
% GE DICOM tags for CEST:
%   - (0018,9022) : RF Saturation Frequency Offset (Hz)
%   - (0018,0080) : Repetition Time TR
%   - (0018,0081) : Echo Time TE
%   - (0018,0082) : Inversion Time (if applicable)
%   - (0018,1314) : Flip Angle
%   - (0018,0087) : Magnetic Field Strength

clearvars; clc; close all;

fprintf('========================================\n');
fprintf('CEST DICOM Offset Verification Tool\n');
fprintf('========================================\n\n');

%% Step 1: Select DICOM Directory
fprintf('Step 1: Select DICOM directory\n');
DIR_CEST = uigetdir(pwd, 'Select CEST DICOM folder (all 63 files)');
if DIR_CEST == 0
    error('No folder selected. Exiting.');
end

dicomFiles = dir(fullfile(DIR_CEST, '*.dcm'));
if isempty(dicomFiles)
    dicomFiles = dir(fullfile(DIR_CEST, '*.IMA'));  % Try Siemens format
end

numFiles = length(dicomFiles);
fprintf('Found %d DICOM files\n\n', numFiles);

if numFiles == 0
    error('No DICOM files found in selected directory!');
end

%% Step 2: Extract Offset Information from DICOM Headers
fprintf('Step 2: Extracting offset information from DICOM headers...\n');

% Initialize storage
offsets_Hz = zeros(numFiles, 1);
B0_field = 0;
valid_offset_found = false(numFiles, 1);

% GE-specific DICOM tags
% Standard tag for RF saturation frequency offset
TAG_SAT_OFFSET = '0018,9022';  % RFEchoTrainLength or saturation offset

for i = 1:numFiles
    filePath = fullfile(DIR_CEST, dicomFiles(i).name);

    try
        % Read DICOM header
        info = dicominfo(filePath);

        % Get B0 field strength (only need once)
        if i == 1 && isfield(info, 'MagneticFieldStrength')
            B0_field = info.MagneticFieldStrength;
        end

        % Try multiple methods to extract saturation offset
        offset_found = false;

        % Method 1: Check for standard RF Saturation tag (0018,9022)
        if isfield(info, 'RFEchoTrainLength')
            % Sometimes GE stores offset info here
            offsets_Hz(i) = info.RFEchoTrainLength;
            offset_found = true;
        end

        % Method 2: Check private GE tags
        if ~offset_found
            % GE stores CEST info in private tags starting with 0019, 0043, etc.
            privateFields = fieldnames(info);

            % Look for fields containing 'Private' or specific patterns
            for f = 1:length(privateFields)
                fieldName = privateFields{f};

                % Check if field contains offset information
                if contains(lower(fieldName), {'offset', 'freq', 'saturation'})
                    try
                        val = info.(fieldName);
                        if isnumeric(val) && abs(val) < 10000  % Reasonable range for Hz
                            offsets_Hz(i) = val;
                            offset_found = true;
                            break;
                        end
                    catch
                        % Skip if error reading field
                    end
                end
            end
        end

        % Method 3: Try reading from sequence/protocol name
        if ~offset_found && isfield(info, 'ProtocolName')
            protocolName = info.ProtocolName;
            % Try to parse offset from protocol name (e.g., "CEST_+3.5ppm")
            tokens = regexp(protocolName, '([+-]?\d+\.?\d*)\s*ppm', 'tokens');
            if ~isempty(tokens)
                ppm_val = str2double(tokens{1}{1});
                if B0_field > 0
                    f0_MHz = 42.577 * B0_field;
                    offsets_Hz(i) = ppm_val * f0_MHz;
                    offset_found = true;
                end
            end
        end

        % Method 4: Check SeriesDescription
        if ~offset_found && isfield(info, 'SeriesDescription')
            seriesDesc = info.SeriesDescription;
            tokens = regexp(seriesDesc, '([+-]?\d+\.?\d*)\s*(Hz|ppm)', 'tokens');
            if ~isempty(tokens)
                val = str2double(tokens{1}{1});
                unit = tokens{1}{2};
                if strcmpi(unit, 'ppm') && B0_field > 0
                    f0_MHz = 42.577 * B0_field;
                    offsets_Hz(i) = val * f0_MHz;
                else
                    offsets_Hz(i) = val;
                end
                offset_found = true;
            end
        end

        valid_offset_found(i) = offset_found;

    catch ME
        fprintf('Warning: Could not read file %d: %s\n', i, dicomFiles(i).name);
        fprintf('  Error: %s\n', ME.message);
    end
end

fprintf('Successfully extracted offsets from %d/%d files\n\n', sum(valid_offset_found), numFiles);

%% Step 3: Analyze Extracted Offsets
fprintf('Step 3: Analyzing extracted offsets...\n');

if B0_field == 0
    warning('B0 field strength not found in DICOM headers. Using default 3T.');
    B0_field = 3.0;
end

f0_MHz = 42.577 * B0_field;  % Larmor frequency in MHz
fprintf('Magnetic Field Strength: %.1f T\n', B0_field);
fprintf('Larmor Frequency: %.3f MHz\n\n', f0_MHz);

% Convert to ppm
offsets_ppm = offsets_Hz / f0_MHz;

% Display results
fprintf('========================================\n');
fprintf('EXTRACTED OFFSETS\n');
fprintf('========================================\n');
fprintf('Index | Filename                    | Offset (Hz) | Offset (ppm)\n');
fprintf('------+-----------------------------+-------------+--------------\n');

for i = 1:numFiles
    fprintf('%5d | %-27s | %11.1f | %12.4f\n', ...
        i, dicomFiles(i).name(1:min(27,end)), offsets_Hz(i), offsets_ppm(i));
end

%% Step 4: Identify WASSR vs CEST Images
fprintf('\n========================================\n');
fprintf('OFFSET CLASSIFICATION\n');
fprintf('========================================\n');

% Find unique offsets
[unique_offsets_Hz, ~, idx] = unique(round(offsets_Hz));
unique_offsets_ppm = unique_offsets_Hz / f0_MHz;

fprintf('Number of unique offsets: %d\n', length(unique_offsets_Hz));
fprintf('Offset range: %.1f to %.1f Hz (%.3f to %.3f ppm)\n\n', ...
    min(offsets_Hz), max(offsets_Hz), min(offsets_ppm), max(offsets_ppm));

% Classify as WASSR (small offsets) vs CEST (large offsets)
% WASSR typically: -300 to +300 Hz (±2.35 ppm @ 3T)
% CEST typically: -1000 to +1000 Hz (±7.8 ppm @ 3T)

wassr_threshold_Hz = 250;  % Adjustable threshold
wassr_indices = find(abs(offsets_Hz) <= wassr_threshold_Hz);
cest_indices = find(abs(offsets_Hz) > wassr_threshold_Hz);

fprintf('WASSR images (|offset| <= %.0f Hz):\n', wassr_threshold_Hz);
fprintf('  Indices: %s\n', mat2str(wassr_indices'));
fprintf('  Count: %d\n', length(wassr_indices));
if ~isempty(wassr_indices)
    fprintf('  Offset range: %.1f to %.1f Hz (%.3f to %.3f ppm)\n', ...
        min(offsets_Hz(wassr_indices)), max(offsets_Hz(wassr_indices)), ...
        min(offsets_ppm(wassr_indices)), max(offsets_ppm(wassr_indices)));
end

fprintf('\nCEST images (|offset| > %.0f Hz):\n', wassr_threshold_Hz);
fprintf('  Indices: %s\n', mat2str(cest_indices'));
fprintf('  Count: %d\n', length(cest_indices));
if ~isempty(cest_indices)
    fprintf('  Offset range: %.1f to %.1f Hz (%.3f to %.3f ppm)\n', ...
        min(offsets_Hz(cest_indices)), max(offsets_Hz(cest_indices)), ...
        min(offsets_ppm(cest_indices)), max(offsets_ppm(cest_indices)));
end

%% Step 5: Plot Offset Distribution
fprintf('\n========================================\n');
fprintf('VISUALIZATION\n');
fprintf('========================================\n');

fig = figure('Position', [100, 100, 1400, 800], 'Name', 'CEST Offset Verification');

% Plot 1: Offsets in Hz
subplot(2,3,1);
stem(1:numFiles, offsets_Hz, 'b', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Image Index', 'FontSize', 11);
ylabel('Offset (Hz)', 'FontSize', 11);
title('Saturation Offsets (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold on;
yline(wassr_threshold_Hz, 'r--', 'LineWidth', 1.5, 'Label', 'WASSR/CEST threshold');
yline(-wassr_threshold_Hz, 'r--', 'LineWidth', 1.5);
hold off;

% Plot 2: Offsets in ppm
subplot(2,3,2);
stem(1:numFiles, offsets_ppm, 'r', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Image Index', 'FontSize', 11);
ylabel('Offset (ppm)', 'FontSize', 11);
title('Saturation Offsets (ppm)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
hold on;
yline(wassr_threshold_Hz/f0_MHz, 'k--', 'LineWidth', 1.5);
yline(-wassr_threshold_Hz/f0_MHz, 'k--', 'LineWidth', 1.5);
hold off;

% Plot 3: Histogram of offsets
subplot(2,3,3);
histogram(offsets_ppm, 20, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Offset (ppm)', 'FontSize', 11);
ylabel('Count', 'FontSize', 11);
title('Offset Distribution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 4: WASSR offsets zoomed in
subplot(2,3,4);
if ~isempty(wassr_indices)
    stem(wassr_indices, offsets_Hz(wassr_indices), 'g', 'LineWidth', 1.5, 'MarkerSize', 8);
    xlabel('Image Index', 'FontSize', 11);
    ylabel('Offset (Hz)', 'FontSize', 11);
    title('WASSR Offsets (Zoomed)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
else
    text(0.5, 0.5, 'No WASSR images identified', 'HorizontalAlignment', 'center');
end

% Plot 5: CEST offsets
subplot(2,3,5);
if ~isempty(cest_indices)
    stem(cest_indices, offsets_ppm(cest_indices), 'm', 'LineWidth', 1.5, 'MarkerSize', 8);
    xlabel('Image Index', 'FontSize', 11);
    ylabel('Offset (ppm)', 'FontSize', 11);
    title('CEST Offsets (ppm)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
else
    text(0.5, 0.5, 'No CEST images identified', 'HorizontalAlignment', 'center');
end

% Plot 6: Suggested offset arrays for code
subplot(2,3,6);
axis off;
text(0.1, 0.9, 'SUGGESTED MATLAB CODE:', 'FontSize', 12, 'FontWeight', 'bold');

if ~isempty(wassr_indices)
    wassr_str = sprintf('wassr_indices = %s;\n', mat2str(wassr_indices'));
    wassr_off_str = sprintf('wassr_offsets_Hz = %s;\n', mat2str(round(offsets_Hz(wassr_indices)'), 20));
else
    wassr_str = 'wassr_indices = [];\n';
    wassr_off_str = 'wassr_offsets_Hz = [];\n';
end

if ~isempty(cest_indices)
    cest_str = sprintf('cest_indices = %s;\n', mat2str(cest_indices', 100));
    cest_off_str = sprintf('cest_offsets_Hz = %s;\n', mat2str(round(offsets_Hz(cest_indices)'), 20));
else
    cest_str = 'cest_indices = [];\n';
    cest_off_str = 'cest_offsets_Hz = [];\n';
end

code_text = sprintf(['B0_field = %.1f;  %% Tesla\n' ...
    'f0_MHz = 42.577 * B0_field;  %% %.3f MHz\n\n' ...
    '%s%s\n%s%s'], ...
    B0_field, f0_MHz, wassr_str, wassr_off_str, cest_str, cest_off_str);

text(0.1, 0.7, code_text, 'FontSize', 9, 'FontName', 'Courier', ...
    'VerticalAlignment', 'top', 'Interpreter', 'none');

sgtitle('CEST Offset Verification - GE 3T Scanner', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(fig, 'DICOM_offset_verification.png');
saveas(fig, 'DICOM_offset_verification.fig');
fprintf('✓ Saved: DICOM_offset_verification.png/fig\n');

%% Step 6: Export Results to CSV
fprintf('\nStep 6: Exporting results...\n');

results_table = table((1:numFiles)', offsets_Hz, offsets_ppm, ...
    'VariableNames', {'ImageIndex', 'Offset_Hz', 'Offset_ppm'});

% Add filename column
filenames = cell(numFiles, 1);
for i = 1:numFiles
    filenames{i} = dicomFiles(i).name;
end
results_table.Filename = filenames;

% Add classification
classification = cell(numFiles, 1);
classification(wassr_indices) = {'WASSR'};
classification(cest_indices) = {'CEST'};
classification(cellfun(@isempty, classification)) = {'Unknown'};
results_table.Classification = classification;

% Reorder columns
results_table = results_table(:, {'ImageIndex', 'Filename', 'Offset_Hz', 'Offset_ppm', 'Classification'});

writetable(results_table, 'DICOM_offsets_extracted.csv');
fprintf('✓ Saved: DICOM_offsets_extracted.csv\n');

%% Step 7: Summary and Recommendations
fprintf('\n========================================\n');
fprintf('SUMMARY & RECOMMENDATIONS\n');
fprintf('========================================\n');

if sum(valid_offset_found) < numFiles * 0.5
    fprintf('⚠ WARNING: Offsets could not be extracted from most files!\n');
    fprintf('   This is common with GE scanners as offset info may be in private tags.\n');
    fprintf('   Recommendations:\n');
    fprintf('   1. Check your acquisition protocol documentation\n');
    fprintf('   2. Contact your scanner technician for offset information\n');
    fprintf('   3. Use the sequence protocol parameters directly\n');
    fprintf('   4. Try using DICOM viewer to inspect private tags\n\n');
else
    fprintf('✓ Offset extraction successful!\n\n');

    fprintf('Use these arrays in your CEST analysis script:\n\n');
    fprintf('B0_field = %.1f;  %% Tesla\n', B0_field);
    fprintf('f0_MHz = 42.577 * B0_field;  %% %.3f MHz\n\n', f0_MHz);

    if ~isempty(wassr_indices)
        fprintf('%% WASSR images:\n');
        fprintf('wassr_indices = %s;\n', mat2str(wassr_indices'));
        fprintf('wassr_offsets_Hz = [...\n');
        fprintf('    %s];\n\n', sprintf('%.1f, ', offsets_Hz(wassr_indices)));
    end

    if ~isempty(cest_indices)
        fprintf('%% CEST images:\n');
        fprintf('cest_indices = %s;\n', mat2str(cest_indices'));
        fprintf('cest_offsets_Hz = [...\n');
        fprintf('    %s...\n', sprintf('%.0f, ', offsets_Hz(cest_indices)));
        fprintf('    ];\n\n');
    end
end

fprintf('========================================\n');
fprintf('VERIFICATION COMPLETE\n');
fprintf('========================================\n');
