%% CEST Phantom Analysis with Automatic Tube Detection - B0 Correction Only
% This script performs:
% 1. Automatic tube detection (24 tubes) using algorithm from phantom_t1.m
% 2. B0 field correction
% 3. %CEST calculation (MTR asymmetry)
% 4. Parametric map generation
% 5. CSV export of results
%
% Author: Automated CEST Analysis Pipeline
% Date: 2025-11-05

clearvars; clc; close all;
tic;

Starting_Directory = pwd;

%% Configuration for 24-tube phantom
% Chemical composition (adjust order after tube detection if needed)
tube_labels = {
    'Iopamidol 20mM pH6.2', 'Iopamidol 20mM pH6.8', 'Iopamidol 20mM pH7.4', ...
    'Iopamidol 50mM pH6.2', 'Iopamidol 50mM pH6.8', 'Iopamidol 50mM pH7.4', ...
    'Creatine 20mM pH6.2', 'Creatine 20mM pH6.8', 'Creatine 20mM pH7.4', ...
    'Creatine 50mM pH6.2', 'Creatine 50mM pH6.8', 'Creatine 50mM pH7.4', ...
    'Taurine 20mM pH6.2', 'Taurine 20mM pH6.8', 'Taurine 20mM pH7.4', ...
    'Taurine 50mM pH6.2', 'Taurine 50mM pH6.8', 'Taurine 50mM pH7.4', ...
    'PLL 0.1% pH6.2', 'PLL 0.1% pH6.8', 'PLL 0.1% pH7.4', ...
    'PBS pH6.2', 'PBS pH6.4', 'PBS pH6.8'
};

% Expected chemical shifts for each compound (ppm) - for reference
chemical_shifts = struct();
chemical_shifts.iopamidol = [4.3, 5.6];  % Amide protons
chemical_shifts.creatine = 1.9;          % Amine protons
chemical_shifts.taurine = 3.25;          % Amine protons
chemical_shifts.PLL = 3.6;               % Lysine amine

%% Step 1: Load CEST Z-spectrum Data
disp('========================================');
disp('STEP 1: Load CEST Z-spectrum Data');
disp('========================================');
disp('Select CEST scan directory containing DICOM files');
DIR_CEST = uigetdir(pwd, 'Select CEST DICOM folder');
if DIR_CEST == 0
    error('No folder selected. Exiting.');
end

% Get list of DICOM files
dicomFiles = dir(fullfile(DIR_CEST, '*.dcm'));
numFiles = length(dicomFiles);
fprintf('Found %d DICOM files\n', numFiles);

if numFiles ~= 63
    warning('Expected 63 DICOM files but found %d. Proceeding with caution.', numFiles);
end

% Read first file to get dimensions
firstFile = fullfile(DIR_CEST, dicomFiles(1).name);
temp = dicomread(firstFile);
[xDim, yDim] = size(temp);

% Preallocate full volume (all 63 images)
dicomVolume = zeros(xDim, yDim, numFiles);

% Load all DICOM files
fprintf('Loading DICOM files...\n');
for i = 1:numFiles
    filePath = fullfile(DIR_CEST, dicomFiles(i).name);
    dicomVolume(:,:,i) = double(dicomread(filePath));
end

%% Step 2: Define Frequency Offsets
disp('========================================');
disp('STEP 2: Defining Frequency Offsets');
disp('========================================');

% Field strength for Hz to ppm conversion
B0_field = 3.0;  % Tesla (change if using different field strength)
f0_MHz = 42.577 * B0_field;  % Larmor frequency for protons (MHz)

% WASSR images: indices 4-14 (11 images)
wassr_offsets_Hz = [240, 192, 144, 96, 48, 0, -48, -96, -144, -192, -240]; % Hz
wassr_indices = 4:14;

% CEST images: indices 15-63 (49 images)
cest_offsets_Hz = [896, 864, 832, 800, 768, 736, 704, 672, 640, 608, 576, 544, ...
                   512, 480, 448, 416, 384, 352, 320, 288, 256, 192, 128, 64, 0, ...
                   -64, -128, -192, -256, -288, -320, -352, -384, -416, -448, ...
                   -480, -512, -544, -576, -608, -640, -672, -704, -736, -768, ...
                   -800, -832, -864, -896]; % Hz
cest_indices = 15:63;

% S0 reference image (index 2)
S0_index = 2;
S0_image = dicomVolume(:,:,S0_index);

% Convert Hz to ppm
wassr_offsets_ppm = wassr_offsets_Hz / f0_MHz;
cest_offsets_ppm = cest_offsets_Hz / f0_MHz;

% Extract CEST volume and ppm offsets
zspecVolume = dicomVolume(:,:,cest_indices);
ppmOffsets = cest_offsets_ppm';

fprintf('✓ Loaded complete dataset:\n');
fprintf('  - S0 image: index %d\n', S0_index);
fprintf('  - WASSR images: %d images (indices %d-%d)\n', length(wassr_indices), wassr_indices(1), wassr_indices(end));
fprintf('  - CEST images: %d images (indices %d-%d)\n', length(cest_indices), cest_indices(1), cest_indices(end));
fprintf('  - CEST PPM range: %.2f to %.2f ppm\n', min(ppmOffsets), max(ppmOffsets));
fprintf('  - Image dimensions: %d x %d\n', xDim, yDim);

%% Step 3: Calculate B0 Map from WASSR Data
disp('========================================');
disp('STEP 3: Calculate B0 Map from WASSR');
disp('========================================');

% Ask user if they want to use WASSR for B0 correction or load external B0 map
choice = questdlg('B0 correction method:', 'B0 Correction', ...
    'Use WASSR data', 'Load external B0 map', 'Skip B0 correction', 'Use WASSR data');

switch choice
    case 'Use WASSR data'
        fprintf('Calculating B0 map from WASSR data...\n');

        % Extract WASSR volume
        wasserVolume = dicomVolume(:,:,wassr_indices);

        % Calculate B0 map from WASSR (find frequency offset with minimum signal)
        B0_map_ppm = zeros(xDim, yDim);

        for i = 1:xDim
            for j = 1:yDim
                wassr_spectrum = squeeze(wasserVolume(i,j,:));

                % Find minimum (water saturation dip)
                [~, min_idx] = min(wassr_spectrum);

                % B0 offset is the ppm where minimum occurs
                B0_map_ppm(i,j) = wassr_offsets_ppm(min_idx);
            end
        end

        % Smooth the B0 map to reduce noise
        B0_map_ppm = imgaussfilt(B0_map_ppm, 2);

        % Clip outliers
        B0_map_ppm(B0_map_ppm < -2) = -2;
        B0_map_ppm(B0_map_ppm > 2) = 2;

        apply_B0_correction = true;
        fprintf('✓ Calculated B0 map from WASSR: mean = %.3f ppm, std = %.3f ppm\n', ...
            mean(B0_map_ppm(:)), std(B0_map_ppm(:)));

    case 'Load external B0 map'
        disp('Select external B0 map DICOM folder');
        DIR_B0 = uigetdir(pwd, 'Select B0 map DICOM folder');

        if DIR_B0 == 0
            warning('No folder selected. Skipping B0 correction.');
            B0_map_ppm = zeros(xDim, yDim);
            apply_B0_correction = false;
        else
            dicomFiles_B0 = dir(fullfile(DIR_B0, '*.dcm'));
            if isempty(dicomFiles_B0)
                warning('No DICOM files found. Skipping B0 correction.');
                B0_map_ppm = zeros(xDim, yDim);
                apply_B0_correction = false;
            else
                B0_Hz = double(dicomread(fullfile(DIR_B0, dicomFiles_B0(1).name)));
                B0_map_ppm = B0_Hz / f0_MHz;

                B0_map_ppm(B0_map_ppm < -2) = -2;
                B0_map_ppm(B0_map_ppm > 2) = 2;

                apply_B0_correction = true;
                fprintf('✓ Loaded external B0 map: mean = %.3f ppm, std = %.3f ppm\n', ...
                    mean(B0_map_ppm(:)), std(B0_map_ppm(:)));
            end
        end

    otherwise  % Skip B0 correction
        disp('⊗ B0 correction skipped');
        B0_map_ppm = zeros(xDim, yDim);
        apply_B0_correction = false;
end

%% Step 4: Automatic Phantom Outline Detection
disp('========================================');
disp('STEP 4: Automatic Phantom Detection');
disp('========================================');

% Use S0 image (index 2)
% S0_image already loaded in Step 2

% Smooth and normalize
Sm = imgaussfilt(S0_image, 4);
Smn = mat2gray(Sm);

% Threshold using Otsu's method
bw_phantom = imbinarize(Smn);

% Keep largest connected component (phantom)
CCp = bwconncomp(bw_phantom);
numPix = cellfun(@numel, CCp.PixelIdxList);
[~, idxMax] = max(numPix);
phantom_outline = false(size(bw_phantom));
phantom_outline(CCp.PixelIdxList{idxMax}) = true;

% Fill holes and smooth
phantom_outline = imfill(phantom_outline, "holes");
phantom_outline = imopen(phantom_outline, strel('disk', 5));

% Visualize
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);
imagesc(S0_image); axis image; colormap gray; colorbar;
title('S0 Image (Reference)', 'FontSize', 12);

subplot(1,3,2);
imagesc(phantom_outline); axis image; colormap gray;
title('Phantom Outline (Binary)', 'FontSize', 12);

subplot(1,3,3);
imshow(S0_image, [], 'InitialMagnification', 'fit');
hold on;
visboundaries(phantom_outline, 'Color','r', 'LineWidth', 2);
title('Overlay: Phantom Boundary', 'FontSize', 12);
hold off;

saveas(gcf, 'CEST_phantom_outline.png');
fprintf('✓ Phantom outline detected\n');

%% Step 5: Automatic Tube Detection (24 tubes)
disp('========================================');
disp('STEP 5: Automatic Tube Detection');
disp('========================================');

% Background subtraction for tube detection
bgd = imgaussfilt(S0_image, 4);
I_sub = S0_image - bgd;

% Binarize tubes
bw_tubes = imbinarize(I_sub);
bw_tubes(~phantom_outline) = 0;  % Mask to phantom interior

% Connected component analysis
CC = bwconncomp(bw_tubes, 8);
numOfPixels = cellfun(@numel, CC.PixelIdxList);

% Remove largest component (background) and small noise
[~, indexOfMax] = max(numOfPixels);
bw_tubes(CC.PixelIdxList{indexOfMax}) = 0;
indexOfSmall = find(numOfPixels < 10);
for i = 1:length(indexOfSmall)
    bw_tubes(CC.PixelIdxList{indexOfSmall(i)}) = 0;
end

% Fill holes
bw_tubes = imfill(bw_tubes, "holes");
bw_tubes = imclearborder(bw_tubes);

% Recompute connected components after cleaning
CC = bwconncomp(bw_tubes, 8);
stats = regionprops(CC, 'Centroid', 'Area');

% Get phantom center for circular ordering
phantomStats = regionprops(phantom_outline, 'Centroid');
phantomCenter = phantomStats.Centroid;

% Get tube centroids
centroids = cat(1, stats.Centroid);

% Compute polar angles relative to phantom center
angles = atan2(centroids(:,2) - phantomCenter(2), ...
               centroids(:,1) - phantomCenter(1));

% Sort by angle (clockwise ordering)
[~, sortedIdx] = sort(angles);

% Check if we detected 24 tubes
numTubes = length(sortedIdx);
if numTubes ~= 24
    warning('Expected 24 tubes but detected %d. Proceeding with detected tubes.', numTubes);
    % Adjust tube labels if necessary
    if numTubes < 24
        tube_labels = tube_labels(1:numTubes);
    end
end

% Create ordered tube masks
tubeMasks = false(xDim, yDim, numTubes);
labeledImage = bwlabel(bw_tubes);

for i = 1:numTubes
    tubeMasks(:,:,i) = (labeledImage == sortedIdx(i));
end

% Visualize detected tubes
figure('Position', [100, 100, 1000, 800]);
imagesc(S0_image);
axis image off;
colormap gray;
hold on;

for i = 1:numTubes
    c = centroids(sortedIdx(i), :);
    text(c(1), c(2), sprintf('%d', i), ...
        'Color', 'red', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.5]);

    % Draw tube boundaries
    visboundaries(tubeMasks(:,:,i), 'Color', 'cyan', 'LineWidth', 1);
end
hold off;
title('Automatic Tube Detection (24 tubes)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'CEST_tube_detection.png');

fprintf('✓ Detected %d tubes with automatic ordering\n', numTubes);

% Save tube masks
save('CEST_tubeMasks.mat', 'tubeMasks', 'tube_labels', 'phantomCenter');

%% Step 6: Apply B0 Correction to Z-spectra
disp('========================================');
disp('STEP 6: B0 Correction');
disp('========================================');

if apply_B0_correction
    fprintf('Applying B0 correction...\n');

    % For each voxel, shift the Z-spectrum by the B0 offset
    zspecVolume_corrected = zeros(size(zspecVolume));

    for i = 1:xDim
        for j = 1:yDim
            if phantom_outline(i,j)
                % Get B0 offset for this voxel
                B0_offset = B0_map_ppm(i,j);

                % Shift ppm offsets by B0 offset
                shifted_ppm = ppmOffsets - B0_offset;

                % Interpolate Z-spectrum to original ppm grid
                original_spectrum = squeeze(zspecVolume(i,j,:));
                corrected_spectrum = interp1(shifted_ppm, original_spectrum, ...
                    ppmOffsets, 'linear', 'extrap');

                zspecVolume_corrected(i,j,:) = corrected_spectrum;
            else
                zspecVolume_corrected(i,j,:) = zspecVolume(i,j,:);
            end
        end
    end

    fprintf('✓ B0 correction applied\n');
    zspecVolume = zspecVolume_corrected;
else
    fprintf('⊗ B0 correction skipped\n');
end

%% Step 7: Calculate %CEST (MTR Asymmetry) for Multiple Offsets
disp('========================================');
disp('STEP 7: Calculate %CEST (MTR Asymmetry)');
disp('========================================');

% Define key CEST offsets to analyze (ppm)
target_offsets = [1.9, 2.0, 3.25, 3.5, 4.3, 5.6];  % Creatine, Taurine, Iopamidol
num_offsets = length(target_offsets);

% Calculate MTR asymmetry maps for each offset
MTR_maps = zeros(xDim, yDim, num_offsets);
actual_offsets = zeros(num_offsets, 1);

for idx = 1:num_offsets
    target_ppm = target_offsets(idx);

    % Find closest positive offset
    pos_idx = find(ppmOffsets > 0);
    [~, closest_idx] = min(abs(ppmOffsets(pos_idx) - target_ppm));
    pos_offset_idx = pos_idx(closest_idx);
    actual_pos_ppm = ppmOffsets(pos_offset_idx);

    % Find corresponding negative offset
    neg_idx = find(ppmOffsets < 0);
    [~, closest_neg_idx] = min(abs(ppmOffsets(neg_idx) + actual_pos_ppm));
    neg_offset_idx = neg_idx(closest_neg_idx);
    actual_neg_ppm = ppmOffsets(neg_offset_idx);

    % Calculate MTR asymmetry: MTR_asym = (Z_neg - Z_pos) / Z_ref
    Z_pos = zspecVolume(:,:,pos_offset_idx);
    Z_neg = zspecVolume(:,:,neg_offset_idx);

    % Use S0 as reference
    Z_ref = S0_image;
    Z_ref(Z_ref == 0) = 1;  % Avoid division by zero

    % MTR asymmetry in percent
    MTR_maps(:,:,idx) = 100 * (Z_neg - Z_pos) ./ Z_ref;
    MTR_maps(:,:,idx) = MTR_maps(:,:,idx) .* phantom_outline;  % Mask outside phantom

    actual_offsets(idx) = actual_pos_ppm;

    fprintf('  Offset %.2f ppm: used ±%.2f ppm (neg: %.2f)\n', ...
        target_ppm, actual_pos_ppm, actual_neg_ppm);
end

fprintf('✓ Calculated %CEST maps for %d offsets\n', num_offsets);

%% Step 8: Extract ROI Statistics for Each Tube
disp('========================================');
disp('STEP 8: Extract Tube-wise Statistics');
disp('========================================');

% Preallocate results table
results = struct();
results.tube_number = (1:numTubes)';
results.tube_label = tube_labels';
results.num_voxels = zeros(numTubes, 1);

% MTR asymmetry means and stds for each offset
for idx = 1:num_offsets
    % Replace periods with underscores in field names (MATLAB doesn't allow dots)
    field_name_mean = sprintf('CEST_%s_ppm_mean', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));
    field_name_std = sprintf('CEST_%s_ppm_std', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));

    results.(field_name_mean) = zeros(numTubes, 1);
    results.(field_name_std) = zeros(numTubes, 1);
end

% Extract statistics for each tube
for t = 1:numTubes
    currentMask = tubeMasks(:,:,t);
    results.num_voxels(t) = sum(currentMask(:));

    for idx = 1:num_offsets
        MTR_values = MTR_maps(:,:,idx);
        tube_values = MTR_values(currentMask);

        % Use same field name format
        field_name_mean = sprintf('CEST_%s_ppm_mean', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));
        field_name_std = sprintf('CEST_%s_ppm_std', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));

        results.(field_name_mean)(t) = mean(tube_values, 'omitnan');
        results.(field_name_std)(t) = std(tube_values, 'omitnan');
    end
end

fprintf('✓ Extracted statistics for %d tubes\n', numTubes);

%% Step 9: Export Results to CSV
disp('========================================');
disp('STEP 9: Export Results');
disp('========================================');

% Convert struct to table
T = struct2table(results);

% Write to CSV
csv_filename = 'CEST_results_B0corrected.csv';
writetable(T, csv_filename);
fprintf('✓ Saved results to: %s\n', csv_filename);

%% Step 10: Generate Parametric Maps
disp('========================================');
disp('STEP 10: Generate Parametric Maps');
disp('========================================');

% Create composite parametric maps for each offset
for idx = 1:num_offsets
    figure('Position', [100, 100, 1400, 500]);

    % Subplot 1: Full parametric map
    subplot(1,3,1);
    MTR_display = MTR_maps(:,:,idx);
    MTR_display(~phantom_outline) = NaN;
    imagesc(MTR_display);
    axis image off;
    colormap(gca, jet);
    caxis([0 10]);  % Typical CEST range 0-10%
    colorbar;
    title(sprintf('CEST Map @ %.2f ppm', actual_offsets(idx)), 'FontSize', 14);

    % Subplot 2: Overlay on anatomical
    subplot(1,3,2);
    imshow(S0_image, []);
    hold on;
    h = imagesc(MTR_display, 'AlphaData', 0.6 * (MTR_display > 0));
    colormap(gca, jet);
    caxis([0 10]);
    title('Overlay on S0 Image', 'FontSize', 14);
    hold off;

    % Subplot 3: Tube-wise bar plot
    subplot(1,3,3);
    field_name_mean = sprintf('CEST_%s_ppm_mean', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));
    field_name_std = sprintf('CEST_%s_ppm_std', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));
    tube_means = results.(field_name_mean);
    tube_stds = results.(field_name_std);
    bar(1:numTubes, tube_means);
    hold on;
    errorbar(1:numTubes, tube_means, tube_stds, 'k.', 'LineWidth', 1);
    xlabel('Tube Number');
    ylabel('%CEST');
    title(sprintf('Tube-wise CEST @ %.2f ppm', actual_offsets(idx)));
    grid on;
    hold off;

    % Save figure
    saveas(gcf, sprintf('CEST_parametric_map_%.2fppm.png', actual_offsets(idx)));
end

fprintf('✓ Generated %d parametric maps\n', num_offsets);

%% Step 11: Summary Plot - All Tubes, All Offsets
disp('========================================');
disp('STEP 11: Generate Summary Plots');
disp('========================================');

figure('Position', [100, 100, 1400, 800]);

% Create heatmap of CEST values (tubes x offsets)
CEST_matrix = zeros(numTubes, num_offsets);
for idx = 1:num_offsets
    field_name_mean = sprintf('CEST_%s_ppm_mean', strrep(sprintf('%.2f', actual_offsets(idx)), '.', '_'));
    CEST_matrix(:, idx) = results.(field_name_mean);
end

imagesc(CEST_matrix');
colormap(jet);
colorbar;
xlabel('Tube Number', 'FontSize', 12);
ylabel('CEST Offset (ppm)', 'FontSize', 12);
yticks(1:num_offsets);
yticklabels(arrayfun(@(x) sprintf('%.2f', x), actual_offsets, 'UniformOutput', false));
title('%CEST Heatmap: All Tubes × All Offsets', 'FontSize', 14, 'FontWeight', 'bold');
caxis([0 10]);

% Add text annotations
for t = 1:numTubes
    for o = 1:num_offsets
        text(t, o, sprintf('%.1f', CEST_matrix(t, o)), ...
            'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 8);
    end
end

saveas(gcf, 'CEST_summary_heatmap.png');

fprintf('✓ Summary plots generated\n');

%% Summary Report
disp('========================================');
disp('ANALYSIS COMPLETE!');
disp('========================================');
fprintf('Total analysis time: %.1f seconds\n', toc);
fprintf('\nGenerated Files:\n');
fprintf('  1. CEST_tubeMasks.mat - Tube mask data\n');
fprintf('  2. %s - Results CSV\n', csv_filename);
fprintf('  3. CEST_phantom_outline.png\n');
fprintf('  4. CEST_tube_detection.png\n');
fprintf('  5. CEST_parametric_map_*.png (×%d)\n', num_offsets);
fprintf('  6. CEST_summary_heatmap.png\n');
fprintf('\nKey Results:\n');
fprintf('  - Number of tubes: %d\n', numTubes);
fprintf('  - Number of CEST offsets analyzed: %d\n', num_offsets);
fprintf('  - B0 correction: %s\n', char(apply_B0_correction * "Applied" + ~apply_B0_correction * "Not Applied"));
fprintf('\nNext steps:\n');
fprintf('  - Review CSV file for quantitative values\n');
fprintf('  - Check parametric maps for quality\n');
fprintf('  - Run full correction version (T1/T2/B1/B0) for final analysis\n');
disp('========================================');
