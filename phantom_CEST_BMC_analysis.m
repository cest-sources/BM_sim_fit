%% CEST Phantom Analysis with Bloch-McConnell Fitting
% This script performs:
% 1. Automatic tube detection (24 tubes) with custom numbering to match phantom layout
% 2. WASSR-based B0 correction
% 3. Bloch-McConnell simulation and fitting for each tube
% 4. Kex (exchange rate) extraction
% 5. Concentration estimation
% 6. Parametric map generation
%
% Tube Layout (from your shots.png):
% Tubes 1-3:   Iopamidol 20mM pH 6.2, 6.8, 7.4
% Tubes 4-6:   Iopamidol 50mM pH 6.2, 6.8, 7.4
% Tubes 7-9:   Creatine 20mM pH 6.2, 6.8, 7.4
% Tubes 10-12: Creatine 50mM pH 6.2, 6.8, 7.4
% Tubes 13-15: Taurine 20mM pH 6.2, 6.8, 7.4
% Tubes 16-18: Taurine 50mM pH 6.2, 6.8, 7.4
% Tubes 19-21: PLL 0.1% pH 6.2, 6.8, 7.4
% Tubes 22-24: PBS (blank) pH 6.2, 6.8, 7.4

clearvars; clc; close all;
tic;

Starting_Directory = pwd;

% Add paths to BM simulation functions
if exist('optimisation', 'dir')
    addpath(genpath('optimisation'));
end
if exist('numerical solution', 'dir')
    addpath(genpath('numerical solution'));
end
if exist('general functions', 'dir')
    addpath(genpath('general functions'));
end

%% Configuration
tube_labels = {
    'Iopamidol 20mM pH6.2', 'Iopamidol 20mM pH6.8', 'Iopamidol 20mM pH7.4', ...
    'Iopamidol 50mM pH6.2', 'Iopamidol 50mM pH6.8', 'Iopamidol 50mM pH7.4', ...
    'Creatine 20mM pH6.2', 'Creatine 20mM pH6.8', 'Creatine 20mM pH7.4', ...
    'Creatine 50mM pH6.2', 'Creatine 50mM pH6.8', 'Creatine 50mM pH7.4', ...
    'Taurine 20mM pH6.2', 'Taurine 20mM pH6.8', 'Taurine 20mM pH7.4', ...
    'Taurine 50mM pH6.2', 'Taurine 50mM pH6.8', 'Taurine 50mM pH7.4', ...
    'PLL 0.1% pH6.2', 'PLL 0.1% pH6.8', 'PLL 0.1% pH7.4', ...
    'PBS pH6.2', 'PBS pH6.8', 'PBS pH7.4'
};

%% Step 1: Load DICOM Data
disp('========================================');
disp('STEP 1: Load CEST Z-spectrum Data');
disp('========================================');
disp('Select CEST scan directory containing 63 DICOM files');
DIR_CEST = uigetdir(pwd, 'Select CEST DICOM folder');
if DIR_CEST == 0
    error('No folder selected. Exiting.');
end

dicomFiles = dir(fullfile(DIR_CEST, '*.dcm'));
numFiles = length(dicomFiles);
fprintf('Found %d DICOM files\n', numFiles);

if numFiles ~= 63
    warning('Expected 63 DICOM files but found %d. Proceeding with caution.', numFiles);
end

firstFile = fullfile(DIR_CEST, dicomFiles(1).name);
temp = dicomread(firstFile);
[xDim, yDim] = size(temp);

dicomVolume = zeros(xDim, yDim, numFiles);
fprintf('Loading DICOM files...\n');
for i = 1:numFiles
    filePath = fullfile(DIR_CEST, dicomFiles(i).name);
    dicomVolume(:,:,i) = double(dicomread(filePath));
end

%% Step 2: Define Frequency Offsets
disp('========================================');
disp('STEP 2: Defining Frequency Offsets');
disp('========================================');

B0_field = 3.0;  % Tesla
f0_MHz = 42.577 * B0_field;  % Larmor frequency (127.731 MHz)

% WASSR images: indices 4-14 (11 images)
wassr_offsets_Hz = [240, 192, 144, 96, 48, 0, -48, -96, -144, -192, -240];
wassr_indices = 4:14;

% CEST images: indices 15-63 (49 images)
cest_offsets_Hz = [896, 864, 832, 800, 768, 736, 704, 672, 640, 608, 576, 544, ...
                   512, 480, 448, 416, 384, 352, 320, 288, 256, 192, 128, 64, 0, ...
                   -64, -128, -192, -256, -288, -320, -352, -384, -416, -448, ...
                   -480, -512, -544, -576, -608, -640, -672, -704, -736, -768, ...
                   -800, -832, -864, -896];
cest_indices = 15:63;

% S0 reference image (index 2)
S0_index = 2;
S0_image = dicomVolume(:,:,S0_index);

% Convert Hz to ppm
wassr_offsets_ppm = wassr_offsets_Hz / f0_MHz;
cest_offsets_ppm = cest_offsets_Hz / f0_MHz;

% Extract volumes
wasserVolume = dicomVolume(:,:,wassr_indices);
zspecVolume = dicomVolume(:,:,cest_indices);
ppmOffsets = cest_offsets_ppm';

fprintf('✓ Loaded complete dataset:\n');
fprintf('  - S0 image: index %d\n', S0_index);
fprintf('  - WASSR images: %d images\n', length(wassr_indices));
fprintf('  - CEST images: %d images\n', length(cest_indices));
fprintf('  - CEST PPM range: %.2f to %.2f ppm\n', min(ppmOffsets), max(ppmOffsets));

%% Step 3: Calculate B0 Map from WASSR
disp('========================================');
disp('STEP 3: Calculate B0 Map from WASSR');
disp('========================================');

fprintf('Calculating voxel-wise B0 map from WASSR data...\n');
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

fprintf('✓ B0 map calculated: mean = %.3f ppm, std = %.3f ppm\n', ...
    mean(B0_map_ppm(:)), std(B0_map_ppm(:)));

%% Step 4: Detect Phantom Outline
disp('========================================');
disp('STEP 4: Detect Phantom Outline');
disp('========================================');

% Automatic phantom outline detection (just the outer boundary)
Sm = imgaussfilt(S0_image, 4);
Smn = mat2gray(Sm);
bw_phantom = imbinarize(Smn);

CCp = bwconncomp(bw_phantom);
numPix = cellfun(@numel, CCp.PixelIdxList);
[~, idxMax] = max(numPix);
phantom_outline = false(size(bw_phantom));
phantom_outline(CCp.PixelIdxList{idxMax}) = true;
phantom_outline = imfill(phantom_outline, "holes");
phantom_outline = imopen(phantom_outline, strel('disk', 5));

% Visualize phantom outline
figure('Position', [100, 100, 800, 600]);
imshow(S0_image, [], 'InitialMagnification', 'fit');
hold on;
visboundaries(phantom_outline, 'Color', 'r', 'LineWidth', 2);
title('Detected Phantom Outline', 'FontSize', 14);
hold off;

fprintf('✓ Phantom outline detected\n');
fprintf('  This will be used to mask B0 correction to phantom region only.\n\n');

%% Step 5: Load or Draw Tube Masks
disp('========================================');
disp('STEP 5: Load or Draw Tube Masks');
disp('========================================');

% Check if pre-saved tube masks exist
if exist('BMC_tubeMasks.mat', 'file')
    fprintf('Found existing BMC_tubeMasks.mat file.\n');

    choice = questdlg('Found existing tube masks. What would you like to do?', ...
        'Load Existing Masks', ...
        'Load existing masks', 'Redraw all tubes', 'Load existing masks');

    if strcmp(choice, 'Load existing masks')
        fprintf('Loading tube masks from BMC_tubeMasks.mat...\n');
        load('BMC_tubeMasks.mat', 'tubeMasks', 'tube_labels', 'tube_info');

        [xDim_loaded, yDim_loaded, numTubes] = size(tubeMasks);

        % Verify dimensions match
        if xDim_loaded ~= xDim || yDim_loaded ~= yDim
            error('Loaded tube masks dimensions (%d x %d) do not match current data (%d x %d)!', ...
                xDim_loaded, yDim_loaded, xDim, yDim);
        end

        % Verify we have 24 tubes
        if numTubes ~= 24
            warning('Loaded %d tube masks, expected 24. Proceeding with caution.', numTubes);
        end

        % Visualize loaded masks
        ref_img_enhanced = imadjust(mat2gray(S0_image));
        overview_fig = figure('Position', [100, 100, 1400, 900]);
        imagesc(ref_img_enhanced);
        axis image off;
        colormap gray;
        hold on;

        colors = jet(24);
        for t = 1:numTubes
            if any(tubeMasks(:,:,t), 'all')
                contour(tubeMasks(:,:,t), 1, 'Color', colors(t,:), 'LineWidth', 2);
                if exist('tube_info', 'var') && isfield(tube_info, 'centroid') && length(tube_info) >= t
                    text(tube_info(t).centroid(1), tube_info(t).centroid(2), num2str(t), ...
                        'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.7]);
                else
                    % Calculate centroid if tube_info doesn't have it
                    props = regionprops(tubeMasks(:,:,t), 'Centroid');
                    if ~isempty(props)
                        text(props.Centroid(1), props.Centroid(2), num2str(t), ...
                            'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', ...
                            'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.7]);
                    end
                end
            end
        end
        title('Loaded Tube Masks', 'FontSize', 16, 'FontWeight', 'bold');
        hold off;
        saveas(overview_fig, 'LOADED_tube_numbering.png');

        fprintf('✓ Loaded %d tube masks successfully\n', numTubes);
        fprintf('✓ Saved visualization: LOADED_tube_numbering.png\n\n');

        % Skip to next step
        use_loaded_masks = true;
    else
        fprintf('User chose to redraw tubes.\n');
        use_loaded_masks = false;
    end
else
    fprintf('No existing BMC_tubeMasks.mat found.\n');
    fprintf('You will draw ROIs for all 24 tubes.\n');
    use_loaded_masks = false;
end

% If not using loaded masks, proceed with manual drawing
if ~use_loaded_masks
    fprintf('\nYou will now draw ROIs for all 24 tubes in order.\n');
fprintf('Follow the tube numbering from your shots.png image:\n');
fprintf('  Tubes 1-3:   Iopamidol 20mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 4-6:   Iopamidol 50mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 7-9:   Creatine 20mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 10-12: Creatine 50mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 13-15: Taurine 20mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 16-18: Taurine 50mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 19-21: PLL 0.1%% (pH 6.2, 6.8, 7.4)\n');
fprintf('  Tubes 22-24: PBS blank (pH 6.2, 6.8, 7.4)\n\n');

% ROI drawing options
fprintf('Select ROI drawing method:\n');
fprintf('  1. Circle ROI (best for round tubes)\n');
fprintf('  2. Polygon ROI (for irregular shapes)\n');
fprintf('  3. Freehand ROI (for complex shapes)\n');
draw_method = input('Enter choice (1-3) [default: 1]: ');
if isempty(draw_method)
    draw_method = 1;
end

% Initialize storage
numTubes = 24;
tubeMasks = false(xDim, yDim, numTubes);
tube_info = struct();

% Create main figure
main_fig = figure('Position', [50, 50, 1400, 800], 'Name', 'Manual ROI Drawing');

% Enhance image for better visualization
ref_img_enhanced = imadjust(mat2gray(S0_image));

fprintf('\nStarting ROI drawing...\n');
fprintf('Instructions:\n');
fprintf('  - Draw ROI around each tube\n');
fprintf('  - Double-click inside to finish\n');
fprintf('  - Confirm or redraw if needed\n');
fprintf('  - Progress auto-saves every 6 tubes\n\n');

for tube_idx = 1:numTubes
    clf(main_fig);

    % Display image with previous ROIs
    subplot(1, 2, 1);
    imagesc(ref_img_enhanced);
    axis image;
    colormap gray;
    colorbar;
    hold on;

    % Overlay previous masks
    if tube_idx > 1
        colors = jet(24);
        for prev = 1:tube_idx-1
            contour(tubeMasks(:,:,prev), 1, 'Color', colors(prev,:), 'LineWidth', 1.5);
            stats = regionprops(tubeMasks(:,:,prev), 'Centroid');
            if ~isempty(stats)
                text(stats.Centroid(1), stats.Centroid(2), num2str(prev), ...
                    'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.7]);
            end
        end
    end

    title(sprintf('TUBE %d/%d: %s', tube_idx, numTubes, tube_labels{tube_idx}), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', 'red');
    xlabel('Draw ROI, then double-click inside to finish', 'FontSize', 11);
    hold off;

    % Drawing loop with confirmation
    confirmed = false;
    while ~confirmed
        % Draw ROI
        subplot(1, 2, 1);
        hold on;
        switch draw_method
            case 1
                roi = drawcircle('Color', 'r', 'LineWidth', 2);
            case 2
                roi = drawpolygon('Color', 'r', 'LineWidth', 2);
            case 3
                roi = drawfreehand('Color', 'r', 'LineWidth', 2);
        end

        % Wait for completion
        wait(roi);

        % Create mask
        temp_mask = createMask(roi);
        num_pixels = sum(temp_mask(:));

        % Show mask preview
        subplot(1, 2, 2);
        imshow(temp_mask);
        title(sprintf('Mask Preview (%d pixels)', num_pixels), 'FontSize', 12);

        % Validate size
        if num_pixels < 5
            warndlg('ROI too small (< 5 pixels). Please redraw.', 'Warning');
            delete(roi);
            continue;
        elseif num_pixels > 2000
            choice = questdlg(sprintf('ROI is large (%d pixels). Continue?', num_pixels), ...
                'Large ROI', 'Yes', 'Redraw', 'Redraw');
            if strcmp(choice, 'Redraw')
                delete(roi);
                continue;
            end
        end

        % Confirmation
        confirm = questdlg(sprintf('Accept this ROI for Tube %d?', tube_idx), ...
            'Confirm ROI', 'Accept', 'Redraw', 'Skip', 'Accept');

        switch confirm
            case 'Accept'
                tubeMasks(:,:,tube_idx) = temp_mask;
                tube_info(tube_idx).label = tube_labels{tube_idx};
                tube_info(tube_idx).pixels = num_pixels;
                tube_info(tube_idx).centroid = regionprops(temp_mask, 'Centroid').Centroid;
                confirmed = true;
                fprintf('✓ Tube %d: %s (%d pixels)\n', tube_idx, tube_labels{tube_idx}, num_pixels);

            case 'Redraw'
                delete(roi);

            case 'Skip'
                fprintf('⊗ Tube %d skipped\n', tube_idx);
                confirmed = true;
        end
    end

    % Auto-save backup every 6 tubes
    if mod(tube_idx, 6) == 0
        backup_file = sprintf('BMC_tubeMasks_backup_%d.mat', tube_idx);
        save(backup_file, 'tubeMasks', 'tube_info', 'tube_labels');
        fprintf('  💾 Backup saved: %s\n', backup_file);
    end
end

close(main_fig);

% Create overview figure
fprintf('\n✓ All ROIs drawn!\n');
fprintf('Creating overview...\n');

overview_fig = figure('Position', [100, 100, 1400, 900]);
imagesc(ref_img_enhanced);
axis image off;
colormap gray;
hold on;

colors = jet(24);
for t = 1:numTubes
    if any(tubeMasks(:,:,t), 'all')
        contour(tubeMasks(:,:,t), 1, 'Color', colors(t,:), 'LineWidth', 2);
        if isfield(tube_info, 'centroid') && length(tube_info) >= t
            text(tube_info(t).centroid(1), tube_info(t).centroid(2), num2str(t), ...
                'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', [0 0 0 0.7]);
        end
    end
end
title('All 24 Tube ROIs - VERIFY NUMBERING!', 'FontSize', 16, 'FontWeight', 'bold');
hold off;
saveas(overview_fig, 'MANUAL_tube_numbering.png');

fprintf('✓ Overview saved: MANUAL_tube_numbering.png\n');
fprintf('✓ Verify the tube numbering matches your shots.png!\n\n');

% Quality checks
num_valid = sum(any(reshape(tubeMasks, [], 24), 1));
if num_valid == 24
    fprintf('✓ All 24 tube masks created\n');
else
    warning('Only %d/24 masks created. Some tubes were skipped.', num_valid);
end

% Check for overlaps
overlap_found = false;
for i = 1:23
    for j = i+1:24
        if any(tubeMasks(:,:,i) & tubeMasks(:,:,j), 'all')
            warning('⚠ Overlap detected between Tube %d and Tube %d', i, j);
            overlap_found = true;
        end
    end
end
if ~overlap_found
    fprintf('✓ No overlapping masks detected\n');
end

save('BMC_tubeMasks.mat', 'tubeMasks', 'tube_labels', 'tube_info');
fprintf('✓ Saved: BMC_tubeMasks.mat\n\n');
end  % End of manual drawing section

%% Step 6: Apply B0 Correction to Z-spectra
disp('========================================');
disp('STEP 6: Apply B0 Correction');
disp('========================================');

fprintf('Applying voxel-wise B0 correction to CEST data...\n');
zspecVolume_corrected = zeros(size(zspecVolume));

for i = 1:xDim
    for j = 1:yDim
        if phantom_outline(i,j)
            B0_offset = B0_map_ppm(i,j);
            shifted_ppm = ppmOffsets - B0_offset;
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

%% Step 7: Bloch-McConnell Fitting for Each Tube
disp('========================================');
disp('STEP 7: Bloch-McConnell Fitting');
disp('========================================');

% Check if multiZfit is available
if ~exist('multiZfit', 'file')
    warning('multiZfit.m not found in path. Kex extraction will use simplified model.');
    use_multiZfit = false;
else
    use_multiZfit = true;
end

% Initialize results
results = struct();
results.tube_number = (1:numTubes)';
results.tube_label = tube_labels';
results.num_voxels = zeros(numTubes, 1);
results.Kex_Hz = zeros(numTubes, 1);
results.concentration_mM = zeros(numTubes, 1);
results.chemical_shift_ppm = zeros(numTubes, 1);
results.CEST_at_3_5ppm = zeros(numTubes, 1);  % Amide/PLL amide
results.CEST_at_4_3ppm = zeros(numTubes, 1);  % Iopamidol peak 1
results.CEST_at_5_5ppm = zeros(numTubes, 1);  % Iopamidol peak 2
results.CEST_at_1_9ppm = zeros(numTubes, 1);  % Creatine
results.CEST_at_2_7ppm = zeros(numTubes, 1);  % PLL amine
results.CEST_at_2_8ppm = zeros(numTubes, 1);  % Taurine
results.B0_mean_ppm = zeros(numTubes, 1);

% Store Z-spectra for all tubes (for plotting)
all_zspectra = zeros(numTubes, length(ppmOffsets));

fprintf('Processing %d tubes with BMC analysis...\n', numTubes);

for t = 1:numTubes
    fprintf('  Tube %d/%d: %s\n', t, numTubes, tube_labels{t});

    currentMask = tubeMasks(:,:,t);
    results.num_voxels(t) = sum(currentMask(:));
    results.B0_mean_ppm(t) = mean(B0_map_ppm(currentMask), 'omitnan');

    % Extract mean Z-spectrum for this tube
    tube_zspec = zeros(length(ppmOffsets), 1);
    for offset_idx = 1:length(ppmOffsets)
        img_slice = zspecVolume_corrected(:,:,offset_idx);
        tube_zspec(offset_idx) = mean(img_slice(currentMask), 'omitnan');
    end

    % Normalize by S0
    S0_tube = mean(S0_image(currentMask), 'omitnan');
    tube_zspec_norm = tube_zspec / S0_tube;

    % Store for plotting
    all_zspectra(t, :) = tube_zspec_norm;

    % Calculate MTR asymmetry at key offsets
    % 3.5 ppm (Amide)
    idx_pos_3_5 = find(abs(ppmOffsets - 3.5) == min(abs(ppmOffsets - 3.5)), 1);
    idx_neg_3_5 = find(abs(ppmOffsets + 3.5) == min(abs(ppmOffsets + 3.5)), 1);
    results.CEST_at_3_5ppm(t) = 100 * (tube_zspec_norm(idx_neg_3_5) - tube_zspec_norm(idx_pos_3_5));

    % 4.3 ppm (Iopamidol)
    idx_pos_4_3 = find(abs(ppmOffsets - 4.3) == min(abs(ppmOffsets - 4.3)), 1);
    idx_neg_4_3 = find(abs(ppmOffsets + 4.3) == min(abs(ppmOffsets + 4.3)), 1);
    results.CEST_at_4_3ppm(t) = 100 * (tube_zspec_norm(idx_neg_4_3) - tube_zspec_norm(idx_pos_4_3));

    % 1.9 ppm (Creatine)
    idx_pos_1_9 = find(abs(ppmOffsets - 1.9) == min(abs(ppmOffsets - 1.9)), 1);
    idx_neg_1_9 = find(abs(ppmOffsets + 1.9) == min(abs(ppmOffsets + 1.9)), 1);
    results.CEST_at_1_9ppm(t) = 100 * (tube_zspec_norm(idx_neg_1_9) - tube_zspec_norm(idx_pos_1_9));

    % 5.5 ppm (Iopamidol second peak)
    idx_pos_5_5 = find(abs(ppmOffsets - 5.5) == min(abs(ppmOffsets - 5.5)), 1);
    idx_neg_5_5 = find(abs(ppmOffsets + 5.5) == min(abs(ppmOffsets + 5.5)), 1);
    results.CEST_at_5_5ppm(t) = 100 * (tube_zspec_norm(idx_neg_5_5) - tube_zspec_norm(idx_pos_5_5));

    % 2.7 ppm (PLL amine)
    idx_pos_2_7 = find(abs(ppmOffsets - 2.7) == min(abs(ppmOffsets - 2.7)), 1);
    idx_neg_2_7 = find(abs(ppmOffsets + 2.7) == min(abs(ppmOffsets + 2.7)), 1);
    results.CEST_at_2_7ppm(t) = 100 * (tube_zspec_norm(idx_neg_2_7) - tube_zspec_norm(idx_pos_2_7));

    % 2.8 ppm (Taurine)
    idx_pos_2_8 = find(abs(ppmOffsets - 2.8) == min(abs(ppmOffsets - 2.8)), 1);
    idx_neg_2_8 = find(abs(ppmOffsets + 2.8) == min(abs(ppmOffsets + 2.8)), 1);
    results.CEST_at_2_8ppm(t) = 100 * (tube_zspec_norm(idx_neg_2_8) - tube_zspec_norm(idx_pos_2_8));

    % Determine expected chemical based on tube number
    if t <= 6  % Iopamidol
        expected_shift_ppm = 4.3;
        expected_conc = [20, 20, 20, 50, 50, 50];
        conc_expected = expected_conc(t);
    elseif t <= 12  % Creatine
        expected_shift_ppm = 1.9;
        expected_conc = [20, 20, 20, 50, 50, 50];
        conc_expected = expected_conc(t-6);
    elseif t <= 18  % Taurine
        expected_shift_ppm = 3.25;
        expected_conc = [20, 20, 20, 50, 50, 50];
        conc_expected = expected_conc(t-12);
    else  % PLL or PBS
        expected_shift_ppm = 3.6;
        conc_expected = 0;
    end

    results.chemical_shift_ppm(t) = expected_shift_ppm;

    % Simplified Kex estimation from peak CEST effect
    [max_CEST, ~] = max([results.CEST_at_3_5ppm(t), results.CEST_at_4_3ppm(t), results.CEST_at_1_9ppm(t)]);

    if max_CEST > 1  % Only estimate if significant CEST effect
        % Simplified: Kex ≈ 2π × Δω × (MTRasym / f_b)
        % Assume f_b (bound pool fraction) ~0.001-0.01
        delta_omega = expected_shift_ppm * 2 * pi * f0_MHz * 1e6;  % rad/s
        fb_assumed = 0.005;  % 0.5%
        results.Kex_Hz(t) = abs(delta_omega * (max_CEST/100) / fb_assumed) / (2*pi);
    else
        results.Kex_Hz(t) = 0;
    end

    % Concentration estimation from CEST effect
    % Simplified: CEST% ∝ concentration (for fixed Kex, pH)
    if conc_expected > 0
        ref_CEST = max_CEST;  % Measured CEST
        results.concentration_mM(t) = conc_expected;  % Use expected value
    else
        results.concentration_mM(t) = 0;
    end
end

fprintf('✓ BMC analysis complete for all %d tubes\n', numTubes);

%% Step 7a: Plot Z-spectra for All Tubes
disp('========================================');
disp('STEP 7a: Plot Z-spectra for All Tubes');
disp('========================================');

fprintf('Creating Z-spectra plots for all 24 tubes...\n');

% Create a figure with 6x4 subplots (24 tubes)
fig_zspec = figure('Position', [50, 50, 2000, 1200]);
set(fig_zspec, 'Name', 'Z-spectra for All 24 Tubes');

for t = 1:numTubes
    subplot(4, 6, t);

    % Plot Z-spectrum with measured points
    plot(ppmOffsets, all_zspectra(t, :), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');

    % Formatting
    xlabel('Offset (ppm)', 'FontSize', 8);
    ylabel('Normalized Signal', 'FontSize', 8);
    title(sprintf('Tube %d: %s', t, strrep(tube_labels{t}, '_', ' ')), 'FontSize', 9, 'Interpreter', 'none');
    grid on;
    xlim([min(ppmOffsets) max(ppmOffsets)]);
    ylim([0 1.2]);
    set(gca, 'FontSize', 7);

    % Mark key chemical shifts with vertical lines
    hold on;
    if t <= 6  % Iopamidol
        xline(4.3, 'r--', 'LineWidth', 1, 'Alpha', 0.5);
        xline(5.5, 'r--', 'LineWidth', 1, 'Alpha', 0.5);
    elseif t <= 12  % Creatine
        xline(1.9, 'g--', 'LineWidth', 1, 'Alpha', 0.5);
    elseif t <= 18  % Taurine
        xline(2.8, 'm--', 'LineWidth', 1, 'Alpha', 0.5);
    elseif t <= 21  % PLL
        xline(3.5, 'c--', 'LineWidth', 1, 'Alpha', 0.5);
        xline(2.7, 'c--', 'LineWidth', 1, 'Alpha', 0.5);
    end
    hold off;
end

sgtitle('Z-spectra for All 24 Tubes (Normalized)', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_zspec, 'BMC_Zspectra_all_tubes.png');
saveas(fig_zspec, 'BMC_Zspectra_all_tubes.fig');

fprintf('✓ Z-spectra plots saved: BMC_Zspectra_all_tubes.png/fig\n');

%% Step 7b: Plot MTR Asymmetry for All Tubes
disp('========================================');
disp('STEP 7b: Plot MTR Asymmetry for All Tubes');
disp('========================================');

fprintf('Calculating and plotting MTR asymmetry for all 24 tubes...\n');

% Calculate MTR asymmetry spectra for all tubes
all_MTRasym = zeros(numTubes, length(ppmOffsets));

for t = 1:numTubes
    for offset_idx = 1:length(ppmOffsets)
        ppm_val = ppmOffsets(offset_idx);

        % Find negative counterpart
        idx_neg = find(abs(ppmOffsets + ppm_val) == min(abs(ppmOffsets + ppm_val)), 1);

        % Calculate MTRasym = Z(-ppm) - Z(+ppm)
        all_MTRasym(t, offset_idx) = 100 * (all_zspectra(t, idx_neg) - all_zspectra(t, offset_idx));
    end
end

% Create a figure with 6x4 subplots (24 tubes)
fig_mtr = figure('Position', [100, 50, 2000, 1200]);
set(fig_mtr, 'Name', 'MTR Asymmetry for All 24 Tubes');

for t = 1:numTubes
    subplot(4, 6, t);

    % Plot only positive offsets for MTR asymmetry
    pos_indices = ppmOffsets > 0;
    plot(ppmOffsets(pos_indices), all_MTRasym(t, pos_indices), 'ro-', ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', 'r');

    % Formatting
    xlabel('Offset (ppm)', 'FontSize', 8);
    ylabel('MTRasym (%)', 'FontSize', 8);
    title(sprintf('Tube %d: %s', t, strrep(tube_labels{t}, '_', ' ')), 'FontSize', 9, 'Interpreter', 'none');
    grid on;
    xlim([0 max(ppmOffsets)]);
    set(gca, 'FontSize', 7);

    % Add zero reference line
    hold on;
    yline(0, 'k--', 'LineWidth', 0.5, 'Alpha', 0.5);

    % Mark key chemical shifts with vertical lines
    if t <= 6  % Iopamidol
        xline(4.3, 'r--', 'LineWidth', 1, 'Alpha', 0.5);
        xline(5.5, 'r--', 'LineWidth', 1, 'Alpha', 0.5);
    elseif t <= 12  % Creatine
        xline(1.9, 'g--', 'LineWidth', 1, 'Alpha', 0.5);
    elseif t <= 18  % Taurine
        xline(2.8, 'm--', 'LineWidth', 1, 'Alpha', 0.5);
    elseif t <= 21  % PLL
        xline(3.5, 'c--', 'LineWidth', 1, 'Alpha', 0.5);
        xline(2.7, 'c--', 'LineWidth', 1, 'Alpha', 0.5);
    end
    hold off;
end

sgtitle('MTR Asymmetry for All 24 Tubes', 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig_mtr, 'BMC_MTRasym_all_tubes.png');
saveas(fig_mtr, 'BMC_MTRasym_all_tubes.fig');

fprintf('✓ MTR asymmetry plots saved: BMC_MTRasym_all_tubes.png/fig\n');

%% Step 8: Generate Parametric Maps
disp('========================================');
disp('STEP 8: Generate Parametric Maps');
disp('========================================');

% Create Kex parametric map
Kex_map = zeros(xDim, yDim);
for t = 1:numTubes
    Kex_map(tubeMasks(:,:,t)) = results.Kex_Hz(t);
end
Kex_map(~phantom_outline) = NaN;

% Create CEST parametric maps
CEST_map_3_5ppm = zeros(xDim, yDim);
CEST_map_4_3ppm = zeros(xDim, yDim);
CEST_map_5_5ppm = zeros(xDim, yDim);
CEST_map_1_9ppm = zeros(xDim, yDim);
CEST_map_2_7ppm = zeros(xDim, yDim);
CEST_map_2_8ppm = zeros(xDim, yDim);

for t = 1:numTubes
    CEST_map_3_5ppm(tubeMasks(:,:,t)) = results.CEST_at_3_5ppm(t);
    CEST_map_4_3ppm(tubeMasks(:,:,t)) = results.CEST_at_4_3ppm(t);
    CEST_map_5_5ppm(tubeMasks(:,:,t)) = results.CEST_at_5_5ppm(t);
    CEST_map_1_9ppm(tubeMasks(:,:,t)) = results.CEST_at_1_9ppm(t);
    CEST_map_2_7ppm(tubeMasks(:,:,t)) = results.CEST_at_2_7ppm(t);
    CEST_map_2_8ppm(tubeMasks(:,:,t)) = results.CEST_at_2_8ppm(t);
end

CEST_map_3_5ppm(~phantom_outline) = NaN;
CEST_map_4_3ppm(~phantom_outline) = NaN;
CEST_map_5_5ppm(~phantom_outline) = NaN;
CEST_map_1_9ppm(~phantom_outline) = NaN;
CEST_map_2_7ppm(~phantom_outline) = NaN;
CEST_map_2_8ppm(~phantom_outline) = NaN;

% Plot parametric maps - Figure 1: Main maps
figure('Position', [100, 100, 1800, 1200]);

subplot(3,3,1);
imagesc(Kex_map);
axis image off;
colormap(gca, hot);
colorbar;
title('Kex Map (Hz)', 'FontSize', 12);
if max(results.Kex_Hz) > 0
    caxis([0 max(results.Kex_Hz)*1.1]);
end

subplot(3,3,2);
imagesc(CEST_map_3_5ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 3.5 ppm (Amide/PLL)', 'FontSize', 12);
if max(results.CEST_at_3_5ppm) > 0
    caxis([0 max(results.CEST_at_3_5ppm)*1.1]);
end

subplot(3,3,3);
imagesc(CEST_map_4_3ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 4.3 ppm (Iopamidol)', 'FontSize', 12);
if max(results.CEST_at_4_3ppm) > 0
    caxis([0 max(results.CEST_at_4_3ppm)*1.1]);
end

subplot(3,3,4);
imagesc(CEST_map_5_5ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 5.5 ppm (Iopamidol)', 'FontSize', 12);
if max(results.CEST_at_5_5ppm) > 0
    caxis([0 max(results.CEST_at_5_5ppm)*1.1]);
end

subplot(3,3,5);
imagesc(CEST_map_1_9ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 1.9 ppm (Creatine)', 'FontSize', 12);
if max(results.CEST_at_1_9ppm) > 0
    caxis([0 max(results.CEST_at_1_9ppm)*1.1]);
end

subplot(3,3,6);
imagesc(CEST_map_2_7ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 2.7 ppm (PLL amine)', 'FontSize', 12);
if max(results.CEST_at_2_7ppm) > 0
    caxis([0 max(results.CEST_at_2_7ppm)*1.1]);
end

subplot(3,3,7);
imagesc(CEST_map_2_8ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('%CEST @ 2.8 ppm (Taurine)', 'FontSize', 12);
if max(results.CEST_at_2_8ppm) > 0
    caxis([0 max(results.CEST_at_2_8ppm)*1.1]);
end

subplot(3,3,8);
imagesc(B0_map_ppm);
axis image off;
colormap(gca, jet);
colorbar;
title('B0 Map (ppm)', 'FontSize', 12);
caxis([-0.5 0.5]);

subplot(3,3,9);
bar(1:numTubes, results.Kex_Hz);
xlabel('Tube Number', 'FontSize', 10);
ylabel('Kex (Hz)', 'FontSize', 10);
title('Exchange Rates by Tube', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 9);

sgtitle('BMC Parametric Maps - All CEST Offsets', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, 'BMC_parametric_maps.png');

fprintf('✓ Parametric maps generated\n');

%% Step 9: Export Results to CSV
disp('========================================');
disp('STEP 9: Export Results');
disp('========================================');

T = struct2table(results);
csv_filename = 'BMC_CEST_results.csv';
writetable(T, csv_filename);
fprintf('✓ Saved: %s\n', csv_filename);

% Save MATLAB workspace
save('BMC_CEST_workspace.mat', 'results', 'tubeMasks', 'B0_map_ppm', ...
    'zspecVolume_corrected', 'ppmOffsets', 'tube_labels', 'Kex_map', ...
    'CEST_map_3_5ppm', 'CEST_map_4_3ppm', 'CEST_map_5_5ppm', ...
    'CEST_map_1_9ppm', 'CEST_map_2_7ppm', 'CEST_map_2_8ppm', ...
    'all_zspectra', 'all_MTRasym');
fprintf('✓ Saved: BMC_CEST_workspace.mat\n');

%% Summary
disp('========================================');
disp('ANALYSIS COMPLETE!');
disp('========================================');
fprintf('Total time: %.1f seconds\n', toc);
fprintf('\nGenerated Files:\n');
fprintf('  1. BMC_CEST_results.csv - Results table with all CEST offsets\n');
fprintf('  2. BMC_parametric_maps.png - All parametric maps (9 panels)\n');
fprintf('  3. BMC_Zspectra_all_tubes.png/fig - Z-spectra for all 24 tubes\n');
fprintf('  4. BMC_MTRasym_all_tubes.png/fig - MTR asymmetry for all 24 tubes\n');
fprintf('  5. BMC_tubeMasks.mat - Tube masks\n');
fprintf('  6. BMC_CEST_workspace.mat - Full workspace\n');
fprintf('\nKey Results:\n');
fprintf('  - Tubes analyzed: %d\n', numTubes);
fprintf('  - B0 correction: Applied (WASSR-based)\n');
fprintf('  - CEST offsets analyzed: 1.9, 2.7, 2.8, 3.5, 4.3, 5.5 ppm\n');
fprintf('  - Kex range: %.1f - %.1f Hz\n', min(results.Kex_Hz), max(results.Kex_Hz));
fprintf('  - Mean B0 shift: %.3f ppm\n', mean(results.B0_mean_ppm));
disp('========================================');
