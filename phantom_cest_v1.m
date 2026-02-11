%% Interactive 24-Tube Phantom Mask Creator
% This script allows you to manually draw ROIs for all 24 tubes
% and saves them in a format compatible with the CEST analysis pipeline
%
% Output: tube_masks_24.mat containing 'tube_masks' (xDim x yDim x 24)

clearvars; clc; close all;

%% Configuration
SAVE_BACKUP = true;  % Save backup after every 6 tubes
OUTPUT_FILE = 'tube_masks_24.mat';
Starting_Directory = pwd;

%% Step 1: Load reference image
disp('=== STEP 1: Select Reference Image ===');
disp('Select the S0 image or any high-contrast CEST image for ROI drawing');

% Option 1: Load from DICOM folder
choice = questdlg('Select image source:', 'Image Source', ...
    'DICOM folder', 'MAT file', 'Image file', 'DICOM folder');

switch choice
    case 'DICOM folder'
        DIR = uigetdir(pwd, 'Select DICOM folder');
        if DIR == 0
            error('No folder selected. Exiting.');
        end
        dicomFiles = dir(fullfile(DIR, '*.dcm'));
        
        % Show list of available images
        fprintf('\nAvailable DICOM images:\n');
        for i = 1:min(10, length(dicomFiles))
            fprintf('  %d: %s\n', i, dicomFiles(i).name);
        end
        
        img_idx = input(sprintf('Select image index (1-%d) [default: 2 for S0]: ', length(dicomFiles)));
        if isempty(img_idx)
            img_idx = 2;  % Default to S0 image
        end
        
        reference_image = double(dicomread(fullfile(DIR, dicomFiles(img_idx).name)));
        fprintf('✓ Loaded: %s\n', dicomFiles(img_idx).name);
        
    case 'MAT file'
        [file, path] = uigetfile('*.mat', 'Select MAT file with image data');
        if file == 0
            error('No file selected. Exiting.');
        end
        data = load(fullfile(path, file));
        
        % Try to find image variable
        vars = fieldnames(data);
        fprintf('Available variables: %s\n', strjoin(vars, ', '));
        var_name = input('Enter variable name for reference image: ', 's');
        reference_image = double(data.(var_name));
        
        % If 3D, select slice
        if ndims(reference_image) == 3
            slice_idx = input(sprintf('Select slice (1-%d): ', size(reference_image, 3)));
            reference_image = reference_image(:,:,slice_idx);
        end
        
    case 'Image file'
        [file, path] = uigetfile({'*.png;*.jpg;*.tif', 'Image Files'}, 'Select image file');
        if file == 0
            error('No file selected. Exiting.');
        end
        reference_image = double(imread(fullfile(path, file)));
        if ndims(reference_image) == 3
            reference_image = rgb2gray(reference_image);
        end
        
    otherwise
        error('No selection made. Exiting.');
end

[xDim, yDim] = size(reference_image);
fprintf('✓ Image dimensions: %d x %d\n', xDim, yDim);

%% Step 2: Enhance image for better visualization
disp('=== STEP 2: Image Enhancement ===');

% Normalize to [0, 1]
ref_img_norm = (reference_image - min(reference_image(:))) / ...
               (max(reference_image(:)) - min(reference_image(:)));

% Apply contrast enhancement
ref_img_enhanced = imadjust(ref_img_norm);

% Option to apply edge detection overlay
enhance_choice = questdlg('Apply edge detection overlay?', 'Enhancement', ...
    'Yes', 'No', 'Preview both', 'No');

switch enhance_choice
    case 'Yes'
        edges = edge(ref_img_enhanced, 'Canny');
        ref_img_display = ref_img_enhanced;
        ref_img_display(edges) = 1;  % Highlight edges
        
    case 'Preview both'
        figure('Position', [100, 100, 1200, 500]);
        subplot(1,2,1);
        imagesc(ref_img_enhanced); axis image; colormap gray; colorbar;
        title('Original Enhanced', 'FontSize', 12);
        
        subplot(1,2,2);
        edges = edge(ref_img_enhanced, 'Canny');
        temp = ref_img_enhanced;
        temp(edges) = 1;
        imagesc(temp); axis image; colormap gray; colorbar;
        title('With Edge Detection', 'FontSize', 12);
        
        use_edges = questdlg('Use edge-enhanced version?', 'Choice', 'Yes', 'No', 'No');
        close(gcf);
        
        if strcmp(use_edges, 'Yes')
            ref_img_display = temp;
        else
            ref_img_display = ref_img_enhanced;
        end
        
    otherwise
        ref_img_display = ref_img_enhanced;
end

%% Step 3: Define phantom layout
disp('=== STEP 3: Phantom Layout Configuration ===');
fprintf('\nYour phantom has 24 tubes with:\n');
fprintf('  - 3 Iopamidol @ 20mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 Iopamidol @ 50mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 Creatine @ 20mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 Creatine @ 50mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 Taurine @ 20mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 Taurine @ 50mM (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 PLL 0.1%% (pH 6.2, 6.8, 7.4)\n');
fprintf('  - 3 PBS blank (pH 6.2, 6.8, 7.4)\n');

% Define tube labels
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

% Initialize storage
tube_masks = false(xDim, yDim, 24);
tube_info = struct();

%% Step 4: Drawing mode selection
disp('=== STEP 4: ROI Drawing Setup ===');
fprintf('\nDrawing options:\n');
fprintf('  1. Circle ROI (best for round tubes)\n');
fprintf('  2. Polygon ROI (for irregular shapes)\n');
fprintf('  3. Freehand ROI (for complex shapes)\n');

draw_method = input('Select method (1-3) [default: 1]: ');
if isempty(draw_method)
    draw_method = 1;
end

drawing_methods = {'drawcircle', 'drawpolygon', 'drawfreehand'};
selected_method = drawing_methods{draw_method};

fprintf('\n✓ Using %s for ROI drawing\n', selected_method);

% Ask if user wants to see previously drawn ROIs
show_previous = questdlg('Show previously drawn ROIs while drawing new ones?', ...
    'Display Option', 'Yes', 'No', 'Yes');

%% Step 5: Interactive ROI drawing
disp('=== STEP 5: Drawing 24 Tube ROIs ===');
fprintf('\nInstructions:\n');
fprintf('  1. Draw ROI around the tube\n');
fprintf('  2. Double-click inside to finish\n');
fprintf('  3. Review and confirm (or redraw)\n');
fprintf('  4. Progress saves every 6 tubes\n\n');

% Create figure
main_fig = figure('Position', [50, 50, 1400, 800], 'Name', 'Tube Mask Creator');

for tube_idx = 1:24
    clf(main_fig);
    
    % Display image with previous masks
    subplot(1, 2, 1);
    imagesc(ref_img_display); 
    axis image; 
    colormap gray; 
    colorbar;
    hold on;
    
    % Overlay previous masks
    if strcmp(show_previous, 'Yes') && tube_idx > 1
        colors = jet(24);
        for prev = 1:tube_idx-1
            contour(tube_masks(:,:,prev), 1, 'Color', colors(prev,:), 'LineWidth', 1.5);
            stats = regionprops(tube_masks(:,:,prev), 'Centroid');
            if ~isempty(stats)
                text(stats.Centroid(1), stats.Centroid(2), num2str(prev), ...
                    'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center');
            end
        end
    end
    
    title(sprintf('TUBE %d/%d: %s', tube_idx, 24, tube_labels{tube_idx}), ...
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
                tube_masks(:,:,tube_idx) = temp_mask;
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
    
    % Save backup every 6 tubes
    if SAVE_BACKUP && mod(tube_idx, 6) == 0
        backup_file = sprintf('tube_masks_backup_%d.mat', tube_idx);
        save(fullfile(Starting_Directory, backup_file), 'tube_masks', 'tube_info', 'tube_labels');
        fprintf('  💾 Backup saved: %s\n', backup_file);
    end
end

close(main_fig);

%% Step 6: Review all masks
disp('=== STEP 6: Review All Masks ===');

review_fig = figure('Position', [100, 100, 1400, 900]);

% Overview with all masks
subplot(2, 3, [1 2 4 5]);
imagesc(ref_img_display); 
axis image; 
colormap gray; 
hold on;

colors = jet(24);
for t = 1:24
    if any(tube_masks(:,:,t), 'all')
        contour(tube_masks(:,:,t), 1, 'Color', colors(t,:), 'LineWidth', 2);
        if isfield(tube_info, 'centroid') && length(tube_info) >= t
            text(tube_info(t).centroid(1), tube_info(t).centroid(2), num2str(t), ...
                'Color', 'yellow', 'FontSize', 11, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'black');
        end
    end
end
title('All 24 Tube Masks', 'FontSize', 16, 'FontWeight', 'bold');
hold off;

% Statistics panel
subplot(2, 3, 3);
axis off;
text(0.1, 0.95, 'Mask Statistics:', 'FontSize', 12, 'FontWeight', 'bold');
num_valid = sum(any(reshape(tube_masks, [], 24), 1));
text(0.1, 0.85, sprintf('Valid masks: %d/24', num_valid), 'FontSize', 10);

if num_valid > 0
    all_pixels = arrayfun(@(x) x.pixels, tube_info);
    text(0.1, 0.75, sprintf('Avg pixels: %.0f', mean(all_pixels)), 'FontSize', 10);
    text(0.1, 0.65, sprintf('Min pixels: %d', min(all_pixels)), 'FontSize', 10);
    text(0.1, 0.55, sprintf('Max pixels: %d', max(all_pixels)), 'FontSize', 10);
end

% List problematic tubes
if num_valid < 24
    missing = find(~any(reshape(tube_masks, [], 24), 1));
    text(0.1, 0.40, 'Missing tubes:', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    y_pos = 0.35;
    for m = missing
        text(0.1, y_pos, sprintf('  Tube %d', m), 'FontSize', 9, 'Color', 'red');
        y_pos = y_pos - 0.05;
    end
end

% Colorbar legend
subplot(2, 3, 6);
imagesc(reshape(1:24, 6, 4)');
colormap(gca, jet);
colorbar('Ticks', 1:24);
title('Tube Index', 'FontSize', 11);
axis off;

sgtitle('Final Mask Review', 'FontSize', 16, 'FontWeight', 'bold');

%% Step 7: Quality checks
disp('=== STEP 7: Quality Checks ===');
fprintf('\n--- Quality Report ---\n');

% Check 1: All masks present
if num_valid == 24
    fprintf('✓ All 24 masks created\n');
else
    fprintf('⚠ Only %d/24 masks created\n', num_valid);
end

% Check 2: Overlap detection
overlap_found = false;
for i = 1:23
    for j = i+1:24
        if any(tube_masks(:,:,i) & tube_masks(:,:,j), 'all')
            fprintf('⚠ Overlap detected: Tube %d and Tube %d\n', i, j);
            overlap_found = true;
        end
    end
end
if ~overlap_found
    fprintf('✓ No overlapping masks\n');
end

% Check 3: Size consistency
if num_valid > 0
    pixel_std = std([tube_info.pixels]);
    if pixel_std / mean([tube_info.pixels]) < 0.3
        fprintf('✓ Mask sizes are consistent (CV = %.1f%%)\n', ...
            100 * pixel_std / mean([tube_info.pixels]));
    else
        fprintf('⚠ High variability in mask sizes (CV = %.1f%%)\n', ...
            100 * pixel_std / mean([tube_info.pixels]));
    end
end

%% Step 8: Save final masks
disp('=== STEP 8: Saving Masks ===');

if num_valid < 24
    choice = questdlg(sprintf('Only %d/24 masks created. Save anyway?', num_valid), ...
        'Incomplete Masks', 'Save', 'Cancel', 'Cancel');
    if strcmp(choice, 'Cancel')
        fprintf('⊗ Save cancelled. Masks not saved.\n');
        return;
    end
end

% Save main file
save(fullfile(Starting_Directory, OUTPUT_FILE), 'tube_masks', 'tube_info', ...
    'tube_labels', 'reference_image', 'xDim', 'yDim');
fprintf('✓ Saved: %s\n', OUTPUT_FILE);

% Save figure
saveas(review_fig, fullfile(Starting_Directory, 'tube_masks_overview.png'));
fprintf('✓ Saved: tube_masks_overview.png\n');

% Save tube info as text file
fid = fopen(fullfile(Starting_Directory, 'tube_info.txt'), 'w');
fprintf(fid, 'CEST Phantom - 24 Tube Mask Information\n');
fprintf(fid, 'Created: %s\n\n', datestr(now));
fprintf(fid, '%-5s %-30s %-10s %-20s\n', 'Tube', 'Label', 'Pixels', 'Centroid (x,y)');
fprintf(fid, '%s\n', repmat('-', 1, 70));
for t = 1:24
    if length(tube_info) >= t && isfield(tube_info(t), 'label')
        fprintf(fid, '%-5d %-30s %-10d (%.1f, %.1f)\n', ...
            t, tube_info(t).label, tube_info(t).pixels, ...
            tube_info(t).centroid(1), tube_info(t).centroid(2));
    else
        fprintf(fid, '%-5d %-30s %-10s %s\n', t, 'NOT DEFINED', '-', '-');
    end
end
fclose(fid);
fprintf('✓ Saved: tube_info.txt\n');

%% Step 9: Create individual mask images (optional)
create_individual = questdlg('Create individual mask images for each tube?', ...
    'Individual Masks', 'Yes', 'No', 'No');

if strcmp(create_individual, 'Yes')
    mask_dir = fullfile(Starting_Directory, 'individual_masks');
    if ~exist(mask_dir, 'dir')
        mkdir(mask_dir);
    end
    
    for t = 1:24
        if any(tube_masks(:,:,t), 'all')
            fig_temp = figure('Visible', 'off');
            imagesc(ref_img_display); 
            axis image; 
            colormap gray;
            hold on;
            contour(tube_masks(:,:,t), 1, 'r', 'LineWidth', 2);
            title(sprintf('Tube %d: %s', t, tube_labels{t}), 'FontSize', 12);
            hold off;
            
            saveas(fig_temp, fullfile(mask_dir, sprintf('tube_%02d.png', t)));
            close(fig_temp);
        end
    end
    fprintf('✓ Individual masks saved to: %s\n', mask_dir);
end

%% Step 10: Test loading
disp('=== STEP 10: Testing Saved File ===');
try
    test_data = load(fullfile(Starting_Directory, OUTPUT_FILE));
    assert(isfield(test_data, 'tube_masks'), 'Missing tube_masks variable');
    assert(size(test_data.tube_masks, 3) == 24, 'Incorrect number of masks');
    fprintf('✓ File loads correctly and contains all required data\n');
catch ME
    warning(ME.identifier,'Error testing saved file: %s', ME.message);
end

%% Summary
fprintf('\n========================================\n');
fprintf('MASK CREATION COMPLETE!\n');
fprintf('========================================\n');
fprintf('Output file: %s\n', OUTPUT_FILE);
fprintf('Valid masks: %d/24\n', num_valid);
fprintf('Total pixels: %d\n', sum(tube_masks(:)));
fprintf('\nUse this file with the CEST analysis pipeline.\n');
fprintf('Load with: load(''%s'');\n', OUTPUT_FILE);
fprintf('========================================\n');