%% Diagnose CEST Data Structure
% This script helps you determine the correct:
% - Total number of DICOM files
% - Which indices correspond to S0, WASSR, CEST
% - How many offsets you actually have
%
% This will help fix the "array bounds exceeded" error

clearvars; clc; close all;

fprintf('========================================\n');
fprintf('CEST Data Structure Diagnostic Tool\n');
fprintf('========================================\n\n');

%% Step 1: Load DICOM files and count them
fprintf('Step 1: Analyzing DICOM files...\n');
DIR_CEST = uigetdir(pwd, 'Select your CEST DICOM folder');
if DIR_CEST == 0
    error('No folder selected. Exiting.');
end

dicomFiles = dir(fullfile(DIR_CEST, '*.dcm'));
if isempty(dicomFiles)
    dicomFiles = dir(fullfile(DIR_CEST, '*.IMA'));
end

numFiles = length(dicomFiles);
fprintf('✓ Found %d DICOM files\n\n', numFiles);

%% Step 2: Try to extract image numbers from filenames
fprintf('Step 2: Analyzing file naming pattern...\n');
fprintf('First 10 files:\n');
for i = 1:min(10, numFiles)
    fprintf('  %3d: %s\n', i, dicomFiles(i).name);
end
fprintf('\n');

%% Step 3: Load a few sample images to check intensity patterns
fprintf('Step 3: Loading sample images to identify S0/WASSR/CEST...\n');

% Load first, middle, and last images
firstFile = fullfile(DIR_CEST, dicomFiles(1).name);
temp = dicomread(firstFile);
[xDim, yDim] = size(temp);

sample_indices = [1, 5, 10, 15, 30, 50, min(63, numFiles)];
sample_images = zeros(xDim, yDim, length(sample_indices));
sample_means = zeros(length(sample_indices), 1);

for i = 1:length(sample_indices)
    idx = sample_indices(i);
    if idx <= numFiles
        filePath = fullfile(DIR_CEST, dicomFiles(idx).name);
        sample_images(:,:,i) = double(dicomread(filePath));
        sample_means(i) = mean(sample_images(:,:,i), 'all');
    end
end

fprintf('Sample image mean intensities:\n');
for i = 1:length(sample_indices)
    if sample_indices(i) <= numFiles
        fprintf('  Index %2d: mean = %.1f\n', sample_indices(i), sample_means(i));
    end
end
fprintf('\n');

%% Step 4: Ask user about acquisition
fprintf('========================================\n');
fprintf('Step 4: Enter your acquisition details\n');
fprintf('========================================\n\n');

fprintf('Based on your scanner protocol, please answer:\n\n');

% S0 images
fprintf('How many S0 (reference, no saturation) images do you have?\n');
fprintf('(Common values: 1, 2, or 4)\n');
numS0 = input('Number of S0 images: ');
if isempty(numS0), numS0 = 1; end

% WASSR images
fprintf('\nHow many WASSR images do you have?\n');
fprintf('(Common values: 7, 11, 13, 15)\n');
numWASSR = input('Number of WASSR images: ');
if isempty(numWASSR), numWASSR = 7; end

% Calculate CEST
numCEST = numFiles - numS0 - numWASSR;

fprintf('\n========================================\n');
fprintf('CALCULATED STRUCTURE:\n');
fprintf('========================================\n');
fprintf('Total DICOM files: %d\n', numFiles);
fprintf('S0 images: %d\n', numS0);
fprintf('WASSR images: %d\n', numWASSR);
fprintf('CEST images: %d (calculated)\n', numCEST);
fprintf('========================================\n\n');

if numCEST < 0
    error('Invalid configuration! CEST images cannot be negative. Check your counts.');
end

%% Step 5: Ask about ordering
fprintf('Step 5: What is the order of images in your DICOM files?\n');
fprintf('Options:\n');
fprintf('  1. S0 first, then WASSR, then CEST (e.g., 1-4: S0, 5-11: WASSR, 12-63: CEST)\n');
fprintf('  2. S0 first, then CEST, then WASSR\n');
fprintf('  3. WASSR first, then S0, then CEST\n');
fprintf('  4. Other (you will specify manually)\n');

order_choice = input('Enter choice (1-4) [default: 1]: ');
if isempty(order_choice), order_choice = 1; end

switch order_choice
    case 1
        % S0, WASSR, CEST
        S0_indices = 1:numS0;
        wassr_indices = (numS0+1):(numS0+numWASSR);
        cest_indices = (numS0+numWASSR+1):numFiles;

    case 2
        % S0, CEST, WASSR
        S0_indices = 1:numS0;
        cest_indices = (numS0+1):(numS0+numCEST);
        wassr_indices = (numS0+numCEST+1):numFiles;

    case 3
        % WASSR, S0, CEST
        wassr_indices = 1:numWASSR;
        S0_indices = (numWASSR+1):(numWASSR+numS0);
        cest_indices = (numWASSR+numS0+1):numFiles;

    case 4
        % Manual
        fprintf('\nEnter indices manually:\n');
        S0_indices = input('S0 indices (e.g., [1 2 3]): ');
        wassr_indices = input('WASSR indices (e.g., 4:10): ');
        cest_indices = input('CEST indices (e.g., 11:63): ');
end

% For single S0, pick the first one
if length(S0_indices) == 1
    S0_index = S0_indices(1);
else
    fprintf('\nYou have %d S0 images. Which one to use as reference?\n', length(S0_indices));
    S0_index = input(sprintf('Enter S0 index (%s) [default: %d]: ', mat2str(S0_indices), S0_indices(1)));
    if isempty(S0_index), S0_index = S0_indices(1); end
end

%% Step 6: Display recommended configuration
fprintf('\n========================================\n');
fprintf('RECOMMENDED CONFIGURATION:\n');
fprintf('========================================\n\n');

fprintf('%% Total files: %d\n', numFiles);
fprintf('%% Structure: %d S0 + %d WASSR + %d CEST\n\n', numS0, numWASSR, numCEST);

fprintf('B0_field = 3.0;  %% Tesla\n');
fprintf('f0_MHz = 42.577 * B0_field;  %% %.3f MHz\n\n', 42.577 * 3.0);

fprintf('%% S0 reference image\n');
fprintf('S0_index = %d;\n\n', S0_index);

fprintf('%% WASSR images: indices %s (%d images)\n', mat2str(wassr_indices), length(wassr_indices));
fprintf('wassr_indices = %s;\n', mat2str(wassr_indices));
fprintf('wassr_offsets_Hz = [');
fprintf('%.0f ', linspace(150, -150, length(wassr_indices)));
fprintf('];  %% PLACEHOLDER - adjust based on your protocol\n\n');

fprintf('%% CEST images: indices %s (%d images)\n', mat2str(cest_indices), length(cest_indices));
fprintf('cest_indices = %s;\n', mat2str(cest_indices));
fprintf('cest_offsets_Hz = [');
fprintf('%.0f ', linspace(900, -900, length(cest_indices)));
fprintf('...  %% PLACEHOLDER - you need %d offset values\n', length(cest_indices));
fprintf('                   ];  %% Make sure you have EXACTLY %d values!\n\n', length(cest_indices));

%% Step 7: Visualize image order
fprintf('========================================\n');
fprintf('VISUALIZATION\n');
fprintf('========================================\n\n');

fig = figure('Position', [100, 100, 1400, 800], 'Name', 'CEST Data Structure');

% Plot 1: Index assignment
subplot(2,2,1);
hold on;
if ~isempty(S0_indices)
    for idx = S0_indices
        rectangle('Position', [idx-0.4, 0, 0.8, 1], 'FaceColor', [0.2 0.8 0.2], 'EdgeColor', 'k');
    end
end
for idx = wassr_indices
    rectangle('Position', [idx-0.4, 0, 0.8, 1], 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'k');
end
for idx = cest_indices
    rectangle('Position', [idx-0.4, 0, 0.8, 1], 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k');
end
hold off;
xlim([0 numFiles+1]);
ylim([0 1.5]);
xlabel('DICOM File Index', 'FontSize', 11);
title('Image Type by Index', 'FontSize', 12, 'FontWeight', 'bold');
legend({'S0 (reference)', 'WASSR', 'CEST'}, 'Location', 'north');
set(gca, 'YTick', []);

% Plot 2: Sample images
subplot(2,2,2);
for i = 1:min(6, length(sample_indices))
    if sample_indices(i) <= numFiles
        subplot(2,2,2);
        hold on;
        bar(i, sample_means(i), 'FaceColor', [0.5 0.5 0.5]);
        hold off;
    end
end
xlabel('Sample Index', 'FontSize', 11);
ylabel('Mean Intensity', 'FontSize', 11);
title('Sample Image Intensities', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Plot 3: First few sample images
for i = 1:min(4, length(sample_indices))
    if sample_indices(i) <= numFiles
        subplot(2,4,4+i);
        imagesc(sample_images(:,:,i));
        axis image off;
        colormap gray;
        title(sprintf('Index %d', sample_indices(i)), 'FontSize', 10);
    end
end

sgtitle('CEST Data Structure Analysis', 'FontSize', 14, 'FontWeight', 'bold');

saveas(fig, 'CEST_structure_analysis.png');
fprintf('✓ Saved: CEST_structure_analysis.png\n\n');

%% Step 8: Count offsets in user's current array
fprintf('========================================\n');
fprintf('Step 8: CHECK YOUR OFFSET ARRAYS\n');
fprintf('========================================\n\n');

fprintf('IMPORTANT: You need to define EXACTLY:\n');
fprintf('  - %d WASSR offsets (for wassr_offsets_Hz)\n', length(wassr_indices));
fprintf('  - %d CEST offsets (for cest_offsets_Hz)\n\n', length(cest_indices));

fprintf('If your current arrays have different counts, that is causing your error!\n\n');

fprintf('Example for YOUR data:\n\n');
fprintf('wassr_offsets_Hz = zeros(1, %d);  %% Fill with your actual WASSR offsets\n', length(wassr_indices));
fprintf('cest_offsets_Hz = zeros(1, %d);   %% Fill with your actual CEST offsets\n\n', length(cest_indices));

%% Step 9: Export configuration file
fprintf('Step 9: Exporting configuration...\n');

config_file = 'CEST_configuration.txt';
fid = fopen(config_file, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'CEST DATA STRUCTURE CONFIGURATION\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Date: %s\n\n', datestr(now));
fprintf(fid, 'Total DICOM files: %d\n', numFiles);
fprintf(fid, 'S0 images: %d (indices %s)\n', numS0, mat2str(S0_indices));
fprintf(fid, 'WASSR images: %d (indices %s)\n', numWASSR, mat2str(wassr_indices));
fprintf(fid, 'CEST images: %d (indices %s)\n\n', numCEST, mat2str(cest_indices));

fprintf(fid, 'COPY THIS TO YOUR ANALYSIS SCRIPT:\n');
fprintf(fid, '-----------------------------------\n\n');
fprintf(fid, 'B0_field = 3.0;  %% Tesla\n');
fprintf(fid, 'f0_MHz = 42.577 * B0_field;\n\n');
fprintf(fid, 'S0_index = %d;\n\n', S0_index);
fprintf(fid, 'wassr_indices = %s;\n', mat2str(wassr_indices));
fprintf(fid, 'wassr_offsets_Hz = [...];  %% ADD YOUR %d OFFSETS HERE\n\n', length(wassr_indices));
fprintf(fid, 'cest_indices = %s;\n', mat2str(cest_indices));
fprintf(fid, 'cest_offsets_Hz = [...];   %% ADD YOUR %d OFFSETS HERE\n\n', length(cest_indices));

fprintf(fid, '========================================\n');
fprintf(fid, 'CRITICAL: The number of values in wassr_offsets_Hz\n');
fprintf(fid, 'and cest_offsets_Hz MUST EXACTLY match the lengths\n');
fprintf(fid, 'of wassr_indices (%d) and cest_indices (%d)!\n', length(wassr_indices), length(cest_indices));
fprintf(fid, '========================================\n');

fclose(fid);

fprintf('✓ Saved: %s\n\n', config_file);

%% Summary
fprintf('========================================\n');
fprintf('DIAGNOSIS COMPLETE\n');
fprintf('========================================\n\n');

fprintf('Next steps:\n');
fprintf('1. Copy the configuration from CEST_configuration.txt\n');
fprintf('2. Update your analysis script with correct indices\n');
fprintf('3. Make sure you have EXACTLY %d WASSR offsets\n', length(wassr_indices));
fprintf('4. Make sure you have EXACTLY %d CEST offsets\n', length(cest_indices));
fprintf('5. Re-run your analysis\n\n');

fprintf('Files generated:\n');
fprintf('  - CEST_structure_analysis.png\n');
fprintf('  - CEST_configuration.txt\n');
fprintf('========================================\n');
