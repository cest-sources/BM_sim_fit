%% Tube Reordering Script
% Use this script if automatic tube detection numbering doesn't match
% your phantom layout (shots.png)
%
% Instructions:
% 1. Run phantom_CEST_BMC_analysis.m until it creates DETECTED_tube_numbering.png
% 2. Open DETECTED_tube_numbering.png and your shots.png side by side
% 3. Create a mapping array below showing which detected tube corresponds to which actual tube
% 4. Run this script to reorder the tube masks
% 5. Continue with phantom_CEST_BMC_analysis.m from Step 7 onwards

clearvars; clc;

%% Load the automatically detected tube masks
if ~exist('BMC_tubeMasks.mat', 'file')
    error('BMC_tubeMasks.mat not found. Run phantom_CEST_BMC_analysis.m first!');
end

load('BMC_tubeMasks.mat', 'tubeMasks', 'tube_labels', 'phantomCenter', 'centroids', 'sortedIdx');

[xDim, yDim, numTubes] = size(tubeMasks);
fprintf('Loaded %d tube masks (%d x %d)\n', numTubes, xDim, yDim);

%% Define the mapping
% EXAMPLE: If detected tube 5 is actually tube 1, detected 6 is tube 2, etc.
% Create your own mapping by comparing the two images:
%
% tube_order = [5, 6, 7, 8, 1, 2, 3, 4, 9, 10, 11, 12, ...
%               13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24];
%
% This means:
%   Detected tube 1 → becomes new tube 5
%   Detected tube 2 → becomes new tube 6
%   ...and so on

fprintf('\n=== MANUAL TUBE REORDERING ===\n');
fprintf('Compare DETECTED_tube_numbering.png with your shots.png\n');
fprintf('For each position in your shots.png (1-24), identify which detected tube number it is.\n\n');

% Option 1: Interactive input
use_interactive = questdlg('How do you want to provide the mapping?', ...
    'Mapping Method', 'Interactive (type each)', 'Enter as array', 'Load from file', 'Interactive (type each)');

switch use_interactive
    case 'Interactive (type each)'
        tube_order = zeros(1, 24);
        fprintf('\nFor each tube position in your shots.png, enter the detected tube number:\n');
        for i = 1:24
            tube_order(i) = input(sprintf('Tube %d in shots.png is detected tube #: ', i));
        end

    case 'Enter as array'
        fprintf('\nEnter mapping as array, e.g.: [5 6 7 1 2 3 ...]\n');
        tube_order = input('tube_order = ');

    case 'Load from file'
        [file, path] = uigetfile({'*.mat;*.txt', 'Mapping Files'}, 'Select mapping file');
        if file ~= 0
            [~, ~, ext] = fileparts(file);
            if strcmp(ext, '.mat')
                data = load(fullfile(path, file));
                tube_order = data.tube_order;
            else
                tube_order = load(fullfile(path, file));
            end
        else
            error('No file selected');
        end
end

%% Validate mapping
assert(length(tube_order) == 24, 'Mapping must have exactly 24 entries!');
assert(all(ismember(tube_order, 1:24)), 'Mapping must contain numbers 1-24 only!');
assert(length(unique(tube_order)) == 24, 'Mapping has duplicate entries!');

%% Apply reordering
tubeMasks_reordered = tubeMasks(:,:,tube_order);
centroids_reordered = centroids(tube_order, :);

%% Visualize reordered tubes
% Load S0 image (assuming from workspace or file)
if exist('BMC_CEST_workspace.mat', 'file')
    temp = load('BMC_CEST_workspace.mat', 'S0_image');
    S0_image = temp.S0_image;
else
    fprintf('Warning: S0_image not found. Using mask visualization only.\n');
    S0_image = sum(tubeMasks, 3);
end

figure('Position', [100, 100, 1400, 600]);

subplot(1,2,1);
imagesc(S0_image);
axis image off;
colormap gray;
hold on;
for i = 1:numTubes
    c = centroids(i, :);
    text(c(1), c(2), sprintf('%d', i), ...
        'Color', 'yellow', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [1 0 0 0.6]);
end
hold off;
title('BEFORE Reordering (Detected)', 'FontSize', 14, 'FontWeight', 'bold');

subplot(1,2,2);
imagesc(S0_image);
axis image off;
colormap gray;
hold on;
for i = 1:numTubes
    c = centroids_reordered(i, :);
    text(c(1), c(2), sprintf('%d', i), ...
        'Color', 'lime', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', [0 0.5 0 0.6]);
    visboundaries(tubeMasks_reordered(:,:,i), 'Color', 'cyan', 'LineWidth', 1);
end
hold off;
title('AFTER Reordering (Matches shots.png)', 'FontSize', 14, 'FontWeight', 'bold');

sgtitle('Tube Reordering Verification', 'FontSize', 16, 'FontWeight', 'bold');
saveas(gcf, 'REORDERED_tube_numbering.png');

fprintf('\n✓ Reordering complete!\n');
fprintf('Check REORDERED_tube_numbering.png to verify.\n\n');

%% Save reordered masks
tubeMasks = tubeMasks_reordered;
centroids = centroids_reordered;

save('BMC_tubeMasks.mat', 'tubeMasks', 'tube_labels', 'phantomCenter', 'centroids', 'sortedIdx', 'tube_order');
fprintf('✓ Saved reordered masks to BMC_tubeMasks.mat\n');
fprintf('✓ Original tube_order mapping saved for reference\n\n');

%% Display mapping
fprintf('=== APPLIED MAPPING ===\n');
fprintf('Original detected → New position in shots.png\n');
for i = 1:24
    fprintf('  Detected tube %2d → Tube %2d (%s)\n', tube_order(i), i, tube_labels{i});
end

fprintf('\nYou can now continue with phantom_CEST_BMC_analysis.m from Step 7\n');
fprintf('Or rerun the entire script - it will load the reordered masks.\n');
