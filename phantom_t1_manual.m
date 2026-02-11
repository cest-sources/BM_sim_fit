clearvars; clc; close all;
tic
Starting_Directory=pwd;

%% Load and process the images
addpath('/Users/cbd/Downloads/CESTPhantomCode');
addpath('/Users/cbd/Downloads/MegaMouse_Code');
addpath('/Users/cbd/Downloads/Matfiles_New/ColorMaps');
cd(Starting_Directory);
disp('Select T1 scan');
DIR = uigetdir(pwd, 'Select T1 map containing dicom files which has .dcm');

% Get list of DICOM files in the folder
dicomFiles = dir(fullfile(DIR, '*.dcm'));
numFiles = length(dicomFiles);

% Read first file to preallocate
firstFile = fullfile(DIR, dicomFiles(1).name);
temp = dicomread(firstFile);
[xDim, yDim] = size(temp);
dicomVolume = zeros(xDim, yDim, numFiles);
mask_of_tubes = zeros(xDim, yDim);

TRs = zeros(numFiles,1);

%% Load all DICOM files and stack into 3D matrix
for i = 1:numFiles
    filePath = fullfile(DIR, dicomFiles(i).name);
    dicomVolume(:,:,i) = dicomread(filePath);
    info = dicominfo(filePath);

    % Extract TR (Repetition Time) in milliseconds
    if isfield(info, 'RepetitionTime')
        TRs(i) = double(info.RepetitionTime);
    else
        warning('No TR found in %s', dicomFiles(i).name);
        TRs(i) = NaN;
    end
end

%% Automatically detect the phantom outline instead of drawing manually
Sa = dicomVolume(:,:,10);  % First slice

% Step 1: Smooth the image to reduce noise
Sm = imgaussfilt(Sa, 4);

% Step 2: Normalize image to [0,1]
Smn = mat2gray(Sm);

% Step 3: Threshold (Otsu)
bw_phantom = imbinarize(Smn);

% Step 4: Keep only the largest connected component (phantom box)
CCp = bwconncomp(bw_phantom);
numPix = cellfun(@numel, CCp.PixelIdxList);
[~, idxMax] = max(numPix);
phantom_outline = false(size(bw_phantom));
phantom_outline(CCp.PixelIdxList{idxMax}) = true;

% Step 5: Fill any holes inside phantom
phantom_outline = imfill(phantom_outline, "holes");

% Optional: Smooth outline
phantom_outline = imopen(phantom_outline, strel('disk', 5));

% Show detected phantom outline
figure; 
imshowpair(Sa, phantom_outline, 'montage');
title('Auto-detected Phantom Outline (right)');

% Overlay outline in red on original
figure; 
imshow(Sa, [], 'InitialMagnification', 'fit'); colormap("gray"); hold on;
visboundaries(phantom_outline, 'Color','r','LineWidth',2);
title('Automatic Phantom Outline');

load("Manual_TubeMasks.mat");

% %%
% %T1 mapping
% 
dim = size(dicomVolume);        % [256 256 10]
numTubes = size(tubeMasks, 3); % 24

% Preallocate
T1map = zeros(dim(1), dim(2), numTubes);
M0    = zeros(dim(1), dim(2), numTubes);
T1_means = zeros(numTubes,1);

disp("Starting T1 mapping ...");

for t = 1:numTubes
    currentMask = tubeMasks(:,:,t) > 0;   % 2D mask for tube t
    
    dicomVolume_masked = dicomVolume;
    for i = 1:dim(3)
        dicomVolume_masked(:,:,i) = dicomVolume(:,:,i) .* currentMask;
    end

    for i = 1:dim(1)
        for j = 1:dim(2)
            if currentMask(i,j)
                % Signal across TRs for voxel (i,j)
                S = squeeze(dicomVolume_masked(i,j,:));  
                [M0(i,j,t), T1map(i,j,t),~,~] = t1fitting_VTR(S, TRs);
            end
        end
    end

    % --- Compute mean signal across the tube for plotting ---
    meanSignal = zeros(length(TRs),1);
    for k = 1:length(TRs)
        tmpImg = dicomVolume(:,:,k);
        meanSignal(k) = mean(tmpImg(currentMask), 'omitnan');  % mean over mask
    end

    % --- Fit mean signal for ROI using t1fitting_VTR ---
    [M0_fit, T1_fit, c_fit, S_fit] = t1fitting_VTR(meanSignal, TRs);

    subplot(4,6,t);
    scatter(TRs, S, 40, 'filled', 'MarkerFaceColor',[0.1 0.7 1]); hold on;
    plot(TRs, S_fit, 'r-', 'LineWidth', 1.5);
    
    xlabel('TR (ms)'); ylabel('Signal');
    title(sprintf('Tube %d: T1 = %.0f ms', t, T1_fit));
    legend('Data','Fit','Location','southeast');

    % Mean T1 only for current tube
    maskedT1 = T1map(:,:,t);
    T1_means(t) = mean(maskedT1(currentMask & maskedT1>0 & maskedT1<5000)); 
    % thresholds remove nonsense fits if needed
end

disp("Done!");

compositeT1 = zeros(size(tubeMasks(:,:,1)));

for t = 1:numTubes
    tubeMask = tubeMasks(:,:,t) > 0;
    T1tube   = T1map(:,:,t);
    compositeT1(tubeMask) = T1tube(tubeMask); % keep voxelwise T1 values
end

% Display
figure;
ax1 = axes;
imagesc(ax1, dicomVolume(:,:,10));  % anatomical background
axis image off;
colormap(ax1, 'gray');

ax2 = axes;
imagesc(ax2, compositeT1, 'AlphaData', compositeT1 > 0); % overlay T1 map
axis image off;
load('T1cm.mat');
colormap(ax2, T1colormap);
clim([0 8000])
ax2.Visible = 'off';

linkprop([ax1 ax2],'Position');
colorbar;
set(gca, 'FontSize', 18);
saveas(gcf, 'T1map_37.tiff');

T1_means_table = array2table(T1_means, ...
    'VariableNames', {'Mean_T1_ms'});
writetable(T1_means_table, 'T1_means.xlsx');

toc