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

%% Automatically find the tubes

St = dicomVolume(:,:,10);

bgd = imgaussfilt(St,4);     % Find the background
figure; imagesc(bgd); colormap("gray"); title("Background");
I = St - bgd;    % Subtract the background
figure; imagesc(I); colormap("gray"); title("Remove Background");

% threshold = mean(I, "all");
% bw = imbinarize (I, threshold);
bw = imbinarize (I);
figure; imagesc(bw); colormap("gray"); title("Binarize");
% Otsu, N., "A Threshold Selection Method from Gray-Level Histograms." IEEE Transactions on Systems, Man, and Cybernetics. Vol. 9, No. 1, 1979, pp. 62–66.

bw(~phantom_outline) = 0;
figure; imagesc(bw); colormap("gray"); title("Apply Phantom Box Outline");

%bw_temp = bw;

CC = bwconncomp(bw,8);
numOfPixels = cellfun(@numel,CC.PixelIdxList);
labeled = labelmatrix(CC);
RGB_label = label2rgb(labeled);
figure; imagesc(RGB_label); title("Connected Components - All");
%figure; imagesc(bw(CC.PixelIdxList)); colormap("gray"); title("Connected Components - All");
[~,indexOfMax] = max(numOfPixels);
bw(CC.PixelIdxList{indexOfMax}) = 0;
[indexOfSmall] = find(numOfPixels<10);
for i = 1:length(indexOfSmall)
    bw(CC.PixelIdxList{indexOfSmall(i)}) = 0;
end
figure; imagesc(bw); colormap("gray"); title("Connected Components - Tubes");

bw = imfill(bw,"holes");     % Fill holes
figure; imagesc(bw); colormap("gray"); title("Fill Holes");

bw = imclearborder(bw);

mask_of_tubes_bw = bw;

% Create a labeled image (each tube gets a unique ID)
labeledTubes = labelmatrix(CC);

% Convert to RGB for visualization
RGB_label = label2rgb(labeledTubes, 'jet', 'k', 'shuffle');

figure;
imagesc(RGB_label); title('Labeled Tubes');

%% Show index numbers on top of tubes
stats = regionprops(CC, 'Centroid');

% Compute phantom center (center of mass of outline)
phantomStats = regionprops(phantom_outline, 'Centroid');
phantomCenter = phantomStats.Centroid;

% Get centroids of all tubes
centroids = cat(1, stats.Centroid);

% Compute angle of each centroid relative to phantom center
angles = atan2(centroids(:,2) - phantomCenter(2), centroids(:,1) - phantomCenter(1));

% Sort tubes by angle (clockwise order)
[~, sortedIdx] = sort(angles);

%% Assign tube indices according to circular order
tubeIndex = struct('CCindex', [], 'CircularIndex', [], 'Centroid', []);
for i = 1:numel(sortedIdx)
    tubeIndex(i).CCindex       = sortedIdx(i);             % original CC index
    tubeIndex(i).CircularIndex = i;                        % new circular order
    tubeIndex(i).Centroid      = centroids(sortedIdx(i),:);
end

%% Create mask labeled with circular indices
circularMask = zeros(size(bw));
for i = 1:numel(sortedIdx)
    circularMask(CC.PixelIdxList{sortedIdx(i)}) = i; % assign circular index
end

%% Plot: Labeled tubes with circular indices
figure; imagesc(dicomVolume(:,:,10)); colormap("gray"); axis image off; 
title('Tube Indexing (Circular Order)');
hold on;
for i = 1:numel(tubeIndex)
    c = tubeIndex(i).Centroid;
    text(c(1), c(2), sprintf('%d', tubeIndex(i).CircularIndex), ...
        'Color','red','FontSize',12,'FontWeight','bold');
end
hold off;
set(gca, 'FontSize', 18);
saveas(gcf, 'TubeIndexing.tiff');

%% Show mask with circular order
figure; imagesc(circularMask); axis image off; colormap("jet");
title("Mask Labeled by Circular Index");

% % Reorder centroids according to circular order
% sortedCentroids = centroids(sortedIdx,:);
% 
% figure; imagesc(niiVolume(:,:,1)); colormap("gray"); axis image off; title('Tube Indexing');
% hold on;
% for k = 1:numel(stats)
%     c = stats(k).Centroid;
%     text(c(1), c(2), sprintf('%d', k), ...
%         'Color','red','FontSize',12,'FontWeight','bold');
% end
% hold off;
% 
% figure;
% ax1 = axes; 
% imagesc(ax1, niiVolume(:,:,1)); axis image off;
% colormap(ax1, "gray");
% ax2 = axes;
% imagesc(ax2, mask_of_tubes_bw, 'alphadata', mask_of_tubes_bw==1); axis image off;
% colormap(ax2, "jet");
% ax2.Visible = 'off';
% linkaxes([ax1,ax2]);
% title("Circle test");
% 
% CC2 = bwconncomp(bw,4);
% bw(CC2.PixelIdxList{1}) = 6;
% figure; imagesc(bw); colormap("gray"); title("Test Index");

    % I = niiVolume(:,:,1);
    % BW = edge(I,'Canny');
    %BW = edge(I,'Sobel');
figure(3);
imagesc(bw); colormap gray; axis image; hold on;

% Label connected components (each tube gets a unique ID)
labeledImage = bwlabel(bw);

% Measure properties (centroid, area, etc.)
tubeProps = regionprops(labeledImage, 'Centroid');

% Preallocate masks: 3D logical array (rows x cols x 24 tubes)
numTubes = 24;
tubeMasks = false([size(bw), numTubes]);

% Loop through all detected tubes
for i = 1:numTubes
    % Make a mask for tube i
    tubeMasks(:,:,i) = (labeledImage == i);

    % Mark centroid with tube number
    % c = tubeProps(i).Centroid;
    % text(c(1), c(2), sprintf('%d', i), ...
    %     'Color','red','FontSize',12,'FontWeight','bold');
end

figure(4);
imagesc(bw); colormap gray; axis image; hold on;

% Loop through all detected tubes
for i = 1:numTubes
    % Make a mask for tube i
    tubeMasks(:,:,i) = (labeledImage == i);

    %Mark centroid with tube number
    c = tubeProps(i).Centroid;
    text(c(1), c(2), sprintf('%d', i), ...
        'Color','red','FontSize',12,'FontWeight','bold');
end
set(gca, 'FontSize', 18);
saveas(gcf, 'TubeIndexing_2.tiff');
save('tubeMasks.mat', 'tubeMasks');

%%
%T1 mapping

dim = size(dicomVolume);        % [256 256 10]
numTubes = size(tubeMasks, 3); % 24

% Preallocate
T1map = zeros(dim(1), dim(2), numTubes);
M0    = zeros(dim(1), dim(2), numTubes);
T1_means = zeros(numTubes,1);

disp("Starting T1 mapping ...");

for t = 1:numTubes
    currentMask = tubeMasks(:,:,t) > 0;  % 2D mask for tube t
    
    dicomVolume_masked = dicomVolume;
    for i = 1:dim(3)
        dicomVolume_masked(:,:,i) = dicomVolume(:,:,i) .* currentMask;
    end
    
    for i = 1:dim(1)
        for j = 1:dim(2)
            if currentMask(i,j)
                % Signal across TRs for voxel (i,j)
                S = squeeze(dicomVolume_masked(i,j,:));  
                [M0(i,j,t), T1map(i,j,t)] = t1fitting_VTR(S, TRs);
            end
        end
    end

    % --- Compute mean signal across the tube for plotting ---
    meanSignal = zeros(length(TRs),1);
    for k = 1:length(TRs)
        tmpImg = dicomVolume_masked(:,:,k);
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
clim([0 5000]);
ax2.Visible = 'off';

linkprop([ax1 ax2],'Position');
colorbar;
set(gca, 'FontSize', 18);
saveas(gcf, 'T1map_37.tiff');

T1_means_table = array2table(T1_means, ...
    'VariableNames', {'Mean_T1_ms'});
writetable(T1_means_table, 'T1_means.xlsx');

toc